import os.path as path
import numpy as np
import math

class Dataset():
    def __init__(self, dir="."):
        """
        Load the output files stored inside a given uegCCD simulation directory.
        By default, load the current directory.
        """
        self.dir = path.expanduser(dir)
        self._load_Output()
        self._load_fort58()
        self._load_fort59()
        self._load_fort80()
        self._load_fort81()
        self._load_fort82()

    def _load_Output(self):
        """
        Load the primary output file.
        """

        file = path.join(self.dir, "Output")

        summary = {
            "correlation energies": [],
            "iterations": [],
            "biggest changes": []
            }
        blocks = []
        cur_blc = -1

        with open(file) as f:
            for line in f:

                l = line.strip(' \n*')

                if not l: 
                    # divider
                    continue
                elif l == '-'*48:
                    # divider
                    continue
                elif "RHF" in l:
                    # header
                    continue
                elif "CCD summary follows" in l:
                    # header
                    continue
                elif "Including" in l:
                    # info line
                    continue
                elif "Biggest change" in l:
                    # header
                    continue

                # Lines with summary info
                elif "CCD has converged" in l:
                    summary["total iterations"] = int(l.split()[4])
                elif "Ec(CCD) is" in l:
                    # correlation energy; primary output
                    summary["final correlation energy"] = float(l.split()[2])
                elif "Final CCD Energy is" in l:
                    summary["final CCD energy"] = float(l.split()[4])
                elif "Max CCD Residual is" in l:
                    summary["max CCD residual"] = float(l.split()[4])
                
                # Line in eigenvalue blocks
                elif "E(SCF) is" in l: 
                    # start of new eigenvalue block
                    cur_blc += 1
                    blocks.append({
                        "E(SCF)" : float(l.split()[2]),
                        "eigenvalues": [],
                        })
                elif "E(SCF) + E(2) is" in l:
                    # end of block
                    assert math.isclose(float(l.split()[4]),
                     blocks[cur_blc]["E(SCF)"] + blocks[cur_blc]["E(2)"]), f"Error reading entry {cur_blc} in file {file}."
                    blocks[cur_blc]["eigenvalues"] = np.array(blocks[cur_blc]["eigenvalues"])
                elif "E(2) is" in l:
                    blocks[cur_blc]["E(2)"] = float(l.split()[2])
                elif l.strip() == '-'*24:
                    # Ignore
                    continue

                # Lines with data tables
                else:
                    # Process tabular data
                    data = l.split()

                    if len(data) == 1:
                        # Eigenvalue
                        eig = float(data[0])
                        blocks[cur_blc]["eigenvalues"].append(eig)
                        
                    elif len(data) == 3:
                        eng = float(data[0])
                        itr = float(data[1])
                        chg = float(data[2])
                        summary["correlation energies"].append(eng)
                        summary["iterations"].append(itr)
                        summary["biggest changes"].append(chg)

                    else:
                        raise RuntimeError(f"File {file} contains unrecognized line:\n{l}")
        
        # Final data structure cleanup
        summary["correlation energies"] = np.array(summary["correlation energies"])
        summary["iterations"] = np.array(summary["iterations"])
        summary["biggest changes"] = np.array(summary["biggest changes"])

        self.eigenvalues = blocks
        self.summary = summary

    def _load_fort58(self):
        """
        Load the output file containing connectivity data.
        """

        file = path.join(self.dir, "fort.58")

        conn = np.genfromtxt(file, dtype=int)

        self.connectivity = conn

    def _load_fort59(self):
        """
        Load the output file containing twist angle data.
        """

        file = path.join(self.dir, "fort.59")

        ntwist = []
        hf = []
        mp2 = []
        ccd = []
        itwist = []
        conn_diff = []

        with open(file) as f:
            for line in f:

                data = line.split()

                if "special twist angle found" in line:
                    # Angle may not always be 3D, so be flexible
                    special = np.array([float(s) for s in data[4:]])

                elif "lowest connectivity diff squared" in line:
                    low_conn_diff = float(data[-1])

                elif "squared averages" in line:
                    # ntwist seems to always have 3 entries,
                    # but let's be flexible by counting backwards
                    sq_avg = {
                        "ntwist": np.array([float(s) for s in data[:-5]]),
                        "HF": float(data[-5]),
                        "MP2": float(data[-4]),
                        "CCD": float(data[-3])
                    }

                elif "eigenvalue averages" in line:
                    eig_avg = np.array([float(s) for s in data[:-2]])

                elif "averages" in line:
                    # ntwist seems to always have 3 entries,
                    # but let's be flexible by counting backwards
                    avg = {
                        "ntwist": np.array([float(s) for s in data[:-4]]),
                        "HF": float(data[-4]),
                        "MP2": float(data[-3]),
                        "CCD": float(data[-2])
                    }

                elif "errors" in line:
                    # ntwist seems to always have 3 entries,
                    # but let's be flexible by counting backwards
                    err = {
                        "ntwist": np.array([float(s) for s in data[:-4]]),
                        "HF": float(data[-4]),
                        "MP2": float(data[-3]),
                        "CCD": float(data[-2])
                    }

                elif "conn" in line:
                    conn = np.array([float(s) for s in data[:-1]])

                else:
                    
                    if len(data) >= 4:
                        # ntwist seems to always have 3 entries,
                        # but let's be flexible by counting backwards
                        ntwist.append([float(s) for s in data[:-3]])
                        hf.append(float(data[-3]))
                        mp2.append(float(data[-2]))
                        ccd.append(float(data[-1]))
                    
                    elif len(data) == 2:
                        # iTwist, ConnectivityDiff2
                        itwist.append(int(data[0]))
                        conn_diff.append(float(data[1]))

                    else:
                        raise RuntimeError(f"File {file} contains unrecognized line:\n{line}")

        self.twist = {
            "ntwist": np.array(ntwist),
            "HF": np.array(hf),
            "MP2": np.array(mp2),
            "CCD": np.array(ccd),
            "averages": avg,
            "squared averages": sq_avg,
            "errors": err,
            "eigenvalue averages": eig_avg,
            "connectivity averages": conn,
            "itwist": np.array(itwist),
            "connectivity diff squared": np.array(conn_diff),
            "special twist angle": special,
            "lowest connectivity diff squared": low_conn_diff,
        }					

    def _load_fort80(self):
        """
        Load the output file containing transition structure factor data for the CC calculation.
        """
        
        file = path.join(self.dir, "fort.80")

        sfactors = []
        curr_fctr = -1

        with open(file) as f:
            for line in f:
                
                data = line.split()

                if line.strip() == "--------":
                    continue
                if "Starting a new structure factor" in line:
                    curr_fctr += 1
                    sfactors.append({
                        "integerized momentum transfer vector": [],
                        "momentum transfer vector magnitude": [],
                        "S(G)": [],
                        "Coulomb potential": [],
                        "non-exchange S(G)": [],
                        "exchange S(G)": []
                    })
                else:
                    if len(data) == 8:
                        sfactors[curr_fctr]["integerized momentum transfer vector"].append(
                            [int(i) for i in data[:3]]
                        )
                        sfactors[curr_fctr]["momentum transfer vector magnitude"].append(float(data[3]))
                        sfactors[curr_fctr]["S(G)"].append(float(data[4]))
                        sfactors[curr_fctr]["Coulomb potential"].append(float(data[5]))
                        sfactors[curr_fctr]["non-exchange S(G)"].append(float(data[6]))
                        sfactors[curr_fctr]["exchange S(G)"].append(float(data[7]))
                    else:
                        raise RuntimeError(f"File {file} contains unrecognized line:\n{line}")

        # Convert data lists to NumPy arrays		
        for fctr in sfactors:
            for quant, val in fctr.items():
                fctr[quant] = np.array(val)

        self.CCD_structure_factors = sfactors		

    def _load_fort81(self):
        """
        Load the output file containing transition structure factor data for the MP2 calculation.
        """
        
        file = path.join(self.dir, "fort.81")

        sfactors = []
        curr_fctr = -1

        with open(file) as f:
            for line in f:
                
                data = line.split()

                if line.strip() == "--------":
                    continue
                if "Starting a new structure factor" in line:
                    curr_fctr += 1
                    sfactors.append({
                        "integerized momentum transfer vector": [],
                        "momentum transfer vector magnitude": [],
                        "S(G)": [],
                        "Coulomb potential": [],
                    })
                else:
                    if len(data) == 6:
                        sfactors[curr_fctr]["integerized momentum transfer vector"].append(
                            [int(i) for i in data[:3]]
                        )
                        sfactors[curr_fctr]["momentum transfer vector magnitude"].append(float(data[3]))
                        sfactors[curr_fctr]["S(G)"].append(float(data[4]))
                        sfactors[curr_fctr]["Coulomb potential"].append(float(data[5]))
                    else:
                        raise RuntimeError(f"File {file} contains unrecognized line:\n{line}")

        # Convert data lists to NumPy arrays		
        for fctr in sfactors:
            for quant, val in fctr.items():
                fctr[quant] = np.array(val)

        self.MP2_structure_factors = sfactors		

    def _load_fort82(self):
        """
        Load the output file containing transition structure factor data for the HF calculation.

        The names for quantities in the resulting HF_structure_factors object are extrememly questionable - CK
        """
        
        file = path.join(self.dir, "fort.82")

        sfactors = []
        curr_fctr = -1

        with open(file) as f:
            for line in f:
                
                data = line.split()

                if line.strip() == "--------":
                    continue
                if "Starting a new structure factor" in line:
                    curr_fctr += 1
                    sfactors.append({
                        "integerized momentum transfer vector": [],
                        "structure factor exchange energy": [],
                        "Coulomb potential exchange energy": [],
                    })
                else:
                    if len(data) == 5:
                        sfactors[curr_fctr]["integerized momentum transfer vector"].append(
                            [int(i) for i in data[:3]]
                        )
                        sfactors[curr_fctr]["structure factor exchange energy"].append(float(data[3]))
                        sfactors[curr_fctr]["Coulomb potential exchange energy"].append(float(data[4]))
                    else:
                        raise RuntimeError(f"File {file} contains unrecognized line:\n{line}")

        # Convert data lists to NumPy arrays		
        for fctr in sfactors:
            for quant, val in fctr.items():
                fctr[quant] = np.array(val)

        self.HF_structure_factors = sfactors	

    def test_summary(self, tol=1e-13):
        """
        Return a dictionary with summarizing information for regression testing purposes.
        Precision of output is controlled by "tol" parameter.

        Values stored in this dictionary should reflect the breadth of uegCCD's capabilities.
        This dictionary is used with pytest's pytest-regression's extension,
        which does not support Python objects like lists and arrays or NumPy types.
        Note that you must use the ".item()" method on NumPy objects to convert
        them to native Python types.
        """

        i_max_change = np.argmax(self.summary["biggest changes"])
        ndigits = int(-np.log10(tol))

        results = {
            # from Output
            'Summary': {
                'Final CCD Energy': round(self.summary["final CCD energy"], ndigits),
                'CCD Correlation Energy': round(self.summary["final correlation energy"], ndigits),
                'Total Iterations': self.summary["total iterations"],
                'Max CCD Residual': round(self.summary["max CCD residual"], ndigits),
                'Biggest CCD Energy Change': round(self.summary["biggest changes"][i_max_change].item(), ndigits),
                'Iteration of Biggest CCD Energy Change': int(self.summary["iterations"][i_max_change].item()),
            },
            'Eigenvalues': {
                'Initial Iteration': {
                    'Number of Eigenvalues': len(self.eigenvalues[0]['eigenvalues']),
                    'SCF Energy': round(self.eigenvalues[0]['E(SCF)'], ndigits),
                    'MP2 Energy': round(self.eigenvalues[0]['E(2)'], ndigits),
                },
                'Final Iteration': {
                    'Number of Eigenvalues': len(self.eigenvalues[-1]['eigenvalues']),
                    'SCF Energy': round(self.eigenvalues[-1]['E(SCF)'], ndigits),
                    'MP2 Energy': round(self.eigenvalues[-1]['E(2)'], ndigits),
                }
            },

            # from fort.58
            'Connectivity':{
                "Matrix Dim 0": self.connectivity.shape[0],
                "Matrix Dim 1": self.connectivity.shape[1],
                "Maximum": self.connectivity.max().item(),
                "Number of Unique Entries": np.unique(self.connectivity).size
            },

            # from fort.59
            'Twist Angles': {
                'Special Twist Angle': {
                    'x': round(self.twist["special twist angle"][0].item(), ndigits),
                    'y': round(self.twist["special twist angle"][1].item(), ndigits),
                    'z': round(self.twist["special twist angle"][2].item(), ndigits),
                },
                'Lowest Connectivity Diff Squared': self.twist["lowest connectivity diff squared"],
                'Averages': {
                    'ntwist x': round(self.twist["averages"]["ntwist"][0].item(), ndigits),
                    'ntwist y': round(self.twist["averages"]["ntwist"][1].item(), ndigits),
                    'ntwist z': round(self.twist["averages"]["ntwist"][2].item(), ndigits),
                    'HF': round(self.twist["averages"]["HF"], ndigits),
                    'MP2': round(self.twist["averages"]["MP2"], ndigits),
                    'CCD': round(self.twist["averages"]["CCD"], ndigits),
                },
                'Errors': {
                    'ntwist x': round(self.twist["errors"]["ntwist"][0].item(), ndigits),
                    'ntwist y': round(self.twist["errors"]["ntwist"][1].item(), ndigits),
                    'ntwist z': round(self.twist["errors"]["ntwist"][2].item(), ndigits),
                    'HF': round(self.twist["errors"]["HF"], ndigits),
                    'MP2': round(self.twist["errors"]["MP2"], ndigits),
                    'CCD': round(self.twist["errors"]["CCD"], ndigits),
                },
                'Eigenvalue Averages': {
                    'min': round(np.min(self.twist["eigenvalue averages"]).item(), ndigits),
                    'max': round(np.max(self.twist["eigenvalue averages"]).item(), ndigits),
                },
                'Connectivity Averages': {
                    'min': round(np.min(self.twist["connectivity averages"]).item(), ndigits),
                    'max': round(np.max(self.twist["connectivity averages"]).item(), ndigits),
                },
            },

            # from fort.80
            'CCD Structure Factors': {
                'Initial Structure Factor': {
                    'Sum_G( S(G)*V(G) )': round(np.dot(
                        self.CCD_structure_factors[0]["S(G)"],
                        self.CCD_structure_factors[0]["Coulomb potential"]
                        ).item(), ndigits),
                    'Momentum Transfer Vector Magnitude': {
                        'min': round(np.min(
                            self.CCD_structure_factors[0]["momentum transfer vector magnitude"]
                            ).item(), ndigits),
                        'max': round(np.max(
                            self.CCD_structure_factors[0]["momentum transfer vector magnitude"]
                            ).item(), ndigits),
                    }
                },
                'Final Structure Factor': {
                    'Sum_G( S(G)*V(G) )': round(np.dot(
                        self.CCD_structure_factors[-1]["S(G)"],
                        self.CCD_structure_factors[-1]["Coulomb potential"]
                        ).item(), ndigits),
                    'Momentum Transfer Vector Magnitude': {
                        'min': round(np.min(
                            self.CCD_structure_factors[-1]["momentum transfer vector magnitude"]
                            ).item(), ndigits),
                        'max': round(np.max(
                            self.CCD_structure_factors[-1]["momentum transfer vector magnitude"]
                            ).item(), ndigits),
                    }	
                }
            },

            # from fort.81
            'MP2 Structure Factors': {
                'Initial Structure Factor': {
                    'Sum_G( S(G)*V(G) )': round(np.dot(
                        self.MP2_structure_factors[0]["S(G)"],
                        self.MP2_structure_factors[0]["Coulomb potential"]
                        ).item(), ndigits),
                    'Momentum Transfer Vector Magnitude': {
                        'min': round(np.min(
                            self.MP2_structure_factors[0]["momentum transfer vector magnitude"]
                            ).item(), ndigits),
                        'max': round(np.max(
                            self.MP2_structure_factors[0]["momentum transfer vector magnitude"]
                            ).item(), ndigits),
                    }
                },
                'Final Structure Factor': {
                    'Sum_G( S(G)*V(G) )': round(np.dot(
                        self.MP2_structure_factors[-1]["S(G)"],
                        self.MP2_structure_factors[-1]["Coulomb potential"]
                        ).item(), ndigits),
                    'Momentum Transfer Vector Magnitude': {
                        'min': round(np.min(
                            self.MP2_structure_factors[-1]["momentum transfer vector magnitude"]
                            ).item(), ndigits),
                        'max': round(np.max(
                            self.MP2_structure_factors[-1]["momentum transfer vector magnitude"]
                            ).item(), ndigits),
                    }	
                }
            },

            # from fort.82
            'HF Structure Factors': {
                'Initial Structure Factor': {
                    'Sum( ExSf*ExV )': round(np.dot(
                        self.HF_structure_factors[0]["structure factor exchange energy"],
                        self.HF_structure_factors[0]["Coulomb potential exchange energy"]
                    ).item(), ndigits),
                },
                'Final Structure Factor': {
                    'Sum( ExSf*ExV )': round(np.dot(
                        self.HF_structure_factors[-1]["structure factor exchange energy"],
                        self.HF_structure_factors[-1]["Coulomb potential exchange energy"]
                    ).item(), ndigits),
                }
            }
        }

        return results
