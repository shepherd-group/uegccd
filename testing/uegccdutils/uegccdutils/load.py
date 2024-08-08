import os.path as path
import numpy as np
import math

class Dataset():
	def __init__(self, dir=""):
		"""
		Load the output files associated 
		"""
		self.dir = path.expanduser(dir)
		self._load_Output()
		self._load_fort59()

	def _load_Output(self):

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

	def _load_fort59(self):

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
		pass

	def test_summary(self):
		"""
		Return a dictionary with summarizing information for testing purposes.

		Values stored in this dictionary should reflect the breadth of uegCCD's capabilities.
		This dictionary is used with pytest's pytest-regression's extension,
		which does not support Python objects like lists and arrays or NumPy types.
		Note that you must use the ".item()" method on NumPy objects to convert
		them to native Python types.
		"""

		i_max_change = np.argmax(self.summary["biggest changes"])

		results = {
			# from Output
			'summary': {
				'Final CCD Energy': self.summary["final CCD energy"],
				'CCD Correlation Energy': self.summary["final correlation energy"],
				'Total Iterations': self.summary["total iterations"],
				'Max CCD Residual': self.summary["max CCD residual"],
				'Biggest CCD Energy Change': self.summary["biggest changes"][i_max_change].item(),
				'Iteration of Biggest CCD Energy Change': self.summary["iterations"][i_max_change].item(),
			},
			'eigenvalues': {
				'Initial Iteration': {
					'Number of Eigenvalues': len(self.eigenvalues[0]['eigenvalues']),
					'SCF Energy': self.eigenvalues[0]['E(SCF)'],
					'MP2 Energy': self.eigenvalues[0]['E(2)']
				},
				'Final Iteration': {
					'Number of Eigenvalues': len(self.eigenvalues[-1]['eigenvalues']),
					'SCF Energy': self.eigenvalues[-1]['E(SCF)'],
					'MP2 Energy': self.eigenvalues[-1]['E(2)']
				}
			},

			# from fort.59
			'twist': {
				'Special Twist Angle': {
					'x': self.twist["special twist angle"][0].item(),
					'y': self.twist["special twist angle"][1].item(),
					'z': self.twist["special twist angle"][2].item(),
				},
				'Lowest Connectivity Diff Squared': self.twist["lowest connectivity diff squared"],
				'Averages': {
					'ntwist x': self.twist["averages"]["ntwist"][0].item(),
					'ntwist y': self.twist["averages"]["ntwist"][1].item(),
					'ntwist z': self.twist["averages"]["ntwist"][2].item(),
					'HF': self.twist["averages"]["HF"],
					'MP2': self.twist["averages"]["MP2"],
					'CCD': self.twist["averages"]["CCD"],
				},
				'Errors': {
					'ntwist x': self.twist["errors"]["ntwist"][0].item(),
					'ntwist y': self.twist["errors"]["ntwist"][1].item(),
					'ntwist z': self.twist["errors"]["ntwist"][2].item(),
					'HF': self.twist["errors"]["HF"],
					'MP2': self.twist["errors"]["MP2"],
					'CCD': self.twist["errors"]["CCD"],
				},
				'Eigenvalue Averages': {
					'min': np.min(self.twist["eigenvalue averages"]).item(),
					'max': np.max(self.twist["eigenvalue averages"]).item(),
				},
				'Connectivity Averages': {
					'min': np.min(self.twist["connectivity averages"]).item(),
					'max': np.max(self.twist["connectivity averages"]).item(),
				},
			},
		}

		return results
