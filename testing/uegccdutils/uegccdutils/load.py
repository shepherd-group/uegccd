import os.path as path
import numpy as np
import math
import pytest

class Dataset():
	def __init__(self, dir=""):
		"""
		Load the output files associated 
		"""
		self.dir = path.expanduser(dir)
		self._load_Output()
		self._load_fort59()
		self._load_fort60()
		self._load_fort61()

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
					postNOcc = False
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

	def _load_fort60(self):
		
		file = path.join(self.dir, "fort.60")

		divider = set('-')
		basis = {
			"i": [],
			"culm": [],
			"culm*2": [],
		}
		init_pass = True

		with open(file) as f:
			for line in f:

				l = line.strip()

				if not l: 
					# empty line
					continue
				elif set(l) == divider:
					# only dashes; a divider
					# use to reset block flags
					basis_block = False
				
				elif "Basis sets available" in l:
					# header
					continue
				elif "Calculating Madelung Constant" in l:
					# header
					continue

				# Lines in a basis set block
				elif "# ecut k-points M" in l:
					# use as start of basis set block
					# so that we know block ends on next divider
					basis_block = True

				# Lines in a Madelung constant block
				elif "kappa taken from" in l:
					# use a start of Madelung block
					mdlg_blc = {
						"kappa": float(l.split()[-1]),
						"reciprocal space iterations": [],
						"real space iterations": []
					}
				elif "term2" in l:
					mdlg_blc["term2"] = float(l.split()[0])
				elif "term4" in l:
					# this line also signifies start of iterations
					# for the reciprocal space term
					mdlg_blc["term4"] = float(l.split()[0])
					reciprocal_block = True
				elif "reciprocal space" in l:
					# clean up iteration data structure
					block_data = np.array(mdlg_blc["reciprocal space iterations"])
					mdlg_blc["reciprocal space iterations"] = block_data
					
					# check final reciprocal space value
					final = float(l.split()[-1])
					assert math.isclose(final, block_data[-1,-1])
					mdlg_blc["reciprocal space"] = final # save separately for easier access
					
					# switch block type
					reciprocal_block = False
					real_block = True
				elif "real space" in l:
					# clean up iteration data structure
					block_data = np.array(mdlg_blc["real space iterations"])
					mdlg_blc["real space iterations"] = block_data

					# check final real space value
					final = float(l.split()[-1])
					assert math.isclose(final, block_data[-1,-1])
					mdlg_blc["real space"] = final # save separately for easier access

					# end block
					real_block = False
				elif "Madelung constant calculated as" in l:
					# Finish with the initial pass
					init_pass = False
					mdlg_blc["Madelung constant"] = l.split()[-1]
					init = mdlg_blc

				# Lines with data tables
				else:
					if basis_block:
						data = l.split()
						basis["i"].append(float(data[0]))
						basis["culm"].append(float(data[1]))
						basis["culm*2"].append(float(data[2]))

					elif reciprocal_block:
						data = l.split()
						data_flt = (float(data[0]), float(data[1]))
						mdlg_blc["reciprocal space iterations"].append(data_flt)

					elif real_block:
						data = l.split()
						data_flt = (float(data[0]), float(data[1]))
						mdlg_blc["real space iterations"].append(data_flt)

					else:
						pass # ignore the rest of the file for now

		# Final data structure cleanup
		# Madelung block data structures are cleaned up while parsing that block
		basis["i"] = np.array(basis["i"])
		basis["culm"] = np.array(basis["culm"])
		basis["culm*2"] = np.array(basis["culm*2"])
		self.heg = {
			"basis": basis,
			"initialization": init,
		}

	def _load_fort61(self):

		file = path.join(self.dir, "fort.61")
		with open(file) as f:
			for line in f:

				data = line.split()

				if "core" in line:
					core = int(data[0])
				elif "active" in line:
					active = int(data[0])
				else:
					raise RuntimeError(f"File {file} contains unrecognized line:\n{line}")

		self.io = {
			'core electrons': core,
			'active electrons': active}

	def test_summary(self):
		"""
		Return a dictionary with summarizing information for testing purposes.

		Values stored in this dictionary should reflect the breadth of uegCCD's capabilities.
		This dictionary is used with pytest's pytest-regression's extension,
		which does not support Python objects like lists and arrays.
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
			# 'eigenvalues': {
			# 	'First Iteration': {},
			# 	'Last Iteration': {},
			# },
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
