import os.path as path
import numpy as np
import math

class Dataset():
	def __init__(self, dir):
		self.dir = path.expanduser(dir)
		self._load_Output()
		self._load_fort59()
		self._load_fort60()
		self._load_fort61()

	def _load_Output(self):

		file = path.join(self.dir, "Output")

		summary = {
			"CCD Energy": [],
			"Iteration": [],
			"Biggest Change": []
			}
		blocks = []
		cur_blc = -1

		with open(file) as f:
			
			summary_block = False

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
					summary["Total Iterations"] = int(l.split()[4])
				elif "Ec(CCD) is" in l:
					summary["Ec(CCD)"] = float(l.split()[2])
				elif "Final CCD Energy is" in l:
					summary["Final CCD Energy"] = float(l.split()[4])
				elif "Max CCD Residual is" in l:
					summary["Max CCD Residual"] = float(l.split()[4])
				
				# Line for identifying eigenvalue blocks
				elif "E(SCF) is" in l: 
					# start of new eigenvalue block
					cur_blc += 1
					postNOcc = False
					blocks.append({
						"E(SCF)" : float(l.split()[2]),
						"Eig Pre-NOcc": [],
						"Eig Post-NOcc": []
						})
				elif "E(SCF) + E(2) is" in l:
					# end of block
					assert math.isclose(float(l.split()[4]),
					 blocks[cur_blc]["E(SCF)"] + blocks[cur_blc]["E(2)"]), f"Error reading entry {cur_blc} in file {file}."
					blocks[cur_blc]["Eig Pre-NOcc"] = np.array(blocks[cur_blc]["Eig Pre-NOcc"])
					blocks[cur_blc]["Eig Post-NOcc"] = np.array(blocks[cur_blc]["Eig Post-NOcc"])
				elif "E(2) is" in l:
					blocks[cur_blc]["E(2)"] = float(l.split()[2])
				elif l.strip() == '-'*24:
					# Switch from NOcc to rest of eigenvalues
					postNOcc = True

				else:
					# Process tabular data
					data = l.split()

					if len(data) == 1:
						# Eigenvalue
						eig = float(data[0])
						if postNOcc:
							blocks[cur_blc]["Eig Post-NOcc"].append(eig)
						else:
							blocks[cur_blc]["Eig Pre-NOcc"].append(eig)

					elif len(data) == 3:
						eng = float(data[0])
						itr = float(data[1])
						chg = float(data[2])
						summary["CCD Energy"].append(eng)
						summary["Iteration"].append(itr)
						summary["Biggest Change"].append(chg)

					else:
						raise RuntimeError(f"File {file} contains unrecognized line:\n{l}")
		
		summary["CCD Energy"] = np.array(summary["CCD Energy"])
		summary["Iteration"] = np.array(summary["Iteration"])
		summary["Biggest Change"] = np.array(summary["Biggest Change"])

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
						"hf": float(data[-5]),
						"mp2": float(data[-4]),
						"ccd": float(data[-3])
					}

				elif "eigenvalue averages" in line:
					eig_avg = np.array([float(s) for s in data[:-2]])

				elif "averages" in line:
					# ntwist seems to always have 3 entries,
					# but let's be flexible by counting backwards
					avg = {
						"ntwist": np.array([float(s) for s in data[:-4]]),
						"hf": float(data[-4]),
						"mp2": float(data[-3]),
						"ccd": float(data[-2])
					}

				elif "errors" in line:
					# ntwist seems to always have 3 entries,
					# but let's be flexible by counting backwards
					err = {
						"ntwist": np.array([float(s) for s in data[:-4]]),
						"hf": float(data[-4]),
						"mp2": float(data[-3]),
						"ccd": float(data[-2])
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

		self.twist = {
			"ntwist": np.array(ntwist),
			"hf": np.array(hf),
			"mp2": np.array(mp2),
			"ccd": np.array(ccd),
			"averages": avg,
			"squared averages": sq_avg,
			"errors": err,
			"eigenvalue averages": eig_avg,
			"connectivity average": conn,
			"itwist": np.array(itwist),
			"connectivity diff squared": np.array(conn_diff),
			"special twist angle": special,
			"lowest connectivity diff squared": low_conn_diff,
		}

	def _load_fort60(self):
		pass

	def _load_fort61(self):
		pass

if __name__ == "__main__":
	dir = "N_00014_0002/"
	ds = Dataset(dir)