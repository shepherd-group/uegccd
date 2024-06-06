import os.path as path
import numpy as np
import math

class uegCCD_Dataset():
	def __init__(self, dir):
		self.dir = path.expanduser(dir)
		self.load_Output()

	def load_Output(self):

		summary = {
			"CCD Energy": [],
			"Iteration": [],
			"Biggest Change": []
			}
		blocks = []
		cur_blc = -1

		file = path.join(self.dir, "Output")

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
				elif "E(2)" in l and not "E(SCF)" in l:
					blocks[cur_blc]["E(2)"] = float(l.split()[2])
				elif "E(SCF) + E(2)" in l:
					# end of block
					assert math.isclose(float(l.split()[4]),
					 blocks[cur_blc]["E(SCF)"] + blocks[cur_blc]["E(2)"]), f"Error reading entry {cur_blc} in file {file}."
					blocks[cur_blc]["Eig Pre-NOcc"] = np.array(blocks[cur_blc]["Eig Pre-NOcc"])
					blocks[cur_blc]["Eig Post-NOcc"] = np.array(blocks[cur_blc]["Eig Post-NOcc"])
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

		self.eigenvalue_blocks = blocks
		self.summary = summary

if __name__ == "__main__":
	dir = "N_00014_0002/"
	ds = uegCCD_Dataset(dir)