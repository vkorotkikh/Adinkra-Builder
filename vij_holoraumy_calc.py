# ******************************************************************************
#
# Name:    Permutation Matrix Calculation
# Author:  Vadim Korotkikh
# Email:   Vadim.Korotkikh
# Date:    Aug 2016
# Version: 1.0
#
# Description:
#
# ******************************************************************************


# ******************************************************************************
# Begin Imports
import math
import sys
import numpy as np
import numpy.matlib
import itertools
from numpy import array
from numpy.linalg import inv
import time


# import matrix_outerprod_calc
import alpha_beta_4x4

# ******************************************************************************
# Main() function.
def main():
	# gen_signm(4)
	pass

# ******************************************************************************
# Do the final Vij calculation
def calculate_vij_matrices(main_tetrad_list):

	""" Remember that the main_tetrad_ark is a list of lists,
		with each list containing four tuples, with tuples being
		matrix number and the matrices itself. """

	vij_possibilities = []
	vij_possibilities = alpha_beta_4x4.illuminator_of_elfes()
	# vij_possibilities = [alpha1i, alpha2i, alpha3i, beta1i, beta2i, beta3i]
	vij_sixset 		= []

	print("							")
	print("	Calculating Vij matrices")
	print("							")
	vij_alphas 		= []
	vij_betas  		= []
	calc_check		= []

	vij_matrices	= []

	anomaly_switch = 0

	for ti, teti in enumerate(main_tetrad_list):
		print("# ********************************")
		print("								     ")
		print("Tetrad i: ", ti)
		# print(teti[0][1][0,:], teti[1][1][0,:], teti[2][1][0,:], teti[3][1][0,:])
		# print(teti[0][1][1,:], teti[1][1][1,:], teti[2][1][1,:], teti[3][1][1,:])
		# print(teti[0][1][2,:], teti[1][1][2,:], teti[2][1][2,:], teti[3][1][2,:])
		# print(teti[0][1][3,:], teti[1][1][3,:], teti[2][1][3,:], teti[3][1][3,:])
		print("								     ")

		temp_combos = []
		alpha_temp	= []
		beta_temp   = []
		vij_tempset = []

		""" Store 6 Vij matrices in temp_vijmat"""
		temp_vijmat		= []

		""" This section does a double loop over the same tetrad to calculate
		the set of 6 Vij matrices for the tetrad.
		So for each matrix in the tetrad its checked against all the possible others,
	 	bypassing the duplicate calculations
		"""
		for i, li in enumerate(teti):
			# print(li[1])
			bigli = li[1]
			tr_bigli = np.transpose(bigli)

			for j, lj in enumerate(teti):
				biglj = lj[1]
				ij_temp = [i, j]
				ij_temp.sort()
				ir = i + 1
				jr = j + 1
				ijstr = str(ir) + str(jr)
				if ij_temp not in temp_combos and i != j:
					# print("Vij matrix i-j vals:", ij_temp)
					print("Vij matrix i-j vals:", ijstr)
					temp_combos.append(ij_temp)
					tr_biglj = np.transpose(biglj)
					# temp_mat = np.dot(tr_bigli, biglj) - np.dot(tr_biglj, bigli)
					""" Vij eq from 1601.00 (3.2) """
					# temp_mat = np.matmul(tr_biglj, bigli) - np.matmul(tr_bigli, biglj)
					temp_mat = np.dot(tr_bigli, biglj) - np.dot(tr_biglj, bigli)
					""" Compare against the 6 possible matrix solutions """
					tf_bool = 0
					for xi, ijx in enumerate(vij_possibilities):
						ijx_neg = np.multiply(ijx, -1)
						# print(xi)
						if np.array_equal(temp_mat, ijx):
							tf_bool = 1
							temp_vijmat.append(temp_mat)
							print("*************$$$$$$$$$$$$$$$$$$ ")
							print("l-solution found:")
							print(ijx)
							tmint = np.int(1)
							if xi < 3:
								tmp_str = "alpha" + str((xi + 1))
								# print(tmp_str)
								# vij_temp.append((tmp_str, tmint))
								vij_tempset.append([tmp_str, ijstr, tmint])
								alpha_temp.append([tmp_str, ijstr, tmint])
							elif xi >= 3:
								tmp_str = "beta" + str((xi - 2))
								# vij_temp.append(xi + 1)
								# vij_temp.append((tmp_str, tmint))
								vij_tempset.append([tmp_str, ijstr, tmint])
								beta_temp.append([tmp_str, ijstr, tmint])
						elif np.array_equal(temp_mat, ijx_neg):
							tf_bool = 1
							temp_vijmat.append(temp_mat)
							print("*************$$$$$$$$$$$$$$$$$$ ")
							print("l-solution found:")
							print(ijx_neg)
							# xint = (xi + 1) * ( -1)
							# vij_temp.append(xint)
							tmint = np.int(-1)
							if xi < 3:
								# tmp_str = "alpha_" + str((xi + 1) * ( -1))
								tmp_str = "alpha" + str((xi + 1))
								# print(tmp_str)
								vij_tempset.append([tmp_str, ijstr, tmint])
								alpha_temp.append([tmp_str, ijstr, tmint])
							elif xi >= 3:
								# tmp_str = "beta_" + str((xi + 1) * ( -1))
								tmp_str = "beta" + str((xi - 2))
								# vij_temp.append(xi + 1)
								vij_tempset.append([tmp_str, ijstr, tmint])
								beta_temp.append([tmp_str, ijstr, tmint])
						else:
							if i != j and tf_bool == 0 and xi >= 5:
								if not(np.array_equal(temp_mat, ijx)) or not np.array_equal(temp_mat, ijx_neg):
									print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ")
									print("Anomaly found:",i,j)
									print(temp_mat)
									anomaly_switch = 1
					tf_bool = 0
							# print(i,j)
							# pass
		# print(vij_tempset)

		# if vij_tempset == checkset:
		# 	print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ")
		# 	print("P3 set Beta match found")
		# 	print("P3 set match=true")
		vij_matrices.append(temp_vijmat)
		calc_check.append(vij_tempset)


		if alpha_temp:
			vij_alphas.append(alpha_temp)
		elif beta_temp:
			vij_betas.append(beta_temp)
		beta_temp 	= []
		alpha_temp 	= []

	print("*************$$$$$$$$$$$$$$$$$$ ")
	print("Vij Matrix Coefficients Results:")
	print("")
	for mvals in calc_check:
		if any(x for x in mvals if x[0].startswith('alpha')) and any(x for x in mvals if x[0].startswith('beta')):
			print("MIXED ALPHA_BETA ERROR")
			print(mvals)
		else:
			print(mvals)

	# for i, iv in enumerate(vij_matrices):
	#
	# 	for j, jv in enumerate(vij_matrices):
	#
	# 		# if i != j:
	# 		if i != j and i < 6 and j < 6:
	#
	# 			cl1 = [np.multiply(iv[0], mat) for mat in jv]
	# 			cl2 = [np.multiply(iv[1], mat) for mat in jv]
	# 			cl3 = [np.multiply(iv[2], mat) for mat in jv]
	# 			cl4 = [np.multiply(iv[3], mat) for mat in jv]
	# 			cl5 = [np.multiply(iv[4], mat) for mat in jv]
	# 			cl6 = [np.multiply(iv[5], mat) for mat in jv]
	#
	# 			# cl1sum = np
	# 			tempsum = [jv
	# 			for mat in iv
	#
	# 		else:
	# 			pass
	print(len(vij_alphas))
	print(len(vij_betas))

	gadget_vals		= []
	one_count 		= 0
	ptre_count		= 0
	ntre_count		= 0
	zero_count		= 0
	if not anomaly_switch:
		for fi, ijf in enumerate(calc_check):

			for xj, ijx in enumerate(calc_check):
				ind_temp = [fi, xj]
				ind_temp.sort()
				# x = [val]
				if ijf[0][0:2] == ijx[0][0:2] and ijf[1][0:2] == ijx[1][0:2] and ijf[2][0:2] == ijx[2][0:2]:
					# als = ijf[0][3] * ijx[0][3]
					gadget_sum = sum([(ijf[z][2] * ijx[z][2]) for z in range(0, len(ijf))])
					if gadget_sum == 2:
						ptre_count += 1
					elif gadget_sum == -2:
						ntre_count += 1
					elif gadget_sum == 6:
						one_count += 1
					elif gadget_sum == 0:
						zero_count += 1
					else:
						print(ijf)
						print(ijx)
						print("Gadget ERROR 1:",gadget_sum, "Tetrad#:",fi,xj)

					div_const = gadget_sum / 6
					# print("****** Gadget calculation ******")
					# print("Calc #:", calc_count)
					# print(div_const)
					# print("G values:", gadget_vals)
					if div_const not in gadget_vals:
						gadget_vals.append(div_const)
				elif ijf[0][0:2] == ijx[0][0:2] and ijf[1][0:2] != ijx[1][0:2]:
					gadget_sum = sum([(ijf[z][2] * ijx[z][2])  for z in [0, 5]])
					if gadget_sum == 2:
						ptre_count += 1
					elif gadget_sum == -2:
						ntre_count += 1
					elif gadget_sum == 6:
						one_count += 1
					elif gadget_sum == 0:
						zero_count += 1
					else:
						print("Gadget ERROR 2:",gadget_sum, "Tetrad#:",fi,xj)

					div_const = gadget_sum / 6
					# print("Calc #:", calc_count)
					if div_const not in gadget_vals:
						gadget_vals.append(div_const)
				elif ijf[0][0:2] != ijx[0][0:2] and ijf[1][0:2] == ijx[1][0:2]:
					# print(ijf, ijx)
					gadget_sum = sum([(ijf[z][2] * ijx[z][2])  for z in [1, 4]])
					if gadget_sum == 2:
						ptre_count += 1
					elif gadget_sum == -2:
						ntre_count += 1
					elif gadget_sum == 6:
						one_count += 1
					elif gadget_sum == 0:
						zero_count += 1
					else:
						print("Gadget ERROR 3:",gadget_sum, "Tetrad#:",fi,xj)

					div_const = gadget_sum / 6
					# print("Calc #:", calc_count)
					if div_const not in gadget_vals:
						gadget_vals.append(div_const)
				elif ijf[0][0:2] != ijx[0][0:2] and ijf[2][0:2] == ijx[2][0:2]:
					gadget_sum = sum([(ijf[z][2] * ijx[z][2]) for z in [2, 3]])
					if gadget_sum == 2:
						ptre_count += 1
					elif gadget_sum == -2:
						ntre_count += 1
					elif gadget_sum == 6:
						one_count += 1
					elif gadget_sum == 0:
						zero_count += 1
					else:
						print("Gadget ERROR 4:",gadget_sum, "Tetrad#:",fi,xj)

					div_const = gadget_sum / 6
					# print("Calc #:", calc_count)
					if div_const not in gadget_vals:
						gadget_vals.append(div_const)
				elif ijf[0][0:2] != ijx[0][0:2] and ijf[1][0:2] != ijx[1][0:2] and ijf[2][0:2] != ijx[2][0:2]:
					gadget_sum = 0
					zero_count += 1

					div_const = gadget_sum / 6
					if div_const not in gadget_vals:
						gadget_vals.append(div_const)
				else:
					print("ERROR**********")
					print(ijf)
					print(ijx)
			print("zero count", zero_count)
			print(" 1/3 count", ptre_count)
			print("-1/3 count", ntre_count)
			print("  1  count", one_count)
			print(gadget_vals)
	else:
		pass

	print("################################################")
	print(" Printing final Gadget values and counts        ")
	print("							")
	print("zero count", zero_count)
	print(" 1/3 count", ptre_count)
	print("-1/3 count", ntre_count)
	print("  1  count", one_count)
	print(gadget_vals)
