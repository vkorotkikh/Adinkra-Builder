# ******************************************************************************
#
# Name:    Permutation Matrix Calculation
# Author:  Vadim Korotkikh
# Email:   Vadim.Korotkikh
# Date:    November 2016
# Version: 1.3
#
# Description:
#
# ******************************************************************************


# ******************************************************************************
# Begin Imports
import sys
import math
import time
import numpy as np
import numpy.matlib
import itertools
from numpy import array
from numpy.linalg import inv


# import matrix_outerprod_calc
import alpha_beta_4x4

# ******************************************************************************
# Main() function.
# def main():
	# gen_signm(4)
	# pass

# ******************************************************************************
# Do the final Vij calculation
def calculate_vij_matrices(main_tetrad_list):

	""" Remember that the main_tetrad_ark is a list of lists,
		with each list containing four tuples, with tuples being
		matrix number and the matrices itself. """

	startt_vijcalc = time.time()
	vij_possibilities = []
	vij_possibilities = alpha_beta_4x4.illuminator_of_elfes()
	vij_sixset 		= []

	print("							")
	print("	Calculating Vij matrices")
	print("							")
	vij_alphas 		= []
	vij_betas  		= []
	calc_check		= []
	matc_check		= []

	# vij_matrices	= []
	""" Use ells_list to store all the combinations of six ells / tilde ells
		+-1
	"""
	ells_list		= []
	alpha_list		= []
	temp_ls			= []
	temp_alphals	= []
	print_yn		= 0
	anomaly_switch 	= 0

	for ti, teti in enumerate(main_tetrad_list):
		if print_yn:
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

		""" vij_tempset is just for the alpha/beta strings, ij index for Vij and
			the ell or ~ell Coefficient

			vij_tempmat contains the alpha/beta string and the alpha/beta matrix itself
		"""
		vij_tempset = []
		vij_tempmat	= []

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
				ij_temp = [i, j]
				ij_temp.sort()
				ir = i + 1
				jr = j + 1
				ijstr = str(ir) + str(jr)
				if ij_temp not in temp_combos and i != j:
					# print("Vij matrix i-j vals:", ij_temp)
					# print("Vij matrix i-j vals:", ijstr)
					biglj = lj[1]
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
							# temp_vijmat.append(temp_mat)
							if print_yn:
								print("*************$$$$$$$$$$$$$$$$$$ ")
								print("l-solution found:")
								print(ijx)
							ellint = np.int(1)
							# temp_ls.append(tmint)
							if xi < 3:
								tmp_str = "alpha" + str((xi + 1))
								# print(tmp_str)
								vij_tempset.append([tmp_str, ijstr, ellint])
								vij_tempmat.append([tmp_str, np.divide(ijx,2), ellint])
								alpha_temp.append([tmp_str, ijstr, ellint])
								""" Saving the ells and alpha or beta into
								new list for tracking of ells/tildells		"""
								temp_ls.append((tmp_str, ellint))
								# temp_alphals.append(ellint)
							elif xi >= 3:
								tmp_str = "beta" + str((xi - 2))
								vij_tempset.append([tmp_str, ijstr, ellint])
								vij_tempmat.append([tmp_str, np.divide(ijx,2), ellint])
								beta_temp.append([tmp_str, ijstr, ellint])
								temp_ls.append((tmp_str, ellint))
						elif np.array_equal(temp_mat, ijx_neg):
							tf_bool = 1
							# temp_vijmat.append(temp_mat)
							if print_yn:
								print("*************$$$$$$$$$$$$$$$$$$ ")
								print("l-solution found:")
								print(ijx_neg)
							# xint = (xi + 1) * ( -1)
							ellint = np.int(-1)
							# temp_ls.append(ellint)
							if xi < 3:
								tmp_str = "alpha" + str((xi + 1))
								# print(tmp_str)
								vij_tempset.append([tmp_str, ijstr, ellint])
								vij_tempmat.append([tmp_str, np.divide(ijx,2), ellint])
								alpha_temp.append([tmp_str, ijstr, ellint])
								temp_ls.append((tmp_str, ellint))
								# temp_alphals.append(ellint)
							elif xi >= 3:
								tmp_str = "beta" + str((xi - 2))
								vij_tempset.append([tmp_str, ijstr, ellint])
								vij_tempmat.append([tmp_str, np.divide(ijx,2), ellint])
								beta_temp.append([tmp_str, ijstr, ellint])
								temp_ls.append((tmp_str, ellint))
						else:
							if i != j and tf_bool == 0 and xi >= 5:
								if not(np.array_equal(temp_mat, ijx)) or not np.array_equal(temp_mat, ijx_neg):
									print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ")
									print("Anomaly found:",i,j)
									print(temp_mat)
									anomaly_switch = 1
					tf_bool = 0

		if temp_ls not in ells_list:
			ells_list.append(temp_ls)

		# if temp_alphals not in alpha_list:
		# 	if temp_alphals:
		# 		alpha_list.append(temp_alphals)
		temp_ls			= []
		# temp_alphals	= []
		# vij_matrices.append(temp_vijmat)

		if vij_tempset not in calc_check:
			calc_check.append(vij_tempset)
			matc_check.append(vij_tempmat)

		# Alpha and beta break up.
		if alpha_temp not in vij_alphas:
			vij_alphas.append(alpha_temp)
		elif beta_temp not in vij_betas:
			vij_betas.append(beta_temp)

		beta_temp 	= []
		alpha_temp 	= []
		""" End of main_tetrad_list for loop """

	for vijset in matc_check:
		vij_holo = [np.multiply(m[1],m[2]).astype(int) for m in vijset]
		print("vij_holo alpha recalc")
		for i in vij_holo:
			print(i)

	# print_vijmat_coef(calc_check)

	# print("Length Vij alphas tetrads:", len(vij_alphas))
	# print("length Vij beta tetrads:", len(vij_betas))
	# print("")
	# print("All possible Alpha ells permutations:")
	# for iset in ells_list:
	# 	if iset[0][0].startswith('alpha'):
	# 		print(iset)
	# print("--- Vij Holoramy matrices calculation time ---")
	# print("--- %s seconds ----" % (time.time() - startt_vijcalc))

	for i, vijset in enumerate(calc_check):
		print(vijset)
		print(matc_check[i])

	gadget_vals			= []
	one_count 			= 0
	ptre_count			= 0
	ntre_count			= 0
	zero_count			= 0

	anomaly_switch		= 1
	if not anomaly_switch:
		for xi, ijf in enumerate(calc_check):
			for i in range(0,5):
				# if [
				if ijf[i][0] == matc_check[xi][i][0] and ijf[i][2] == matc_check[xi][i][2]:
					pass
				else:
					print("Some sort of error")
					sys.exit("Equivalency eror check")

			for xj, ijx in enumerate(calc_check):

				div_const = 0
				# x = [val]
				if ijf[0][0:2] == ijx[0][0:2] and ijf[1][0:2] == ijx[1][0:2] and ijf[2][0:2] == ijx[2][0:2]:
					# als = ijf[0][3] * ijx[0][3]
					g_sum = np.sum([np.matmul(ijf[i][1], ijx[i][1]) for i in range(0,5)])
					print("G_sum")
					print(g_sum)
					# g_sum = np.sum( np.matmul(ijf[0][1], ijx[0][1]) + np.matmul(ijf[1][1], ijx[1][1]) + np.matmul(ijf[2][1], ijx[2][1]))

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
						print("Gadget ERROR 1:",gadget_sum, "Tetrad#:",xi,xj)
					div_const = gadget_sum / 6

					# if div_const not in gadget_vals:
					# 	gadget_vals.append(div_const)
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
						print("Gadget ERROR 2:",gadget_sum, "Tetrad#:",xi,xj)

					div_const = gadget_sum / 6
					# print("Calc #:", calc_count)
					# if div_const not in gadget_vals:
					# 	gadget_vals.append(div_const)
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
						print("Gadget ERROR 3:",gadget_sum, "Tetrad#:",xi,xj)

					div_const = gadget_sum / 6
					# print("Calc #:", calc_count)
					# if div_const not in gadget_vals:
					# 	gadget_vals.append(div_const)
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
						print("Gadget ERROR 4:",gadget_sum, "Tetrad#:",xi,xj)

					div_const = gadget_sum / 6
					# print("Calc #:", calc_count)
					# if div_const not in gadget_vals:
					# 	gadget_vals.append(div_const)
				elif ijf[0][0:2] != ijx[0][0:2] and ijf[1][0:2] != ijx[1][0:2] and ijf[2][0:2] != ijx[2][0:2]:
					gadget_sum = 0
					zero_count += 1

					div_const = gadget_sum / 6
					# if div_const not in gadget_vals:
					# 	gadget_vals.append(div_const)
				else:
					print("ERROR**********")
					print(ijf)
					print(ijx)
				if div_const not in gadget_vals:
					gadget_vals.append(div_const)
			# print("zero count", zero_count)
			# print("-1/3 count", ntre_count)
			# print(" 1/3 count", ptre_count)
			# print("  1  count", one_count)
			# print("		")
			print(gadget_vals)
			if xi == 768 or xi == 1536 or xi == 6144:
				print("-- Execution time --")
				print("---- %s seconds ----" % (time.time() - xi_calctime))

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

# ******************************************************************************
# Print the Vij alpha/beta matrix results and the ell, ~ell Coefficients (1, -1)
def print_vijmat_coef(vij_calc_check):

	print("*************$$$$$$$$$$$$$$$$$$ ")
	print("Vij Matrix Coefficients Results:")
	print("")
	for mvals in vij_calc_check:
		if any(x for x in mvals if x[0].startswith('alpha')) and any(x for x in mvals if x[0].startswith('beta')):
			print("MIXED ALPHA_BETA ERROR")
			print(mvals)
		else:
			print(mvals)
			# pass
