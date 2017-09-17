# ******************************************************************************
# Name:    Adinkra 4 x 4 analysis
# Author:  Vadim Korotkikh
# Email:   va.korotki@gmail.com
# Date:    June 2017
#
# Description: Added later
#
# ******************************************************************************
# Library Imports

from numpy.linalg import inv
import numpy as np
import itertools
import time
import sys

# Load scripts
import adinkra_nxn_constructor

# ********************************
def core_scheduler():

	matlist			= adinkra_nxn_constructor.gen_product_matrices(4)
	adinkra_list	= adinkra_nxn_constructor.makeall_adinkras(4,4)

	print(len(matlist))
	selfmat_check(matlist)
	breakdown(matlist)

	adinkra_testing(adinkra_list)
	# pass

# ********************************
# Experimenting with Adinkras build
def adinkra_testing(adinkra_list):

	print("")
	print("Adinkra Composition Testing")
	# for i in range(0,4):
	# 	print(adinkra_list[i])
	nomatch	= 0

	perms_rep	= ['1', '2', '3', '4']
	prep		= perms_rep
	perms_list  = list(itertools.permutations(range(4), 4))

	perms_res	= {}
	perms_tups	= []

	for pm in perms_list:
		rstr = prep[pm[0]] + prep[pm[1]] + prep[pm[2]] + "=" + prep[pm[3]]
		# perms_res['%s' % rstr] = [0]
		perms_res['%s' % rstr]	= []
		# perms_tups.append()

	num_check	= []
	print("Length Adinkra list:", len(adinkra_list))
	for idink, tet in enumerate(adinkra_list):
		# if idink >= 1000:
		# 	pass
		bol_check	= 0
		for ival, pm in enumerate(perms_list):
			tmat1  = np.dot(tet[pm[0]][1], tet[pm[1]][1])
			tmat2  = np.dot(tmat1, tet[pm[2]][1])
			nsign  = np.multiply(tet[pm[3]][1], -1)
			if np.array_equal(tmat2, nsign):
				bol_check	= 1
				rstr = prep[pm[0]] + prep[pm[1]] + prep[pm[2]] + "=" + prep[pm[3]]
				# print("%s relationship found" % rstr)
				# print("Adinkra #", idink)
				if idink not in num_check:
					num_check.append(idink)

				if rstr in perms_res:
					perms_res[rstr].append(idink)
				else:
					print("ERROR")
			elif bol_check == 0 and ival == 23:
				nomatch += 1


	print("Adinkras without ABCD prop:", nomatch)
	print("Adinkra num_check:", len(num_check))
	overlap_rec	= []
	novrlap_rec	= []

	overlap_sets	= []

	for key, value in perms_res.items():
		lvalue	= len(value)
		print(key,":",lvalue)
		for key2, vals in perms_res.items():
			if key != key2:
				overlap = [i for i, j in zip(value, vals) if i == j]

				# if overlap_sets:
				# 	overlap_sets.append(set(overlap))

				print("Overlap", key, "&", key2, "Length", len(overlap))
				if len(overlap) == lvalue:
					ovc_str	= key + " & " + key2
					dbl_str	= key2 + " & " + key
					if ovc_str not in overlap_rec:
						if dbl_str not in overlap_rec:
							overlap_rec.append(ovc_str)
				elif len(overlap) == 0:
					ovc_str	= key + " & " + key2
					dbl_str	= key2 + " & " + key
					if ovc_str not in novrlap_rec:
						if dbl_str not in novrlap_rec:
							novrlap_rec.append(ovc_str)
	print("Number of '100%' overlap groups:", len(overlap_rec))
	for ox in overlap_rec:
		print(ox)
	print("Number of '0%' overlap groups:", len(novrlap_rec))
	for ox in novrlap_rec:
		print(ox)


	# 	tmat1	= np.dot(tet[0][1], tet[1][1])
	# 	tmat2	= np.dot(tmat1, tet[2][1])
	# 	if np.array_equal(tmat2, tet[3][1]):
	# 		# print("ABC=D relationship found")
	# 		# print("Adinkra #", idink)
	# 		# print("L matrices:", tet[0][0], tet[1][0], tet[2][0], tet[3][0])
	# 		abc_d += 1
	# 	else:
	# 		pass
	#
	# print("Number of ABC=D ", abc_d)
	# print("Number of ABD=C ", abd_c)

# ********************************
# Breakdown the L matrix matlist via mathematical operations
def selfmat_check(matlist):

    self_inv	= []
    self_trnp	= []
    self_skew	= []
    self_neg	= []
    orth_mat	= []
    orth_err	= []

    for i, im in enumerate(matlist):

        invsmat = inv(im)
        trnpmat = np.transpose(im)
        negmat  = np.multiply(im, -1)

        if np.array_equal(im, invsmat):
        	self_inv.append(i)
        else:
        	pass

        if np.array_equal(im, trnpmat):
        	self_trnp.append(i)
        else:
        	pass

        if np.array_equal(im, np.multiply(trnpmat, -1)):
        	self_skew.append(i)
        else:
        	pass

        if np.array_equal(im, negmat):
        	self_neg.append(i)
        else:
        	pass

        if np.array_equal(invsmat, trnpmat):
        	orth_mat.append(i)
        else:
        	orth_err.append(i)

    """ Check if skew symmetric matrices and inv / transpose matrices
    overlap	"""
    print("Inverse - Skew Symmetric overlap check")
    testinv 	= [i for i in self_skew if i in self_inv]
    print(testinv)
    print("Transpose - Skew Symmetric overlap check")
    testtrnp 	= [i for i in self_skew if i in self_trnp]
    print(testtrnp)
    print("Negative check")
    print("Number of Negative matrices:", len(self_neg))

    print("Transpose check")
    print("Number of Transpose  matrices:", len(self_trnp))
    print("Inverse check")
    print("Number of Inverse  matrices:", len(self_inv))

    print("Skew check")
    print("Number of Skew-symmetric matrices:", len(self_skew))
    print("Real Orthogonal Matrix check")
    print("Number of R. Orthogonal Matrices:", len(orth_mat))
    print("Number of Non-Orth Matrices:", len(orth_err))
    print("Individual Matrix check finished")
    print("")

# ********************************
# Breakdown the L matrix matlist via mathematical operations
def breakdown(matlist):


	invs_check 	= []
	trnp_check 	= []
	skew_check 	= []
	orth_check	= []

	nonskewlist		= []
	skew_recheck	= []


	for i, im in enumerate(matlist):

		invsmat = inv(im)
		trnpmat = np.transpose(im)
		# print(im)
		for j, jm in enumerate(matlist):
			if i != j:
				temp = [i, j]
				if np.array_equal(invsmat, jm):
					# print("INVERSE")
					# print("i, j:", i, j)
					# temp = [i, j]
					temp.sort()
					if temp not in invs_check:
						invs_check.append(temp)
					else:
						pass
				if np.array_equal(trnpmat, jm):
					# temp = [i, j]
					temp.sort()
					if temp not in trnp_check:
						trnp_check.append(temp)
					else:
						pass
				if np.array_equal(trnpmat, np.multiply(jm, -1)):
					temp.sort()
					if temp not in skew_check:
						skew_check.append(temp)
						if i not in nonskewlist:
							nonskewlist.append(i)
						if j not in nonskewlist:
							nonskewlist.append(j)
					else:
						pass
			elif i == j:
				if np.array_equal(invsmat, trnpmat):
					orth_check.append(i)
				else:
					pass

				if np.array_equal(trnpmat, np.multiply(jm, -1)):

					if i not in skew_recheck:
						skew_recheck.append(i)
					else:
						pass
				else:
					pass

	print("Inverse check")
	print("Number of inverse duplicate matrices:", len(invs_check))

	print("Transpose check")
	print("Number of Transpose duplicate matrices:", len(trnp_check))
	print("Skew check")
	print("Number of Skew-symmetric duplicate matrices:", len(skew_check))
	print("nonskewlist:", len(nonskewlist))
	print(missing_elements(nonskewlist))
	print("")
	print("Num of self Skew-symmetric matrices:", len(skew_recheck))
	print(skew_recheck)
	print("Number of Real Orthogonal Matrices:", len(orth_check))

	# print(trnp_check)
	trnp_check.sort()
	# print("")
	print(trnp_check)
	print("")
	print(invs_check)
	test_list = []
	for xmat in invs_check:
		if [mat for mat in trnp_check if xmat == mat]:
			test_list.append(xmat)
	print("Checking if Inverse and Transpose sets are the same")
	print("Number of Inverse - Transpose matches:", len(test_list))


# ********************************
def missing_elements(L):
    # start, end = L[0], L[-1]
	start, end = 0, 383
	return sorted(set(range(start, end + 1)).difference(L))

# ********************************
# Run main()
if __name__ == "__main__":
	start_time = time.time()

	core_scheduler()
	print("-- Execution time --")
	print("---- %s seconds ----" % (time.time() - start_time))
