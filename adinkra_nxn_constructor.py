# ******************************************************************************
# Name:    Adinkra NxN Constructor
# Author:  Vadim Korotkikh
# Email:   va.korotki@gmail.com
# Date:	   February 2017
#
# Description: Algorithm for generating n-colour sized Adinkras of n-node
# matrix representations

from numpy.linalg import inv
import numpy as np
import itertools
import time
import sys

# Makes colums/rows that compose an Identity matrix
def unit_vector(n,i):
	vec		= [0 for j in range(n)]
	vec[i] 	= 1
	return np.array(vec)

# Creates the 24 unsigned permutation matrices from permutation group S4
def gen_permutations(n):
	# bracket_n = list(range(n))
	perms_list = list(itertools.permutations(range(n), n))
	ts = [unit_vector(n,i) for i in range(n)]
	temp = []
	for mi in perms_list:
		columns = [ts[mi[i]] for i in range(n)]
		temp_mat = np.column_stack(columns)
		matint = temp_mat.astype(int)
		temp.append(matint)
	return temp

# Generate all sign permutations of an nxn Identity Matrix
def gen_signm(n):

	n 			= int(n)
	items		= [1] * n
	sign_mat 	= []
	for signs in itertools.product([-1,1], repeat=len(items)):
		temp = np.array([a*sign for a,sign in zip(items,signs)],dtype=int)
		# ptemp.append(temp)
		sign_mat.append(np.diag(temp))
	return sign_mat

# Creates the (384) sign permutation matrices
def gen_product_matrices(n):
	legal_matrices = []
	# from adinkra_tetrad_calc import gen_signm
	sign_pmat = gen_signm(n)
	uperm_mat = gen_permutations(n)
	for x in sign_pmat:
		# t1 = np.asmatrix(x)
		for y in uperm_mat:
			# t2 = np.asmatrix(y)
			# legal_matrices.append(np.matmul(t1,t2))
			legal_matrices.append(np.dot(x,y))
	return legal_matrices

# ********************************
def pairing(matli, matlj):
	ri = np.transpose(matli)
	rj = np.transpose(matlj)
	# tmat = np.dot(matli,rj) + np.dot(matlj,ri)
	rtmat = np.dot(ri,matlj) + np.dot(rj,matli)

	# if np.array_equal(ri, inv(matli)) and np.array_equal(rj, inv(matlj)):
	return (np.count_nonzero(rtmat) == 0)

# Make all Adinkras?
def make_adinkras(k, legal_matrices):
	# """
	# Make all legal good lists of matrices of size k inside legal_matrices
	# """
	""" k is the color number of Adinkra, k=4 is a 4 color, 4 L Matrix Adinkra
		ie a tetrads
	"""
	# main_list	= [[] for i in range(len(legal_matrices))]
	# print(len(main_list))
	""" Preallocate lists """
	xtest_pack	= [None] * 12
	fourpack	= [None] * 4

	if k == 1:
		return [list(l) for l in legal_matrices]
	else:
		adinkra_list = []
		print("Length lmats", len(legal_matrices))

		for i, mat in enumerate(legal_matrices):			# Find all matrix pairs for mat
			# good_mats = [m for m in legal_matrices if pairing(mat,m)]

			# test_mats = [ind[0] for ind in enumerate(legal_matrices) if pairing(mat, ind[1])]
			xtest_pack = [ind for ind in enumerate(legal_matrices) if pairing(mat, ind[1])]
			# main_list[i] = xtest_pack

			for val in xtest_pack:
				# main_list[val[0]] = [(i,mat)]
				fourpack	= [nmat for nmat in xtest_pack if pairing(val[1], nmat[1])]
				# print(len(fourpack))
				for xval in fourpack:
					for lastx in [nmat for nmat in fourpack if pairing(xval[1], nmat[1])]:
						# temp = [(i,mat), val, xval, lastx]
						# temp = [i, val[0], xval[0], lastx[0]]
						adinkra_list.append([(i,mat), val, xval, lastx])
			# for im in good_mats:
			# for m in good_mats:
			# 	good_mats_redux = [mx for mx in good_mats if pairing(m,mx)]
			# 	print("Redux", len(good_mats_redux))
			# 	for rm in good_mats_redux:
			# 		flist = [rx for rx in good_mats_redux if pairing(rm, rx)]
			# 		print("flist:", len(flist))
			# print(len(good_mats))
			# biglist  += [[mat]+mlist for mlist in make_adinkras(k-1, good_mats)]
		return adinkra_list


# ********************************
def makeall_adinkras(k,n):

	# print(len(test_list))
	# return make_matrices(n, gen_product_matrices(n))
	main_tetrad = make_adinkras(k,gen_product_matrices(n))

	return main_tetrad

# ********************************
# Run main()
# start_time = time.time()
#
# makeall_adinkras(4,4)
# print("-- Execution time --")
# print("---- %s seconds ----" % (time.time() - start_time))
