# ******************************************************************************
#
# Name:    Outer Product calculation
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

def main():
	# verify_billionspaper()
	pass

# ******************************************************************************
# Main() function.
# def illuminator_of_elles():
def illuminator_of_elfes():
	""" These are the alpha and beta matrices multiplied by 2i
		alpha and beta are results of outerproducts inbetween
		Pauli matrices and Identity matrix They are defined in equtions (4.5)
		in Isaac Chappell II, S. James Gates, Jr - 2012
	"""

	# Alpha - simplified by taking out the i by multiplying the outerproduct by 2i
	alpha1i = np.matrix([[0, 0, 0, 2], [0, 0, 2, 0], [0, -2, 0, 0], [-2, 0, 0, 0]])
	alpha2i = np.matrix([[0, 2, 0, 0], [-2, 0, 0, 0], [0, 0, 0, 2], [0, 0, -2, 0]])
	alpha3i = np.matrix([[0, 0, 2, 0], [0, 0, 0, -2], [-2, 0, 0, 0], [0, 2, 0, 0]])

	# Betas - simplified by taking out the i by multiplication of outerprod by 2i
	beta1i = np.matrix([[0, 0, 0, 2], [0, 0, -2, 0], [0, 2, 0, 0], [-2, 0, 0, 0]])
	beta2i = np.matrix([[0, 0, 2, 0], [0, 0, 0, 2], [-2, 0, 0, 0], [0, -2, 0, 0]])
	beta3i = np.matrix([[0, 2, 0, 0], [-2, 0, 0, 0], [0, 0, 0, -2], [0, 0, 2, 0]])

	# print("alpha 1")
	# print(alpha1i)
	# print("")
	# print("alpha 2")
	# print(alpha2i)
	# print("")
	# print("alpha 3")
	# print(alpha3i)
	# print("")
	# print("beta 1")
	# print(beta1i)
	# print("")
	# print("beta 2")
	# print(beta2i)
	# print("")
	# print("beta 3")
	# print(beta3i)
	# print("")

	# abperm_comb = [ np.multiply(alpha1i,-1), np.multiply(alpha2i,-1), np.multiply(alpha3i,-1), np.multiply(beta1i,-1), np.multiply(beta2i,-1), np.multiply(beta3i,-1)]

	abperm_comb = [alpha1i, alpha2i, alpha3i, beta1i, beta2i, beta3i]
	return abperm_comb

# ******************************************************************************
# alpha & beta calculation
def verify_alphas():
	pass


# ******************************************************************************
# Calculate all the possible +1 and -1 l coefficient combinations
def sign_combimutations():
	combo_l = []

	lsigns = [1,1,1,1,1,1]
	for signs in itertools.product([-1,1], repeat=len(lsigns)):
		temp = np.array([a*sign for a,sign in zip(lsigns,signs)])
		combo_l.append(temp)
		# print(temp)
	# print(len(combo_l))


# main()
