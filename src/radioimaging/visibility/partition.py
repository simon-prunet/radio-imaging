#!/usr/bin/env python

"""\
Code for finding the ideal partitions for some set of visibilities
""" 

__author__ = "Sunrise Wang"
__email__ = "sunrise.wang@oca.eu, sunrisewng@gmail.com"

import math

def compute_ells(cdf, icdf, delta, N, alpha):
	if N > 1:
		last_ells = compute_ells(cdf, icdf, delta, N - 1, alpha)
		last_ells.append(icdf(alpha + cdf(last_ells[-1] - delta)) - delta)

		return last_ells
	else:
		return [icdf(alpha) - delta]

def cost(cdf, icdf, delta, N, alpha):
	assert(N > 1)

	last_ell = compute_ells(cdf, icdf, delta, N - 1, alpha)[-1]

	return cdf(last_ell - delta) + alpha - 1

def get_bin_sizes(cdf, ells, delta):
	bin_sizes = []

	bin_sizes.append(cdf(ells[0] + delta).item())
	print(bin_sizes)
	for i, ell in enumerate(ells[1:]):
		bin_sizes.append(cdf(ell + delta) - cdf(ells[i] - delta))

	bin_sizes.append(1 - cdf(ells[-1] - delta))

	return bin_sizes

def get_partitions(cdf, icdf, delta, N, alpha):
	ells = compute_ells(cdf, icdf, delta, N - 1, alpha)

	return ells, get_bin_sizes(cdf, ells, delta)


def dichotomy(cdf, icdf, delta, N, tolerance=1e-10):
	a = 1/N
	b = 1

	while b-a > tolerance:
		cand = (a+b)/2

		c = cost(cdf, icdf, delta, N, cand)

		if math.isnan(c):
			b = cand
		elif c > 0:
			b = cand
		else:
			a = cand

	alpha = (a+b)/2

	return alpha