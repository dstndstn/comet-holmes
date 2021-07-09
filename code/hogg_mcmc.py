# hogg_mcmc.py
#   Hogg's MCMC sandbox.
#
# intellectual property:
#   Copyright 2009 David W. Hogg.  All rights reserved.  Watch this space
#   for an open-source license to be inserted later.
#
# comments:
#   - Written for clarity, not speed.  The code is intended to be human-
#     readable.
#
# bugs:
#

from numpy import *
from pylab import *
# pylab must go before numpy.random import
import numpy.random as random
OUTPUTPERIOD = 10000

# the Metropolis-Hastings routine
# inputs:
# - nlinks      - number of links of chain to produce
# - data        - input to "lnposterior" function
# - params      - starting params (can be any shape or type)
# - prior_info  - input to "lnposterior" function (can be anything)
# - lnposterior - function lnp = lnposterior(data, params, prior_info)
# - step_info   - input (and output) to "proposal" function (can be anything)
# - proposal    - function (newparams, new_info) = proposal(params, step_info)
# - [beta]      - inverse temperature; default to 1.0
# outputs (chain, step_info):
# - chain       - the chain, duh
# - step_info   - the output from the last call to proposal(); this permits
#                   the proposal function to keep some statistics or be
#                   adaptive
def metropolis_chain(nlinks, data, params, prior_info,
                     lnposterior, step_info, proposal, beta=1.0,
		     output_period=OUTPUTPERIOD):
	lnp = lnposterior(data, params, prior_info)
	print 'starting lnp=', lnp
	oldlnp = lnp
	oldparams = params
	bestparams = oldparams
	bestlnp = oldlnp
	print 'doing', nlinks, 'links of MCMC...'
	naccept = 0
	chain = []
	link = 0
	while link < nlinks:
		(newparams, step_info) = proposal(oldparams, step_info)
		lnp = lnposterior(data, newparams, prior_info)
		if (beta * (lnp - oldlnp)) > log(random.uniform()):
			# keep new parameters
			chain.append((lnp, newparams))
			oldparams = newparams
			oldlnp = lnp
			naccept += 1
			if lnp > bestlnp:
				bestlnp = lnp
				bestparams = newparams
		else:
			# keep old parameters
			chain.append((oldlnp, oldparams))
		link += 1
		if ((link-1) % output_period == 0):
			print link, bestlnp, bestparams
	return (chain, step_info, naccept)

# make gaussian-blob proposals
# step_info contains w, v outputs from eigh()
# HACK:  inner transpose stuff wrong
def gaussian_proposal(params, step_info):
	npar = len(params)
	(w, v, Q, npropose) = step_info
	dparams = zeros(npar)
	for i in range(npar):
		dparams += random.uniform() * w[i] * v[:,i]
	npropose += 1
	newstep_info = (w, v, Q, npropose)
	newparams = params.copy() + Q * dparams
	return (newparams, newstep_info)

# make axis-aligned proposals
def axis_aligned_proposal(params, step_info):
	(dparams, npropose) = step_info
	indx = find(dparams > 0)
	npar = len(indx)
	indx = indx[random.randint(0, npar)]
	npropose[indx] += 1
	newstep_info = (dparams, npropose)
	newparams = params.copy()
	newparams[indx] += dparams[indx] * random.normal()
	return (newparams, newstep_info)

# HACK: not yet written
def equal_posterior_proposal(lnp, params, chain):
	nlink = len(chain)
	(newlnp, newparams) = chain[random.randint(0, nlink)]
	return (newlnp, newparams)

# some functional testing
if __name__ == '__main__':
	print 'hello world'
