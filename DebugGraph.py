#!/usr/bin/env python

import matplotlib
matplotlib.use('PDF') # Must be before importing matplotlib.pyplot or pylab!
matplotlib.rc('font', size=18, **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)
import matplotlib.pyplot as pl

FFylab = r'$ FF $'
FFxlab = r'$ q^{2} $'
pl.plot([1,2,3,4],[4,5,6,8],label='pie')
pl.xlabel(FFxlab)
pl.ylabel(FFylab)
pl.legend()
pl.savefig('./debugging.pdf')
pl.clf()


FFylab = r'$ FF $'
FFxlab = r'$ q^{2} $'
pl.plot([1,2,3,4],[4,5,6,8],label='pie')
pl.xlabel(FFxlab)
pl.ylabel(FFylab)
pl.legend()
pl.savefig('./debugging2.pdf')
pl.clf()
