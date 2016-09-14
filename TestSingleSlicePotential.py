#!/bin/env python
import numpy
import pylab
import CalcStatistics
from PIMCHelp import *
import random
import numpy
from PIMC import *


tau=0.5
#tau=0.1
lam=0.5
numTimeSlices=5
numParticles=2

n = numpy.array(arange(0.0, 100.0,1.0))
E = 0.5 + n
beta = tau * numTimeSlices
expBetaE = numpy.exp(-beta*E)
print 6.0*sum(E*expBetaE)/sum(expBetaE)
print 'Eclassical = %1.5f' % (12.0/beta)

Path=PathClass(numpy.zeros((numTimeSlices,numParticles,3),float),tau,lam)
Path.SetPotential(HarmonicOscillator)
Path.SetCouplingConstant(0.0)
print PIMC(100000,Path,SingleSliceMove)

