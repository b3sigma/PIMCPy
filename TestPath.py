#!/bin/env python
import numpy
import pylab
import CalcStatistics
from VMCHelp import *
import random
import numpy
from PIMC import *


numParticles=2
numTimeSlices=5
tau=0.1
lam=0.5
Path=PathClass(numpy.zeros((numTimeSlices,numParticles,3),float),tau,lam)
#Path=PathClass(pylab.load("CanonicalPath.txt"),tau,lam)
Path.SetPotential(ZeroFunction)
Path.SetCouplingConstant(0.0)
print "My 0'th slice of the 1'st particle in the z dimension is ",Path.beads[0,1,2]



