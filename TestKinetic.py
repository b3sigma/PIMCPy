#!/bin/env python
import numpy
import pylab
import CalcStatistics
from PIMCHelp import *
import random
import numpy
from PIMC import *


tau=0.1
lam=0.5
Path=PathClass(ReadArray("CanonicalPath.txt"),tau,lam)
Path.SetPotential(ZeroFunction)
Path.SetCouplingConstant(0.0)
print "The value of the kinetic action is ", Path.KineticAction(1,2)



