import numpy
import pylab
import stats
from VMCHelp import *
import random
import stats
CalcStatistics=stats

class PathClass:
  def __init__(self, NumParticles,NumTimeSlices,tau,lam):
     self.NumParticles=NumParticles
     self.NumTimeSlices=NumTimeSlices
     self.tau=tau
     self.lam=lam
     self.beads=numpy.zeros((self.NumTimeSlices,self.NumParticles,3),float)
     self.sigma=numpy.sqrt((self.lam*self.tau))
     self.oneOverFourLambdaTau=1.0/(4.0*lam*tau)
     self.oneOverFourLambdaTau2=1.0/(4.0*lam*tau*tau)
     print "I have setup the path with a temperature of ",1.0/(tau*NumTimeSlices)
  def SetPotential(self,externalPotentialFunction,internalPotentialFunction):
     self.VeeHelper=internalPotentialFunction
     self.VextHelper=externalPotentialFunction
  def Vee(self,R):
     return self.VeeHelper(self,R)
  def Vext(self,R):
     return self.VextHelper(self,R)
  def KineticAction(self,slice,ptcl):
     disp=self.beads[slice,ptcl]-self.beads[slice+1,ptcl]
     dist2=numpy.dot(disp,disp)
     return dist2*self.oneOverFourLambdaTau
  def PotentialAction(self,slice):
     return self.tau*(self.Vee(slice)+self.Vext(slice))
  def ShiftMove(self):
    slicesToShift=random.randint(0,self.NumTimeSlices-1)
    l=range(slicesToShift,len(self.beads))+range(0,slicesToShift)
    self.beads=self.beads[l].copy()
  def KineticEnergy(self):
    KE=0.0
    for ptcl in range(0,self.NumParticles):
        for slice in range(0,self.NumTimeSlices-1):
           disp=self.beads[slice,ptcl]-self.beads[slice+1,ptcl]
           dist2=numpy.dot(disp,disp)
           KE=KE+3.0/(2.0*self.tau)-dist2*self.oneOverFourLambdaTau2
        disp=self.beads[self.NumTimeSlices-1,ptcl]-self.beads[0,ptcl]
        dist2=numpy.dot(disp,disp)
        KE=KE+3.0/(2.0*self.tau)-dist2*self.oneOverFourLambdaTau2
    return KE/(self.NumTimeSlices+0.0)
  def PotentialEnergy(self):
     PE=0.0
     for slice in range(0,self.NumTimeSlices):
           PE = PE+self.Vee(slice)+self.Vext(slice)
     return PE/(self.NumTimeSlices+0.0)
  def Energy(self):
     return self.PotentialEnergy()+self.KineticEnergy()


def ZeroFunction(self,slice):
   return 0.0

def HarmonicOscillator(self,slice):
       return 0.5*(numpy.dot(self.beads[slice,0],self.beads[slice,0])+numpy.dot(self.beads[slice,1],self.beads[slice,1]))

def CalcDensity(self,DensityHistogram):
    for slice in range(0,self.NumTimeSlices): 
        dist=numpy.sqrt(dot(self.beads[slice,0],self.beads[slice,0]))
        DensityHistogram.add(self.beads[slice,0,0])
        dist=numpy.sqrt(dot(self.beads[slice,1],self.beads[slice,1]))
        DensityHistogram.add(self.beads[slice,1,0])

def PairCorrelationFunction(self,PairHistogram):
    for slice in range(0,self.NumTimeSlices):
        disp=self.beads[slice,0]-self.beads[slice,1]
        dist=sqrt(dot(disp,disp) )
        PairHistogram.add(dist)

PairHistogram=Histogram(0.1,10.0,100)
DensityHistogram=Histogram(-5.0,5.0,100)
def PIMC(numSteps,Path,myMove):
   observableSkip=10
   EnergyTrace=[]
   accepted=0.0
   attempted=0.0
   for steps in range(0,numSteps):
         (accepted,attempted)=myMove(Path,accepted,attempted)
         if steps % observableSkip==0 and steps>100:
             EnergyTrace.append(Path.Energy())
             PairCorrelationFunction(Path,PairHistogram)
             CalcDensity(Path,DensityHistogram)
   print CalcStatistics.Stats(numpy.array(EnergyTrace))
   pylab.plot(EnergyTrace)
   PairHistogram.plotMeNorm("pairme.png")
   DensityHistogram.plotMe("density.png")
   pylab.savefig("broken.png")
   print "Accepted Percentage: ",accepted/attempted



def SingleSliceMove(Path,accepted,attempted):
   attempted=attempted+1.0
   Path.ShiftMove()
   ptclToMove=random.randint(0,Path.NumParticles-1)
   sliceToMove=1
   oldAct=-Path.KineticAction(sliceToMove,ptclToMove)
   oldAct-=Path.KineticAction(sliceToMove-1,ptclToMove)
   oldAct-=Path.PotentialAction(sliceToMove)
   delta=1.0*numpy.array([random.random()-0.5,random.random()-0.5,random.random()-0.5])
   Path.beads[1,ptclToMove]=Path.beads[1,ptclToMove]+delta
   newAct=-Path.KineticAction(sliceToMove,ptclToMove)
   newAct-=Path.KineticAction(sliceToMove-1,ptclToMove)
   newAct-=Path.PotentialAction(sliceToMove)
   if not(newAct-oldAct>numpy.log(random.random())):
      Path.beads[sliceToMove,ptclToMove]=Path.beads[sliceToMove,ptclToMove]-delta
   else:
      accepted=accepted+1.0
   return (accepted,attempted)


def Bisect(Path,ptclToMove,maxStepSize):
   stepSize=maxStepSize
   logSampleProb=0.0
   while stepSize>=2:
        sigma=Path.sigma*numpy.sqrt(stepSize/2)
        for i in numpy.arange(0,maxStepSize,stepSize):
            midVec=(Path.beads[i,ptclToMove]+Path.beads[i+stepSize,ptclToMove])/2.0
            delta=numpy.random.normal(0.0,sigma,3)
            logSampleProb+=numpy.dot(delta,delta)/(2.0*sigma*sigma)
            Path.beads[i+stepSize/2,ptclToMove]=midVec+delta
        stepSize=stepSize/2
   return logSampleProb

def NewToOld(Path,ptclToMove,maxStepSize):
   stepSize=maxStepSize
   logSampleProb=0.0
   while stepSize>=2:
        sigma=Path.sigma*numpy.sqrt(stepSize/2)
        for i in numpy.arange(0,maxStepSize,stepSize):
            midVec=(Path.beads[i,ptclToMove]+Path.beads[i+stepSize,ptclToMove])/2.0
            delta=Path.beads[i+stepSize/2,ptclToMove]-midVec
            logSampleProb+=numpy.dot(delta,delta)/(2.0*sigma*sigma)
        stepSize=stepSize/2
   return logSampleProb

def BisectionMove(Path,accepted,attempted):
   attempted=attempted+1
   maxStepSize=8
   Path.ShiftMove()
   ptclToMove=random.randint(0,1)
#   print "PTCL TO MOVE",ptclToMove                                                  
   oldCoords=Path.beads[0:maxStepSize,ptclToMove].copy()
   T_newToOld=NewToOld(Path,ptclToMove,maxStepSize)


   oldKA=0.0
   for slice in range(0,maxStepSize):
     disp=Path.beads[slice,ptclToMove]-Path.beads[slice+1,ptclToMove]
     dist2=numpy.dot(disp,disp)
     oldKA=oldKA-dist2*Path.oneOverFourLambdaTau
     oldKA=oldKA-Path.PotentialAction(slice)

   T_oldToNew=Bisect(Path,ptclToMove,maxStepSize)


   newKA=0.0
   for slice in range(0,maxStepSize):
      disp=Path.beads[slice,ptclToMove]-Path.beads[slice+1,ptclToMove]
      dist2=numpy.dot(disp,disp)
      newKA=newKA-dist2*Path.oneOverFourLambdaTau
      newKA=newKA-Path.PotentialAction(slice) #Path.V_external(slice)+Path.V_internal(slice)                                                   
   if not((newKA-oldKA-(T_newToOld-T_oldToNew))>numpy.log(random.random())):
     Path.beads[0:maxStepSize,ptclToMove]=oldCoords.copy()
   else:
     accepted=accepted+1
#     print "ERRROR"                                                                 
   return (accepted,attempted)

SingleSliceBisection=BisectionMove

def CoulumbPotential(self,slice): 
   c=self.c
   disp=self.beads[slice,0]-self.beads[slice,1]
   dist=numpy.sqrt(numpy.dot(disp,disp))
   return self.c/(dist+0.1)

def DoubleWell(self,slice):
    a=2.0
    x1=self.beads[slice,0,0]
    y1=self.beads[slice,0,1]
    z1=self.beads[slice,0,2]
    x2=self.beads[slice,1,0]
    y2=self.beads[slice,1,1]
    z2=self.beads[slice,1,2]
    return 1.0*(((x1-a)*(x1-a)+y1*y1+z1*z1)*((x1+a)*(x1+a)+y1*y1+z1*z1)+((x2-a)*(x2-a)+y2*y2+z2*z2)*((x2+a)*(x2+a)+y2*y2+z2*z2))

