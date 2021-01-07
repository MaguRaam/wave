import numpy
from matplotlib import pyplot, cm
from math import log

def convergence_plot(Ngpts,l2,linfty,filename):
  fig,ax = pyplot.subplots()
  pyplot.loglog(N, linfty,'-b',label = 'Linfty error',marker = 'o')
  pyplot.loglog(N, l2,'--r',label = 'L2 error',marker = 'o')
  leg = ax.legend()
  pyplot.title("Convergence in space")
  pyplot.xlabel("N")
  pyplot.ylabel("Error")
  pyplot.savefig(filename)

def convergence_rate(error):
  for i in range(1,numpy.size(error)):
    covergence_rate = (log(error[i-1]) - log(error[i]))/log(2)
    print(covergence_rate)

N  = numpy.array([16,32,64,128,256])
linfty = numpy.array([1.134483e-02,1.125525e-03,7.935688e-05,5.105060e-06,3.214497e-07])  
l2     = numpy.array([4.413436e-03,4.148342e-04,2.926878e-05,1.888073e-06,1.189517e-07])
convergence_plot(N,l2,linfty,"convergence")

#L2 rate of convergence
print("L2 Rate of convergence:")
convergence_rate(l2)

#linfty rate of convergence:
print("Linfty Rate of convergence")
convergence_rate(linfty)
