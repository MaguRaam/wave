import numpy as np 
import matplotlib.pyplot as plt 
from scipy.interpolate import UnivariateSpline

t,K,D = np.loadtxt('ke.csv', unpack = 'true', delimiter=',')
t_Re100,dK_dt_Re100 = np.loadtxt('Re800.csv', unpack = 'true', delimiter=',')

plt.figure(1)
plt.plot(t,K)
plt.xlabel("time")
plt.ylabel("KE")
plt.title("Kinetic Energy")
plt.grid()
plt.savefig("ke.pdf")


plt.figure(2)
plt.plot(t,D)
plt.xlabel("time")
plt.ylabel("DR")
plt.title("Dissipation Rate")
plt.grid()
plt.savefig("dr.pdf")

# Calculate the gradient 

K = K/((2.0*np.pi)**3)
dK_dt = np.gradient(K, t, edge_order=2)


plt.figure(3)
plt.plot(t, -dK_dt, label='Numerical')
plt.plot(t_Re100,dK_dt_Re100, label='Reference (Brachet et. al.)')
plt.xlabel("time")
plt.ylabel('-dK/dt')
plt.title("Kinetic Energy Dissipation Rate")
#plt.ylim((0.0,0.02))
plt.legend()
plt.grid()
plt.title(r'$Re = 100$')
plt.savefig("dK_dt.pdf")

#np.savetxt('dK_dt.csv', np.transpose([t,-dK_dt]),delimiter=',')
