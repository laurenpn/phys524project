import numpy as np 
import matplotlib.pyplot as plt
from math import sqrt
import cmath
from scipy.special import spherical_jn as jn, spherical_yn as yn


# variables
V_0 = 60 			# MeV
hbarc = 197.3269 	# MeV * fm
mass = 938 			# MeV/c**2
mu = mass/2			# MeV/c**2
R = 1.45 			# fm, width of potential
Rmax = 20 			# fm for example, anything in this range
pi = cmath.pi


# parameters
numpoints = 5000            # sampling density        
h = Rmax/numpoints 			# fm, step size
xlist = np.linspace(0, Rmax, numpoints)


# functions
def k(E):
    return sqrt(mass*E/hbarc**2)

def V(x): 
	if x < R: 
		return -V_0
	else: return 0 

def K(k,r,l):
	return k**2 - mass*V(r)/hbarc**2 - (l**2+l)/r**2    

# Numerov algorithm:

def u_l(E,l):
    k = sqrt(mass*E/hbarc**2)
    xi = 0
    u_l = []
    
    for x in xlist:
        # Small values of x, u_l ~ x**(l+1)
        if x < 0.05:
            u_l.append(x**(l+1))
        # otherwise use numerov
        else:
            u_rph = (1/(1+K(k,x+h,l)*h**2/12))*(u_l[xi-1]*(2 - K(k,x,l)*5*h**2/6) - u_l[xi-2]*(1 + K(k,x-h,l)*h**2/12))
            u_l.append(u_rph)
        xi += 1
        
    return u_l
V_0 = 0
plt.plot(xlist, u_l(10,0), label="no potential")
V_0 = 60
plt.plot(xlist, u_l(10,0), label="with potential")
plt.legend()
plt.show()
"""
plt.plot(xlist, u_l(1,0), label="1 MeV")
plt.plot(xlist, u_l(10,0), label="10 MeV")
plt.plot(xlist, u_l(100,0), label="100 MeV")
plt.plot(xlist, u_l(200,0), label="200 MeV")

plt.xlabel('fm')
plt.ylabel('Radial wavefunction')
plt.legend()
plt.show()
"""
# Calculating the phase shifts: 

def phase_shift(E,l):
    k = sqrt(mass*E/hbarc**2)
    # r1 and r2 >> R, way outside the potential
    r1 = xlist[4444]
    r2 = xlist[4454]
    
    U = u_l(E,l)
    U1 = U[4444]
    U2 = U[4454]
    al = (U1*r2)/(U2*r1)
    x1 = jn(l, k*r1) - al*jn(l, k*r2)
    x2 = yn(l, k*r1) - al*yn(l, k*r2)
    # return the tangent of delta_0
    if np.arctan2(x1, x2) < 0:
        return np.arctan2(x1, x2) + pi
    else:
        return np.arctan2(x1, x2)

# Comparison with exact solution for delta_0:

def analytic_d0(E):
    K = sqrt(mass*E/hbarc**2)
    k_0 = sqrt(mass*(E+V_0))/hbarc
    if np.arctan2(K*np.tan(k_0*R), k_0) - K*R < 0:
        return np.arctan2(K*np.tan(k_0*R), k_0) - K*R + pi
    else:
        return np.arctan2(K*np.tan(k_0*R), k_0) - K*R
#analytic = [analytic_d0(E) for E in range(1,200)]


"""
del_0 = [phase_shift(E, 0) for E in range(1,200)]
del_1 = [phase_shift(E, 1) for E in range(1,200)]
del_2 = [phase_shift(E, 2) for E in range(1,200)]

plt.plot(del_0, label="S")
plt.plot(del_1, label="P")
plt.plot(del_2, label="D")
#plt.title('L=0 phase shifts vs. energy, V_0=20')
plt.ylabel('Radians')
plt.xlabel('Energy')
plt.legend()
plt.plot()
"""

# Evaluating the total cross section:
def sigma_tot(E,l):
    s = 0
    k = sqrt(mass*E/hbarc**2)
    for i in range(0,l+1):
        s += (2*l+1)*(np.sin(phase_shift(E,l)))**2
        
    return (s*4*pi)/k**2

def sigma(E,l):
    k = sqrt(mass*E/hbarc**2)
    d = phase_shift(E,l)
    s = (2*l+1)*(np.sin(d))**2
    return (s*4*pi)/k**2

# analytic cross section
def sigma_an(E):
    k = sqrt(mass*E/hbarc**2)
    return ((4*pi)*(np.sin(analytic_d0(E)))**2)/k**2

print(np.sin(phase_shift(20,1)))

analytic_sigma = [sigma_an(E) for E in range(1,201)]
"""
sigma0 = [(4*pi/k(E)**2)*(np.sin(phase_shift(E,0)))**2 for E in range(1,201)]
sigma1 = [(4*pi/k(E)**2)*(np.sin(phase_shift(E,1)))**2 for E in range(1,201)]
sigma2 = [(4*pi/k(E)**2)*(np.sin(phase_shift(E,2)))**2 for E in range(1,201)]
#sigmatot = [(sigma(E,0)+sigma(E,1)+sigma(E,2)) for E in range(1,201)]


plt.plot(sigma0, label="S")
plt.plot(sigma1, label="P")
plt.plot(sigma2, label="D")
#plt.plot(sigmatot, label="total")
#plt.plot(sigma0, label="numerical")
#plt.plot(analytic_sigma, label="analytic")

plt.plot(analytic_sigma)
plt.xlabel('Energy')
plt.ylabel('V=60')
plt.legend()
plt.show()


"""







