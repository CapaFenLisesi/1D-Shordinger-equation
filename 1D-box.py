"""
Particle in a 1D-box problem, we have let 2m = L = hbar =1.

Adjust the initial value of E to obtain the higher eigenvalues, the wave
function is not good at higher energies.

reference:
http://www.wired.com/2016/03/can-solve-quantum-mechanics-classic-particle-box-problem-code/
"""

import numpy as np
import matplotlib.pyplot as plt

# number of slices divided, affect most accuracy.
N = 1000
# the length of each slice
dx = 1 / N
# the intial value of E, increase to obtain larger eigenvalues
E = 150
# increment of E
dE = 0.02
# initialize psi to zero, psi[0] to psi[N] has N+1 values.
psi = np.zeros(N+1)
# initialize psiDot, intialize psiDot[0] = 0 will lead to trivial result.
psiDot = np.ones(N+1)
# the error allowed
epsPsi = 0.001

# let the loop run
psi[N] = 1
# check if psi(L) = 0
while(abs(psi[N])>epsPsi):
    E = E + dE
    # i should from 0, to N-1, otherwise psi[i+1] would not exist
    for i in range(0,N):
        psi[i+1] = psiDot[i]*dx + psi[i]
        psiDot[i+1] = psiDot[i] -E*psi[i]*dx

# unnormalized totoal probability
p = sum(psi**2)*dx
# the normalization coefficient
c = np.sqrt(1/p)

print(E)
x = np.linspace(0,1,N+1)
plt.plot(x,c*psi,'r')
# the analytical expression
plt.plot(x,np.sqrt(2)*np.sin(4*np.pi*x))
plt.show()
