# Integral for Fermi-Dirac. See M. Zennaro 2016 for reference

import os
import numpy as np
from scipy.integrate import quad

# These numbers are only used to get an idea of what y to start and end at
# Choose y_min as the y at which a=1e-5 in the neutrino context, and y_max as a=1e2

k = 8.617330350e-5         # boltzman in eV/K
T_g = 2.73                 # photon temp today in K
Gamma_nu = 0.71649         # neutrino/photon temp today (including non-instantaneous decoupling effects)
T_nu = Gamma_nu * T_g      # neutrino temp today
N_nu = 3
m_nu_i = [0.02]*N_nu       # nu mass in eV. give them all the same mass for now
M_nu = sum(m_nu_i)

a_min = 1e-6
a_max = 1e2
log_y_min = np.log10(M_nu/N_nu / (k*T_nu) * a_min)
log_y_max = np.log10(M_nu/N_nu / (k*T_nu) * a_max)


def I(x,y):
    # Integrand of (6)
    return x**2 * np.sqrt(x**2 + y**2) / (1. + np.exp(x))

def dIdy(x,y):
    return x**2 * y / np.sqrt(x**2 + y**2) / (1. + np.exp(x))
    
def d2Idy2(x,y):
    return x**2 / y / (x**2 + y**2) * dIdy(x,y)
    
def solve_integral(integrand, y):
    """ y is the parameter, as defined in (6), assumed to be an np.ndarray. """
    solnarr = np.empty(len(y))
    for i in range(len(y)):
        yy = y[i]
        soln = quad(integrand, 0, np.inf, args=(yy))
        solnarr[i] = soln[0]
    return solnarr

Na = 10000
Ncols = 4
y_arr = np.empty(Na) 
table = np.empty((Na, Ncols))

y_arr = np.logspace(log_y_min,log_y_max,Na)

table[:,0] = y_arr
table[:,1] = solve_integral(I, y_arr)          # F
table[:,2] = solve_integral(dIdy, y_arr)       # F'
table[:,3] = solve_integral(d2Idy2, y_arr)     # F'' 

filename = "Ftable"
array_init = "double Ftable[%d][%d]" % (Ncols, Na)
size_defn = "#define Fsize (%d)\n" % Na

# output to .c file
with open(filename+".c", "w+") as file:
    
    file.write("// y F(y) F'(y) F''(y)\n")
    file.write(array_init)
    file.write(" = {\n\t{")
    np.savetxt(file, table.T, delimiter = ", ", newline = "},\n\t{")
    
    # remove the \t{ from the final line
    file.seek(0, os.SEEK_END)                 # seek to end of file
    file.seek(file.tell() - 2, os.SEEK_SET)   # go backwards 2 bytes (chars)
    
    file.write("};\n")

# make corresponding .h file
with open("../api/fastpm/"+filename+".h", "w+") as file:
    file.write(size_defn)
    file.write("extern ")
    file.write(array_init)
    file.write(";\n")
