#Integral for Fermi-Dirac. Using M. Zennaro

import os
import numpy as np
from scipy.integrate import quad

#these numbers are ONLY used to get an idea of what y to start at. Because I'll use specifically for neutrinos here, I will choose y_min as the y at which a=1e-5 in the neutrino context, and y_max as a=1.

k = 8.617330350e-5         #boltzman in eV/K
T_g = 2.73                 #photon temp today in K
#omega_g0 = 2.47e-5        #this can be determined directly from T_g only! Could type in formula.   

Gamma_nu = 0.71649         #neutrino/photon temp today (including non-instantaneous decoupling effects)
T_nu = Gamma_nu * T_g      #neutrino temp today

N_nu = 3

m_nu_i = [1]*N_nu          #nu mass in eV. give them all the same mass for now
M_nu = sum(m_nu_i)

#for choosing range of y to tabulate 
#maybe testcosmology.c will change paras to change a, so have given some leway around recomb and today
a_min = 1e-6
a_max = 10
log_y_min = np.log10(M_nu/N_nu / (k*T_nu) * a_min)  #because will log space
log_y_max = np.log10(M_nu/N_nu / (k*T_nu) * a_max)


def I(x,y):
    #Integrand of (6)
    return x**2 * np.sqrt(x**2 + y**2) / (1. + np.exp(x))

def dIdy(x,y):
    return x**2 * y / np.sqrt(x**2 + y**2) / (1. + np.exp(x))
    
def d2Idy2(x,y):
    return x**2 / y / (x**2 + y**2) * dIdy(x,y)
    
def solve_integral(integrand, y):
    """y is the para, assumed to be an np.ndarray"""
    #as defined in (6)
    #quad does not vectorize nicely
    solnarr = np.empty(len(y))
    for i in range(len(y)):
        yy = y[i]
        soln = quad(integrand, 0, np.inf, args=(yy))
        solnarr[i] = soln[0]
    return solnarr
    

#def omega_nu_times_E2(a):
    #(11)*h**2 : omega_nu(a) * E(a)**2
#    return 15. / np.pi**4 * N_nu * Gamma_nu**4 * omega_g0 / a**4 * F(M_nu/N_nu / (k*T_nu) * a)

Na = 10000
Ncols = 4
y_arr = np.empty(Na) 
table = np.empty((Na, Ncols))             #want to make dictionary like structure in .c

#put 0 at start of table. maybe end would be cleaner for indexing? meh.
y_arr = np.logspace(log_y_min,log_y_max,Na) #is linspace the best idea? maybe log? maybe do log interp later?

table[:,0] = y_arr
table[:,1] = solve_integral(I, y_arr)          #F
table[:,2] = solve_integral(dIdy, y_arr)       #F'
table[:,3] = solve_integral(d2Idy2, y_arr)     #F'' 

filename = "Ftable"
array_init = "double Ftable[%d][%d]" % (Ncols, Na)   #swapped for .T
size_defn = "#define Fsize (%d)\n" % Na

#output to .c file
with open(filename+".c", "w+") as file:
    
    file.write("// y F(y) F'(y) F''(y)\n")
    file.write(array_init)
    file.write(" = {\n\t{")
    np.savetxt(file, table.T, delimiter = ", ", newline = "},\n\t{")   #Yu said loop and file.write instead of this incase numpy change what this does in the future.
    
    #remove the \t{ from the final line (complicated when not in binary bode)
    file.seek(0, os.SEEK_END)                 # seek to end of file
    file.seek(file.tell() - 2, os.SEEK_SET)   # go backwards 2 bytes (chars)
    
    file.write("};")

#make .h file with matching variable names to .c
#makefile has certain directories to take a .h, otherwise it'll be unhappy, so put it in the api dir with all the other .h's
with open("../api/fastpm/"+filename+".h", "w+") as file:
    file.write(size_defn)
    file.write("extern ")
    file.write(array_init)
    file.write(";")

    #maybe need this too 
    #ifndef nuenergy_h
    #define nuenergy_h


    #endif /* nuenergy_h */

    
#probably just add to make file instead of all this include and .h stuff. This is just a quick fix.