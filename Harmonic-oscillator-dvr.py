import math
import numpy as np
from numpy import linalg as LA
#----------------------------------------------
def diag(M):
    # This function returns the eigenvalues and
    # eigenvectors of a complex Hermitian or 
    # a real symmetric matrix.

    E,V = LA.eigh(M)
    return E,V
#---------------------------------------------
def xgrid():
    # This function returns the grid points x_i
    x_ij = np.zeros((tot_grid,tot_grid))
    for row in range(tot_grid):
        for colm in range(tot_grid):
            x_ij[row,colm]= math.sqrt(float(row+1)/(2.0*mass*freq))* \
                float(row==colm-1) + x_eq * float(row==colm)+ math.sqrt\
                    (float(row)/(2.0*mass*freq))*float(row==colm+1)
    
    x_i,vect = diag(x_ij)

    return x_i,vect
#----------------------------------------------  
def weight():
    # This function returns the weight for each grid point x_i.
    w_i = np.zeros(tot_grid)
    x_i,vect = xgrid()
    for i in range(tot_grid):
        w_i[i]=((mass*freq/math.pi)**(-0.25)*math.exp(0.5*mass*freq*(x_i[i]-x_eq)**2)*vect[0,i])**2

    return w_i
#-----------------------------------------------
def second_derivative():
    # This function returns the second derivavtive matrix.
    dif2mat = np.zeros((tot_grid, tot_grid))
    x_i,vect= xgrid()
    for row in range(tot_grid):
        for colm in range(tot_grid):
            for n in range(tot_grid):
                dif2mat[row,colm] += vect[n,row] * (n+0.5) * vect[n,colm]* -2.0 * mass * freq
            dif2mat[row,colm] +=  mass**2 * freq**2 *(x_i[row]-x_eq)**2 * float(row==colm)           

    return dif2mat
#------------------------------------------------
def Hamiltonian():
    # This function returns the energy and the wavefunction. hbar=1
    H_ij = np.zeros((tot_grid,tot_grid))
    x_i,vect = xgrid()
    dif2mat =second_derivative()

    # Compute the Hamiltonian.
    for row in range(tot_grid):
        for colm in range(tot_grid):
            H_ij[row,colm] = -0.5 * (1.0/mass) * dif2mat[row,colm] + \
                0.5* mass * freq**2 * (x_i[row])**2 * float(row==colm)
    
    energy,coef = diag(H_ij)

    return energy,coef
#-------------------------------------------------
if __name__ == "__main__":
    
    # Input: Units are given in atomic units.
    mass  = 1.0/0.007
    freq = 0.007
    tot_grid = 21
    x_eq =0.0

    energy,coef = Hamiltonian()
    x_i,vect = xgrid()
    dx = weight()

    gs_wf= open("wf-E0.txt", "w")
    erg = open("energy.txt", "w")
    for i in range(tot_grid):
        wf = (coef[i,0])**2/dx[i]
        gs_wf.write(str(x_i[i]) + " " + str(wf) + "\n")
        erg.write(str(energy[i]) + "\n")
    gs_wf.close()
    erg.close()