import math
import numpy as np
from numpy import linalg as LA
#------------------------------------
def diag(M):
    # This function returns the eigenvalues and
    # eigenvectors of a complex Hermitian or 
    # a real symmetric matrix.

    E,V = LA.eigh(M)

    return E,V
#--------------------------------------
def xgrid():
    # This function returns the grid points x=cos(theta). 

    x_ij = np.zeros((tot_grid,tot_grid))
    for colm in range(1,tot_grid):
        row = colm-1
        x_ij[row,colm] = math.sqrt(float(colm**2-mag_quant_num**2)/float(4*colm**2 -1))
        x_ij[colm,row]=x_ij[row,colm]

    x_i,vect = diag(x_ij)

    return x_i,vect
#---------------------------------------
def weight():
    # This function returns the grid point weight.

    dtheta = np.zeros(tot_grid)
    theta = np.zeros(tot_grid)
    x_i,vect=xgrid()

    for i in range(tot_grid):
        theta[i] = math.acos(x_i[i])
    
    for i in range(tot_grid):
        dtheta[i] = (math.sqrt(float(2.0**(mag_quant_num+1) * math.factorial(mag_quant_num))/ \
            float(math.factorial(math.factorial(2*mag_quant_num+1)))) * \
            (math.sin(theta[i]))**(-1*mag_quant_num)*vect[0,i])**2

    return dtheta
#--------------------------------------- 
def second_derivative():
    # This function returns the second derivavtive matrix.

    dif2_leg=np.zeros((tot_grid,tot_grid))
    x_i,vect = xgrid()

    for row in range(tot_grid):
        for colm in range(tot_grid):
            for n in range(tot_grid):
                dif2_leg[row,colm] -= vect[n,row] * float(n*(n+1)) * vect[n,colm]

    return dif2_leg
#---------------------------------------
def Hamiltonian():

    dif2_leg=second_derivative()
    H_ij = np.zeros((tot_grid,tot_grid))

    for row in range(tot_grid):
        for colm in range(tot_grid):
            H_ij[row,colm] = dif2_leg[row,colm] * -0.5 / (reduce_mass * r**2)

    energy,coef = diag(H_ij)

    return energy,coef
#---------------------------------------
if __name__ == "__main__":
    # Input: Unit are given in atomic units.

    # Total no. of grid points.
    tot_grid = 15
    # Magnetic quantum number is set to be zero.
    mag_quant_num = 0
    # Hamiltonian parameter.
    reduce_mass =9267.5656
    r = 2.98728

    x_i,vect = xgrid()
    dtheta = weight()
    energy,coef = Hamiltonian()

    gs_wf= open("wf-E0.txt", "w")
    erg = open("energy.txt", "w")
    for i in range(tot_grid):
        wf = (coef[i,0])**2/dtheta[i]
        gs_wf.write(str(x_i[i]) + " " + str(wf) + "\n")
        erg.write(str(energy[i]) + "\n")
    gs_wf.close()
    erg.close()