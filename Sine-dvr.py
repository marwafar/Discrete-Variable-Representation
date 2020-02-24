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
    # This function returns the grid point x and the dx.

    dx = float(x_final-x_init)/float(tot_grid-1)
    x_i = np.zeros(tot_grid)

    for i in range(tot_grid):
        x_i[i] = x_init + float(i)*dx

    return x_i,dx
#-----------------------------------------------
def second_derivavtive():
    # This function returns the second derivavtive.

    x_i,dx = xgrid()
    dif2_mat = np.zeros((tot_grid,tot_grid))

    const = math.pi/float(tot_grid+1)
    for row in range(1,tot_grid+1):
        for colm in range(1,tot_grid+1):
            if row!=colm:
                dif2_mat[row-1,colm-1] = -1.0*(math.pi/dx)**2 *2.0*(-1.0)**(row-colm) \
                    /float(tot_grid+1)**2 * math.sin(row*const)* math.sin(colm*const) / \
                        (math.cos(row*const) - math.cos(colm*const))**2
            else:
                dif2_mat[row-1,colm-1] = -1.0*(math.pi/dx)**2*((1.0/3.0)+1.0/ \
                    (6.0*(tot_grid+1)**2) - 1.0/(2.0*(tot_grid+1)**2 * math.sin(row*const)**2))

    return dif2_mat

#-----------------------------------------------
if __name__ == "__main__":

    x_init =1.0
    x_final = 3.0
    tot_grid = 3

    