import math
import numpy as np
from numpy import linalg as LA
#-----------------------------------
def diag(M):
    # This function returns the eigenvalues and the eigenvectors.

    E,V = LA.eigh(M)
    return E,V
#------------------------
def second_derivative_FFT():
    # This function returns the second derivative 
    dif2mat=np.zeros((tot_ygrid,tot_ygrid))

    dtheta = (2.0*math.pi)/(tot_ygrid+1)
    K=math.pi/dtheta

    for row in range(tot_ygrid):
        for colm in range(tot_ygrid):
            dif2mat[row,colm]=float(row==colm)*(0.5*red_mass)*(K**2/3.0)\
                *(1+(2.0/tot_ygrid)**2)
            if row!=colm:
                dif2mat[row,colm]=0.5*red_mass*(2.0*K**2/(tot_ygrid)**2)*\
                    ((-1)**(colm-row)/(math.sin(math.pi*(colm-row)/(tot_ygrid))**2))
    return dif2mat
#---------------------------
def FFT_grid():
    #This function returns the FFT grid point.
    theta = np.zeros(tot_ygrid)
    dtheta = (2.0*math.pi)/(tot_ygrid+1)
    i=0
    for y in range(-tot_ygrid/2,tot_ygrid/2):
        theta[i]=y*dtheta
        i+=1
    return theta,dtheta
#------------------------------
def Hamiltonian():
    # This function compute the energy and eigenvectors.
    H_ij = np.zeros((tot_ygrid,tot_ygrid))
    dif2mat = second_derivative_FFT()
    theta,dtheta=FFT_grid()

    for row in range(tot_ygrid):
        for colm in range(tot_ygrid):
            H_ij[row,colm] = dif2mat[row,colm] + float(row==colm) * 0.5*W_0\
                *(1.0-math.cos(theta[colm]))
    energy,coef = diag(H_ij)
    return energy,coef
#-------------------------------
if __name__ == "__main__":

    tot_ygrid=200
    red_mass = 0.002806/27.21138386
    W_0=3.56/27.21138386

    energy,coef=Hamiltonian()
    theta,dtheta = FFT_grid()

    FFT_erg=open("energy-FFT.txt", "w")
    FFT_wf=open("wavefunction.txt", "w")

    for i in range(tot_ygrid):
        wf=coef[i,0]**2/dtheta
        FFT_wf.write(str(theta[i])+ " " + str(wf)+"\n")
        FFT_erg.write(str(energy[i]) + "\n")
    
    FFT_erg.close()
    FFT_wf.close()
