import numpy as np
import scipy.io.matlab
from scipy.linalg import expm, solve
import scipy as sc
from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
#import pandas as pd


'''
Q1 Part C

Header for all arrays:
    [ U38, U39, Np39, Pu39 ]

inputs are:
    Fission = [ o_f_U38, o_f_U39, o_f_Np39, o_f_Pu39 ]
    Capture = " (same format just swap with sigma_c for each isotope)
    Decay   = [0 , -lambda_U39, +lambda_U39-lambda_Np39, +lambda_U39-lambda_Np39]

'''

def create_transmutation_matrices(fission, capture, decay):
    # For HW3 [U38, U39, Np39, Pu39]
    A = [[-(fission[0] + capture[0]), 0, 0, 0], [capture[0],-(fission[1] + capture[1]), 0, 0], [0, 0, -(fission[2] + capture[2]), 0], [0, 0, 0, -(fission[3] + capture[3])]]
    B = [[0,0,0,0], [0, -decay[1],0,0], [0,decay[1],-decay[2],0], [0,0,decay[2],-decay[3]]]

    return A,B


'''
Q1 Part D
Run BU eqs
use delta t not per year***********************

'''
def function_do_burnup(N0, A, B, flux, years):
    time = 0
    N = []
    N.append(N0)
    while time < years:
        N.append(np.dot(((np.array(A) * flux) + B), N[time])) # check matrix math*********************************
        time += 1
    print(f"This is the new comp of isotopes after {years} years: {N}")

    U38 = []
    U39 = []
    Np39 = []
    Pu39 = []
    index = 0
    #N = [float(val) for val in N]
    print(f"this is the comp: {N}")
    while index < years+1:
        U38.append(abs(float(N[index][0])))
        U39.append(abs(float(N[index][1])))
        Np39.append(abs(float(N[index][2])))
        Pu39.append(abs(float(N[index][3])))
        index += 1


    
    #print("This is U38")
    print(U39)
    

    #plt.plot(U38)
    #plt.title(f"U238 after {years} years")
    #plt.yscale('log')
    #plt.show()
    #plt.plot(U39)
    #plt.title(f"U239 after {years} years")
    #plt.yscale('log')
    #plt.show()
    #plt.plot(Np39)
    #plt.title(f"Np239 after {years} years")
    #plt.yscale('log')
    #plt.show()
    #plt.plot(Pu39)
    #plt.title(f"Pu239 after {years} years")
    #plt.yscale('log')
    #plt.show()


    fig, axs = plt.subplots(2, 2, figsize=(10, 8))
    axs[0, 0].plot(U38, marker='o', linestyle='-', color='blue', linewidth=2)
    axs[0, 0].set_title(f"U238 after {years} years")
    axs[0, 0].set_yscale('log')


    axs[0, 1].plot(U39, marker='o', linestyle='-', color='green', linewidth=2)
    axs[0, 1].set_title(f"U239 after {years} years")
    axs[0, 1].set_yscale('log')

    axs[1, 0].plot(Np39, marker='o', linestyle='-', color='red', linewidth=2)
    axs[1, 0].set_title(f"Np239 after {years} years")
    axs[1, 0].set_yscale('log')

    axs[1, 1].plot(Pu39, marker='o', linestyle='-', color='purple', linewidth=2)
    axs[1, 1].set_title(f"Pu239 after {years} years")
    axs[1, 1].set_yscale('log')

# Adjust layout to prevent overlap
plt.tight_layout()
plt.show()



# initalize arrays. All data from the matlab "data-Q1"
o_c = [.2071, .4287, 1.2872, .2789]
o_f = [.0369, .05372, .4823, 1.6214]
decay = [4.915703903743227e-18, 4.926419193745202e-4, 3.405151448232792e-6, 9.11013271367728e-13] *31557600
#N0 = [2.73e-2, 0, 0, 6.81e-3]
N0 = [2.73, 0, 0, 0.681]

# create and print trans. matrices
A,B = create_transmutation_matrices(o_f, o_c, decay)
print(f"This is matrix A: {A}")
print(f"This is matrix B: {B}")


# calculate the BU
# converted 3e15 n/cm^2-s to years or 3e-9 n/bn-s
function_do_burnup(N0, A, B, 31557600 * 3e-9, 10)
