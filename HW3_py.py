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
    
    while index < years+1:
        U38.append(N[index][0])
        U39.append(N[index][1])
        Np39.append(N[index][2])
        Pu39.append(N[index][3])
        index += 1

    #print("This is U38")
    #print(U38)

    plt.plot(U38)
    plt.title(f"U238 after {years} years")
    plt.show()
    plt.plot(U39)
    plt.title(f"U239 after {years} years")
    plt.show()
    plt.plot(Np39)
    plt.title(f"Np239 after {years} years")
    plt.show()
    plt.plot(Pu39)
    plt.title(f"Pu239 after {years} years")
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










'''
Question 2

'''

col_vec = np.zeros((1506,1))

#Part A
# replace U238 and Pu239 with found values
col_vec[1439] = 2.73e-2 #U
col_vec[1460] = 6.81e-3 #Pu
print(f"This is the changed value for U238 in the 0`s vector: {col_vec[1439]}")
print(f"This is the changed value for Pu239 in the 0`s vector: {col_vec[1460]}")

#Part B
def load_matrix_from_mat(file):
    """
    Loads a 1506x1506 matrix from an Excel file.
    Assumes the matrix is on the given sheet in the Excel file.
    """

    '''
    fission = pd.read_excel(f_file)
    o_f = fission.to_numpy() 
    capture = pd.read_excel(c_file)
    o_c = capture.to_numpy() 
    decay = pd.read_excel(d_file)
    d = decay.to_numpy() 
    '''
    fission = loadmat(file)

    o_c = [.2071, .4287, 1.2872, .2789]
    o_f = [.0369, .05372, .4823, 1.6214]
    d = [4.915703903743227e-18, 4.926419193745202e-4, 3.405151448232792e-6, 9.11013271367728e-13] *31557600
    A = [[-(o_f[0] + o_c[0]), 0, 0, 0], [o_c[0],-(o_f[1] + o_c[1]), 0, 0], [0, 0, -(o_f[2] + o_c[2]), 0], [0, 0, 0, -(o_f[3] + o_c[3])]]
    B = [[0,0,0,0], [0, -d[1],0,0], [0,d[1],-d[2],0], [0,0,d[2],-d[3]]]

    return A,B


# Load matrices from Excel files
A1,B1 = load_matrix_from_mat("data-Q2.mat")

N_results = function_do_burnup(N0, A1, B1, 31557600 * 3e-9, 10)
