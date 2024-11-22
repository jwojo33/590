import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.linalg import solve

def create_transmutation_matrices():

    data = scipy.io.loadmat('data-Q1.mat')
    fission = data['fission']
    capture = data['capture']
    decay = data['lambda']
    # Matrix A: Represents loss due to fission and capture
    A = np.array([
        [-(capture[0, 0] + fission[0, 0]), 0, 0, 0],
        [capture[0, 0], -(capture[1, 0] + fission[1, 0]), 0, 0],
        [0, 0, -(capture[2, 0] + fission[2, 0]), 0],
        [0, 0, 0, -(capture[3, 0] + fission[3, 0])]
    ])

    # Matrix B: Represents the impact of decay
    B = np.array([
        [0, 0, 0, 0],
        [0, -decay[1, 0], 0, 0],
        [0, decay[1, 0], -decay[2, 0], 0],
        [0, 0, decay[2, 0], -decay[3, 0]]
    ])

    return A, B

def chbv(H, x):

    # Coefficients and poles of the partial fraction expansion
    alpha0 = 0.183216998528140087e-11
    alpha = [
        0.557503973136501826e+02 - 1j * 0.204295038779771857e+03,
        -0.938666838877006739e+02 + 1j * 0.912874896775456363e+02,
        0.469965415550370835e+02 - 1j * 0.116167609985818103e+02,
        -0.961424200626061065e+01 - 1j * 0.264195613880262669e+01,
        0.752722063978321642e+00 + 1j * 0.670367365566377770e+00,
        -0.188781253158648576e-01 - 1j * 0.343696176445802414e-01,
        0.143086431411801849e-03 + 1j * 0.287221133228814096e-03
    ]
    theta = [
        -0.562314417475317895e+01 + 1j * 0.119406921611247440e+01,
        -0.508934679728216110e+01 + 1j * 0.358882439228376881e+01,
        -0.399337136365302569e+01 + 1j * 0.600483209099604664e+01,
        -0.226978543095856366e+01 + 1j * 0.846173881758693369e+01,
        0.208756929753827868e+00 + 1j * 0.109912615662209418e+02,
        0.370327340957595652e+01 + 1j * 0.136563731924991884e+02,
        0.889777151877331107e+01 + 1j * 0.166309842834712071e+02
    ]
    
    p = 7
    theta = [-t for t in theta]
    alpha = [-a for a in alpha]
    
    I = np.eye(H.shape[0], dtype=complex)  # Ensure complex dtype for compatibility
    y = alpha0 * x.astype(complex)  # Initialize `y` with complex type
    
    if np.isrealobj(H) and np.isrealobj(x):
        for i in range(p):
            y += solve(H - theta[i] * I, alpha[i] * x)
        y = y.real  # Convert back to real if inputs were real
    else:
        theta += [np.conj(t) for t in theta]
        alpha += [0.5 * np.conj(a) for a in alpha]
        for i in range(2 * p):
            y += solve(H - theta[i] * I, alpha[i] * x)
    
    # Filter small values to zero
    y[np.abs(y) < 1e-30] = 0
    return y

def do_burnup(N0, A, B, flux, years):
    # Compute the burnup matrix H
    H = (A * flux + B) * 3600 * 24 * 365  # Convert years to seconds
    x = N0

    # Use the chbv function to calculate the final nuclide concentrations
    N = chbv(H, x)
    return N

# Initial fuel
rhorod = 15  # g/cm^3
Zr = 0.1 * rhorod
Ac = 0.9 * rhorod
U238 = 0.8 * Ac
Pu239 = 0.2 * Ac

# Initial concentrations in atoms/barn-cm
N0U238 = U238 / 238.0508 * 6.022e23 * 1e-24
N0Pu239 = Pu239 / 239.052162 * 6.022e23 * 1e-24
N0 = np.array([N0U238, 0, 0, N0Pu239])

flux = 3e-9  # neutrons/barn-s
years = 10  # Number of years

# Create transmutation matrices
A, B = create_transmutation_matrices()

# Perform burnup over time
Nnew = []
for i in range(years):
    N = do_burnup(N0, A, B, flux, 1)  # Burnup for one year
    N0 = N
    Nnew.append(N)

Nnew = np.array(Nnew).T  # Transpose for plotting

# Time array
time = np.arange(1, years + 1)

# Plot results
plt.figure(1)
plt.plot(time, Nnew[0, :], '-or', label='U238')
plt.xlabel('Time [years]')
plt.ylabel('Concentration [atoms/barn-cm]')
plt.title('U238')
plt.yscale('log')
plt.show()

plt.figure(2)
plt.plot(time, Nnew[1, :], '-og', label='U239')
plt.xlabel('Time [years]')
plt.ylabel('Concentration [atoms/barn-cm]')
plt.title('U239')
plt.yscale('log')
plt.show()

plt.figure(3)
plt.plot(time, Nnew[2, :], '-ob', label='Np239')
plt.xlabel('Time [years]')
plt.ylabel('Concentration [atoms/barn-cm]')
plt.title('Np239')
plt.yscale('log')
plt.show()

plt.figure(4)
plt.plot(time, Nnew[3, :], '-om', label='Pu239')
plt.xlabel('Time [years]')
plt.ylabel('Concentration [atoms/barn-cm]')
plt.title('Pu239')
plt.yscale('log')
plt.show()
