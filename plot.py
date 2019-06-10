import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

methods = [
    # 'BAOAB',
    # 'memory_BAOAB',
    'middle_IsoK_NHC',
    'middle_NHC',
    'middle_NHL',
    'middle_SINR',
    'middle_SubK_NHC',
    'xi_respa_NHC',
]

# L = 5         # Ratio between slow and fast periods of oscillation
# first = 5
# last = 100
# step = 5

L = 10         # Ratio between slow and fast periods of oscillation
first = 5
last = 100
step = 5

m = 1         # Body mass
kT = 1        # Boltzmann constant times temperature

ksoft = 4*np.pi**2*m/L**2
kstiff = (L**2 - 1)*ksoft/2
K = np.array([[ksoft + kstiff, -kstiff], [-kstiff, ksoft + kstiff]])

def potential(sigma):
    return 1/2*np.multiply(K, sigma[0:2, 0:2]).sum()

Z2 = np.zeros((2,2))
I2 = np.identity(2)
sigma_exact = kT*np.block([[np.linalg.inv(K), Z2], [Z2, m*I2]])

frames = {}
for n in [1, 2, 4]:
    data = {'nrespa': range(first, last+1, step)}
    for method in methods:
        values = []
        files = glob.glob(f'L10/harmonic-{method}-n{n}*.sigma')
        if files:
            for nrespa in range(first, last+1, step):
                file = os.path.join(f'L{L}', f'harmonic-{method}-n{n}-respa{nrespa}.sigma')
                sigma = np.loadtxt(file)
                U = potential(sigma)
                values.append(U - 1.0)
            data[method] = values
    frames[f'n={n}'] = pd.DataFrame(data)

for title, frame in frames.items():
    frame.plot(x='nrespa')
    # plt.yscale('log')

plt.show()
