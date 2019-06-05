import numpy as np
import os

for file in os.listdir('.'):
  if file.endswith(".csv"):
    data = np.loadtxt(file, skiprows=1)
    n = data.shape[0]
    sigma = data.transpose().dot(data)/n
    output = file.replace(".csv",".sigma")
    print(f"Saving {output}")
    np.savetxt(output, sigma)
