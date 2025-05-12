# plot f(r)=erf(r)/r with matplotlib
# save figure as svg file in current directory

import numpy as np
from matplotlib import pyplot as plt
from scipy.special import erf

# Define the function f(r) = erf(r)/r
def f(r):
    return erf(r) / r

# Avoid division by zero at r=0 by using the limit: erf(r)/r → 2/sqrt(π) as r→0
def safe_f(r):
    with np.errstate(divide='ignore', invalid='ignore'):
        result = np.where(r != 0, f(r), 2/np.sqrt(np.pi))
    return result

# Create an array of r values from 0 to 5
r = np.linspace(-5, 5, 500)

# Calculate f(r)
y = safe_f(r)

# Create the plot
plt.figure(figsize=(8, 6))
plt.plot(r, y, label=r'$f(r) = \mathrm{erf}(r)/r$', color='blue')

# Add labels and title
plt.xlabel('r', fontsize=12)
plt.ylabel('f(r)', fontsize=12)
plt.title('Plot of f(r) = erf(r)/r', fontsize=14)

# Add grid and legend
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=12)

# Adjust the x-axis to better show the behavior near zero
plt.xlim(-5, 5)

# Save the figure as SVG
plt.savefig('erf_ratio_plot.svg', format='svg', bbox_inches='tight')

# Show the plot
plt.show()