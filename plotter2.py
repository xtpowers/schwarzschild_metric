import numpy as np
import matplotlib.pyplot as plt
from output_file import points

point_len = len(points)
r_sol = np.empty(point_len)
theta_sol = np.empty(point_len)
phi_sol = np.empty(point_len)
counter2 = np.empty(point_len)

counter = 0
for point in points:
    r_sol[counter] = point[0]
    theta_sol[counter] = point[1]
    phi_sol[counter] = point[2]
    counter2[counter] = counter
    counter += 1

plt.figure(figsize=(8, 6))
#plt.plot(counter2, r_sol, label='r(tau)')
plt.plot(counter2, r_sol * 1e-8, label='r(tau)')
#plt.plot(counter2, theta_sol, label='theta(tau)')
#plt.plot(counter2, phi_sol, label='phi(tau)')
plt.xlabel("tau")
plt.ylabel("Values")
#plt.yscale("log")  # Set y-axis to logarithmic scale
plt.legend()
plt.show()