import matplotlib.pyplot as plt
import numpy as np
from output_file import points
from exported_horizons import radius_points


# point after 30 days to check against horizons data
# end_point = points[1800]

differences = np.empty(30)
counter2 = np.empty(30)
for i in range(0, 30, 1):
    differences[i] = np.fabs(
        points[i * 60][0] - radius_points[i]
    ) / (0.5 * (
        points[i * 60][0] + radius_points[i]
    ))
    counter2[i] = i

print(differences)
"""
[0.00000000e+00 2.82724433e-05 6.92283513e-05 6.44633974e-05
 3.57619849e-05 2.73834753e-04 6.85918457e-04 1.30292233e-03
 2.15133113e-03 3.25391437e-03 4.63033169e-03 6.29764879e-03
 8.27077675e-03 1.05628456e-02 1.31855220e-02 1.61492780e-02
 1.94636197e-02 2.31372790e-02 2.71783760e-02 3.15945541e-02
 3.63930917e-02 4.15809941e-02 4.71650662e-02 5.31519693e-02
 5.95482611e-02 6.63604227e-02 7.35948697e-02 8.12579494e-02
 8.93559239e-02 9.78949369e-02]
 """
plt.figure(figsize=(8, 6))
plt.plot(counter2, differences, label='r percent difference')
plt.xlabel("days")
plt.ylabel("Values")
plt.legend()
plt.show()