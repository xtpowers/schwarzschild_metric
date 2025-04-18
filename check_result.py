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


plt.figure(figsize=(8, 6))
plt.plot(counter2, differences, label='r percent difference')
plt.xlabel("days")
plt.ylabel("Values")
plt.legend()
plt.show()