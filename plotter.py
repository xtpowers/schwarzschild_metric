import matplotlib.pyplot as plt
import numpy as np
from output_file import points
# Fixing random state for reproducibility
np.random.seed(19680801)

point_len = len(points)
r = np.empty(point_len)
theta = np.empty(point_len)
phi = np.empty(point_len)

counter = 0
for point in points:
    r[counter] = point[0]
    theta[counter] = point[2] #flipped
    phi[counter] = point[1] #flipped
    counter += 1

# Convert to Cartesian coordinates
x = r * np.sin(phi) * np.cos(theta)
y = r * np.sin(phi) * np.sin(theta)
z = r * np.cos(phi)

# Create the 3D scatter plot
fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection='3d')

ax.scatter(0,0,0,linewidths=3, alpha=1)
ax.scatter(x, y, z, c=np.arange(len(x)), cmap='viridis', alpha=0.25)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

ax.set_xlim(-1e8,1e8)
ax.set_ylim(-1e8,1e8)
ax.set_zlim(-1e8,1e8)
xlim = ax.get_xlim3d()
ylim = ax.get_ylim3d()
zlim = ax.get_zlim3d()
ax.set_box_aspect((xlim[1]-xlim[0], ylim[1]-ylim[0], zlim[1]-zlim[0]))
plt.show()