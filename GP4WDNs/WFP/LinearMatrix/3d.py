from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
nbins = 50

c='r'
z=30
ys = np.random.normal(loc=100, scale=10, size=2000)
hist, bins = np.histogram(ys, bins=nbins)
xs = (bins[:-1] + bins[1:])/2
ax.bar(xs, hist, zs=z, zdir='y', color=c, ec=c, alpha=0.8)


c='g'
z=20
ys = np.random.normal(loc=50, scale=10, size=2000)
hist, bins = np.histogram(ys, bins=nbins)
xs = (bins[:-1] + bins[1:])/2
ax.bar(xs, hist, zs=z, zdir='y', color=c, ec=c, alpha=0.8)


c='b'
z=10
ys = np.random.normal(loc=10, scale=10, size=2000)
hist, bins = np.histogram(ys, bins=nbins)
xs = (bins[:-1] + bins[1:])/2
ax.bar(xs, hist, zs=z, zdir='y', color=c, ec=c, alpha=0.8)


c='y'
z=0
ys = np.random.normal(loc=10, scale=10, size=2000)
hist, bins = np.histogram(ys, bins=nbins)
xs = (bins[:-1] + bins[1:])/2
ax.bar(xs, hist, zs=z, zdir='y', color=c, ec=c, alpha=0.8)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.show()
