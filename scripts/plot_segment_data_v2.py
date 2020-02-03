import sys,math,os,subprocess,random
#from PIL import Image
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from scipy import stats
from matplotlib.lines import Line2D
import numpy as np

from scipy.interpolate import UnivariateSpline

# Define some points:
#theta = np.linspace(-3, 2, 40)
#points = np.vstack( (np.cos(theta), np.sin(theta)) ).T

X=[0,1,1.5,4,3,2,3.5,8,10,10,11.5,9,7,5,3,3]
Y=[0,-1,4.5,2,7,5,3,3,2,3,5.5,4,6,8,10,10]
points = np.vstack( (X,Y) ).T

# add some noise:
# points = points + 0.05*np.random.randn(*points.shape)

# Linear length along the line:
distance = np.cumsum( np.sqrt(np.sum( np.diff(points, axis=0)**2, axis=1 )) )
distance = np.insert(distance, 0, 0)/distance[-1]


k=1
s=5
# Build a list of the spline function, one for each dimension:
splines = [UnivariateSpline(distance, coords, k=k, s=s) for coords in points.T]

print( [spl.get_knots() for spl in splines] )
# print( splines.get_coeffs() )

# Computed the spline for the asked distances:
alpha = np.linspace(0, 1, 10)
points_fitted = np.vstack( spl(alpha) for spl in splines ).T

# Graph:
plt.plot(*points.T, 'k', label='original points');
plt.plot(*points_fitted.T, '-r', label='fitted spline k='+str(k)+', s='+str(s));
plt.axis('equal'); plt.legend(); plt.xlabel('x'); plt.ylabel('y');
plt.show()
