import opexebo
import numpy as np

acorr = np.loadtxt('heat_maps_real_ac_py/n1ac.csv', delimiter=',');

grid_stats = opexebo.analysis.grid_score(acorr)

print(grid_stats)