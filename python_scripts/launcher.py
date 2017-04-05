import numpy as np
from aux_functions import map_def


# Map loading
map_walls = map_def()
# Number of walls
double_num_walls, aux = map_walls.shape
num_walls = int(double_num_walls/2)

# Simulation time
TF = 20
t_fine = 1/30  # timing accuracy in the movie (framerate=30)

# SFM Parameters
tau = 0.5
A = 2000
B = 0.08
Aw = 2000
Bw = 0.08
k1 = 1.2*10**5
k2 = 2.4*10**5
# HSFM Parameters
kd = 500
ko = 1
k1g = 200  # forward group cohesion force strength
k2g = 200  # sideward group cohesion force strength
d_o = 0.5  # sideward maximum distance from the center of mass
d_f = 1    # forward maximum distance from the center of mass

# Individual characteristics
# Radius
rm = 0.25  # minimum radius
rM = 0.35  # maximum radius
# Mass
mm = 60  # minimum mass
mM = 90  # maximum mass
# Desired speed
v0m = 1  # minimum velocity
v0M = 1  # maximum velocity