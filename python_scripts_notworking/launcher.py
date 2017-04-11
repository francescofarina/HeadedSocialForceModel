import numpy as np
from HeadedSocialForceModel.python_scripts.aux_functions import map_def, parameters_load, initialization, waypoints_updater


# Simulation time
TF = 20
t_fine = 1/30  # timing accuracy in the movie (framerate=30)

tau, A, B, Aw, Bw, k1, k2, kd, ko, k1g, k2g, d_o, d_f = parameters_load()

# Initial conditions
# Number of individuals in each group
# Define n_i the number of individuals in group i, then
# n_groups = [n_1, n_2, ..., n_N];
n_groups = [6, 6]
# Total number of individuals
N = sum(n_groups)

# s{i} contains the starting point of group i
s = {}
s[0] = [2.5, 0]
s[1] = [2.5, 25]

# waypoints sequence
e_seq = {}
# e_seq{i} contains the points through which the members of group i have to pass
e_seq[0] = np.array([[3, 10], [2.5, 1000]]).transpose()
e_seq[1] = np.array([[2, 10], [2.5, -1000]]).transpose()

e_n = {}        # Number of waypoints
e_ind = {}      # auxiliary index

for i in range(len(n_groups)):
    aux, e_n[i] = e_seq[i].shape                        # number of waypoints of group i
    e_ind[i] = np.zeros((n_groups[i], 1), dtype=int)    # current waypoint for group i

# "am" represents the amplitude of the starting zone. Fix the starting
# points at least at "am" meters from the walls
am = 2

# Individual characteristics
# Radius
rm = 0.25  # minimum radius
rM = 0.35  # maximum radius
# Mass
mm = 60  # minimum mass
mM = 90  # maximum mass
# Desired speed
v0m = 1  # minimum speed
v0M = 1.2  # maximum speed

# Initialization
map_walls, num_walls, r, m , J, v0, v, th, omg, group_membership, X0, p = initialization(n_groups, N, rm, rM, mm, mM, v0m, v0M, s, am)


# Assign the actual position as the current waypoint
e_act = {}
for i in range(len(n_groups)):
    e_act[i] = np.zeros((n_groups[i], 2))
for i in range(N):
    e_act[int(group_membership[i])][i - sum(n_groups[0:int(group_membership[i])])] = p[i][0:2]

waypoints = waypoints_updater(e_seq, e_n, e_ind, e_act, N, n_groups, group_membership)

state = {'X': X0, 'waypoints': waypoints}

def pass_items():
    return N, n_groups, map_walls, num_walls, r, m, J, v0, group_membership, state



