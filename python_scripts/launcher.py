import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
import config
from aux_functions import parameters_load, initialization, waypoints_updater
from HSFM_functions import HSFM_system

# Simulation time
TF = 20
t_fine = 1/30  # time accuracy

tau, A, B, Aw, Bw, k1, k2, kd, ko, k1g, k2g, d_o, d_f, alpha = parameters_load()

# Initial conditions
# Number of individuals in each group
# Define n_i the number of individuals in group i, then
# n_groups = [n_1, n_2, ..., n_N];
n_groups = [4, 4]
# Total number of individuals
N = sum(n_groups)

# s{i} contains the starting point of group i
s = {}
s[0] = [2.5, 0]
s[1] = [2.5, 25]

# waypoints sequence
e_seq = {}
# e_seq{i} contains the points through which the members of group i have to pass
e_seq[0] = np.array([s[0], [4, 10], [2.5, 25]]).transpose()
e_seq[1] = np.array([s[1], [1, 10],  [2.5, 0]]).transpose()

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
map_walls, num_walls, r, m, J, v0, v, th, omg, group_membership, X0, p = initialization(n_groups, N, rm, rM, mm, mM, v0m, v0M, s, am)


# Assign the actual position as the current waypoint
e_act = {}
for i in range(len(n_groups)):
    e_act[i] = np.zeros((n_groups[i], 2))
for i in range(N):
    e_act[int(group_membership[i])][i - sum(n_groups[0:int(group_membership[i])])] = p[i][0:2]


config.waypoints = waypoints_updater(e_seq, e_n, e_ind, e_act, N, n_groups, group_membership)


# System evolution
sol = ode(HSFM_system).set_integrator('dopri5')
t_start = 0.0
t_final = TF
delta_t = t_fine
# Number of time steps: 1 extra for initial condition
num_steps = int(np.floor((t_final - t_start)/delta_t) + 1)
sol.set_initial_value(X0, t_start)
sol.set_f_params(N, n_groups, map_walls, num_walls, r, m, J, v0, group_membership)

t = np.zeros((num_steps, 1))
X = np.zeros((num_steps, N*6))
t[0] = t_start
X[0] = X0
k = 1
while sol.successful() and k < num_steps:
    sol.integrate(sol.t + delta_t)
    t[k] = sol.t
    X[k] = sol.y
    print('time ', sol.t)
    k += 1


# Plotting
colors = {0: "#0066CC", 1: "#006600"}
plt.figure()
# Plot of the walls
for i in range(num_walls):
    plt.plot(map_walls[2*i,:], map_walls[2*i+1,:],'k')
# Starting points
plt.plot(X[0].__getitem__(slice(0, None, 6)), X[0].__getitem__(slice(1, None, 6)), 'ro')
# Trajectories
for i in range(N):
    plt.plot(X[:, 6*i],X[:, 6*i+1], color=colors[int(group_membership[i])])
plt.axis('equal')
plt.savefig('trajectories.eps')


