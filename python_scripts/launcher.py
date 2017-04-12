import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from HeadedSocialForceModel.python_scripts.aux_functions import parameters_load, initialization, waypoints_updater
from HeadedSocialForceModel.python_scripts.HSFM_functions import HSFM_forces


global waypoints


def HSFM_system(X, t, N, n_groups, map_walls, num_walls, r, m, J, v0, group_membership):
    global waypoints
    tau, A, B, Aw, Bw, k1, k2, kd, ko, k1g, k2g, d_o, d_f, alpha = parameters_load()

    # Positions and velocities
    position = np.zeros((N, 2))
    vel = np.zeros((N, 2))
    for i in range(N):
        position[i,:] = [X[6 * i - 6], X[6 * i - 5]]
        vel[i,:] = [X[6 * i - 3] * np.cos(X[6 * i - 4]), X[6 * i - 3] * np.sin(X[6 * i - 4])]

    e, e_act, e_ind, e_n, e_seq = waypoints.waypoint_update(position, 0.5)

    # Acting forces
    F0, Fe, ang = HSFM_forces(X, e, N, map_walls, num_walls, r, m, v0)
    FT = F0 + Fe

    # Magnitude of F0
    F_nV=(np.sqrt(np.sum(np.abs(F0)**2, 1)))

    #  desired theta
    thr=np.mod(ang, 2*np.pi).flatten()

    # actual theta
    th = np.mod(X.__getitem__(slice(2, None, 6)), 2*np.pi)

    # angle to rotate
    ang = th - thr
    td = np.vstack((ang, ang+2*np.pi, ang-2*np.pi))
    I = np.argmin(np.abs(td), 0)

    dX = np.zeros((6*N,1)).flatten()

    # center o of mass of each group
    ci = {}
    for k in range(len(n_groups)):
        ci[k] = np.array([0, 0])

    for i in range(N):
        ci[int(group_membership[i])] = ci[int(group_membership[i])] + position[i]

    for k in range(len(n_groups)):
        ci[k] = ci[k] / n_groups[k]

    for i in range(N):
        a = td[I[i], i]
        kl = 0.3
        kth = J[i] * kl * F_nV[i]
        kom = J[i] * (1+alpha) * np.sqrt(kl * F_nV[i] / alpha)

        p_i = ci[int(group_membership[i])]-position[i]

        dX[6*i] = X[6*i+3] * np.cos(X[6*i+2]) - X[6*i+4] * np.sin(X[6*i+2])
        dX[6*i+1] = X[6*i+3] * np.sin(X[6*i+2]) + X[6*i+4] * np.cos(X[6*i+2])
        dX[6*i+2] = X[6*i+5]

        # Here we substitute the step function in the definition of the group
        # cohesion forces with a sigmoid
        uf_group = k1g * (1+np.tanh(5*(np.abs(np.dot(p_i, [np.cos(X[6*i+2]), np.sin(X[6*i+2])])-d_f)-3))) * \
                   np.dot(p_i / np.linalg.norm(p_i), [np.cos(X[6*i+2]), np.sin(X[6*i+2])])
        uo_group = k2g * (1+np.tanh(5*(np.abs(np.dot(p_i, [-np.sin(X[6*i+2]), np.cos(X[6*i+2])])-d_o)-3))) * \
                   np.dot(p_i / np.linalg.norm(p_i), [-np.sin(X[6*i+2]), np.cos(X[6*i+2])])

        dX[6*i+3] = 1 / m[i] * (np.dot(FT[i], [np.cos(X[6*i+2]), np.sin(X[6*i+2])]) + uf_group)
        dX[6*i+4] = 1 / m[i] * (ko*np.dot(Fe[i], [-np.sin(X[6*i+2]), np.cos(X[6*i+2])]) - kd * X[6*i+4] + uo_group)
        dX[6*i+5] = 1 / J[i] * (-kth * a - kom * X[6*i+5])

    return dX

# Simulation time
TF = 20
t_fine = 1/30  # timing accuracy in the movie (framerate=30)

tau, A, B, Aw, Bw, k1, k2, kd, ko, k1g, k2g, d_o, d_f, alpha = parameters_load()

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
map_walls, num_walls, r, m, J, v0, v, th, omg, group_membership, X0, p = initialization(n_groups, N, rm, rM, mm, mM, v0m, v0M, s, am)


# Assign the actual position as the current waypoint
e_act = {}
for i in range(len(n_groups)):
    e_act[i] = np.zeros((n_groups[i], 2))
for i in range(N):
    e_act[int(group_membership[i])][i - sum(n_groups[0:int(group_membership[i])])] = p[i][0:2]


waypoints = waypoints_updater(e_seq, e_n, e_ind, e_act, N, n_groups, group_membership)

# TODO: Integration through ODE and debugging
t = np.linspace(0, 0.1, 10)
sol = odeint(HSFM_system, X0, t, args=(N, n_groups, map_walls, num_walls, r, m, J, v0, group_membership))
#############################################

