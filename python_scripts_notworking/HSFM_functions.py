import numpy as np
from HeadedSocialForceModel.python_scripts.aux_functions import parameters_load, waypoints_updater
from HeadedSocialForceModel.python_scripts.launcher import pass_items


def HSFM_system(t, state):
    X = state['X']
    waypoints = state['waypoints']
    tau, A, B, Aw, Bw, k1, k2, kd, ko, k1g, k2g, d_o, d_f = parameters_load()
    N, n_groups, map_walls, num_walls, r, m, J, v0, group_membership, X0 = pass_items() #TOGLIERE X0

    # Positions and velocities
    position = np.zeros((N, 2))
    vel = np.zeros((N, 2))
    for i in range(N):
        position[i,:] = [X[6 * i - 6], X[6 * i - 5]]
        vel[i,:] = [X[6 * i - 3] * np.cos(X[6 * i - 4]), X[6 * i - 3] * np.sin(X[6 * i - 4])]

    e, e_act, e_ind, e_n, e_seq = waypoints.waypoint_update(position, 0.5)

    # Acting forces
    F0, Fe, ang = HSFM_forces(X, e, position, vel)
    FT = F0 + Fe

    # Magnitude of F0
    F_nV=(np.sqrt(np.sum(np.abs(F0)**2, 1)))

    #  desired theta
    thr=np.mod(ang, 2*np.pi).flatten()

    # actual theta
    th = np.mod(X.__getitem__(slice(2, None, 6)), 2*np.pi)

    # angle to rotate
    ang = th - thr
    td=np.vstack((ang, ang+2*np.pi, ang-2*np.pi))
    print(td)
    # [~,I]=min(abs(td),[],2);
    #
    # dX=zeros(6*N,1);





def HSFM_forces(X, e, position, vel):
    tau, A, B, Aw, Bw, k1, k2, kd, ko, k1g, k2g, d_o, d_f = parameters_load()
    N, n_groups, map_walls, num_walls, r, m, J, v0, group_membership, X0 = pass_items()  # TOGLIERE X0

    fi0 = np.zeros((N, 2))  # velocity force
    # Interindividual forces
    fij1 = np.zeros((N, 2))   # repulsive
    fij2 = np.zeros((N, 2))   # compression
    fij3 = np.zeros((N, 2))   # friction
    # Obstacles
    fiw1 = np.zeros((N, 2))   # repulsive
    fiw2 = np.zeros((N, 2))   # compression
    fiw3 = np.zeros((N, 2))   # friction
    ang = np.zeros((N,1))
    for i in range(N):
        fi0[i,:] = m[i] * (v0[i] * e[i,:] - vel[i,:]) / tau
        vect = e[i,:]
        ang[i] = np.arctan2(vect[1],vect[0])
        for j in range(N):
            if i != j:
                rij = r[i] + r[j]
                dij = np.linalg.norm(position[i] - position[j])
                nij = (position[i] - position[j]) / dij
                fij1[i] = fij1[i] + A * np.exp((rij - dij) / B) * nij
                if dij < rij:
                    fij2[i] = fij2[i] + k1 * (rij - dij) * nij
                    tij = np.array([-nij[1], nij[0]])
                    dvij = np.dot((vel[j] - vel[i]), tij)
                    fij3[i] = fij3[i] + k2 * (rij - dij) * dvij * tij

            # Walls forces
            for w in range(num_walls):
                xp = position[i,0]
                yp = position[i,1]
                rp = np.array(position[i])
                ra = max([map_walls[2*w, 0], map_walls[2*w + 1,0]], [map_walls[2*w, 1], map_walls[2*w + 1,1]])
                ra = np.array(ra)
                rb = max([map_walls[2*w, 0], map_walls[2*w + 1,0]], [map_walls[2*w, 1], map_walls[2*w + 1,1]])
                rb = np.array(rb)
                xa = ra[0]
                ya = ra[1]
                xb = rb[0]
                yb = rb[1]
                # a point on AB can be parametrized as s(t)=ra+t(tb-ta), t in [0,1]
                # distance from s to p is phi(t)=||s(t)-p||
                # d(phi^2) gives the t which minimizes the distance from p to the
                # line in which AB lives. Since t in [0,1], t_star=min(max(0,t),1);
                # and the distance from p to AB is ||s(t_star)-p||

                t = ((xp - xa) * (xb - xa) + (yp - ya) * (yb - ya)) / ((xb - xa) ** 2 + (yb - ya) ** 2)
                t_star = min(max(0, t), 1)
                rh = ra + t_star * (rb - ra)
                diw = np.linalg.norm(rp - rh)
                niw = (rp - rh) / diw
                tiw = np.array([-niw[0], niw[1]])
                fiw1[i] = fiw1[i] + Aw * np.exp((r[i] - diw) / Bw) * niw
                if diw < r[i]:
                    fiw2[i] = fiw2[i] + k1 * (r(i) - diw) * niw
                    fiw3[i] = fiw3[i] - k2 * (r(i) - diw) * (vel[i] * tiw) * tiw

    # Force due to the desire to move as v0
    F1 = fi0

    # Other forces
    F2 = fij1 + fij2 + fij3 + fiw1 + fiw2 + fiw3

    return F1, F2, ang



N, n_groups, map_walls, num_walls, r, m, J, v0, group_membership, X = pass_items() #TOGLIERE X0


HSFM_system(1, X)



