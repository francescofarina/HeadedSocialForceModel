import numpy as np
from aux_functions import parameters_load
import config

def HSFM_forces(X, e, N, map_walls, num_walls, r, m, v0):
    tau, A, B, Aw, Bw, k1, k2, kd, ko, k1g, k2g, d_o, d_f, alpha = parameters_load()

    # Positions and velocities
    position = np.zeros((N, 2))
    vel = np.zeros((N, 2))
    for i in range(N):
        position[i, :] = [X[6 * i], X[6 * i+1]]
        vel[i, :] = [X[6 * i+3] * np.cos(X[6 * i+2]), X[6 * i +3] * np.sin(X[6 * i +2])]

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
        ang[i] = np.arctan2(vect[1], vect[0])
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
                ra = min([map_walls[2*w, 0], map_walls[2*w + 1, 0]], [map_walls[2*w, 1], map_walls[2*w + 1,1]])
                ra = np.array(ra)
                rb = max([map_walls[2*w, 0], map_walls[2*w + 1, 0]], [map_walls[2*w, 1], map_walls[2*w + 1,1]])
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
                    fiw2[i] = fiw2[i] + k1 * (r[i] - diw) * niw
                    fiw3[i] = fiw3[i] - k2 * (r[i] - diw) * (vel[i] * tiw) * tiw

    # Force due to the desire to move as v0
    F1 = fi0

    # Other forces
    F2 = fij1 + fij2 + fij3 + fiw1 + fiw2 + fiw3

    return F1, F2, ang


def HSFM_system(t, X, N, n_groups, map_walls, num_walls, r, m, J, v0, group_membership):
    tau, A, B, Aw, Bw, k1, k2, kd, ko, k1g, k2g, d_o, d_f, alpha = parameters_load()

    # Positions and velocities
    position = np.zeros((N, 2))
    vel = np.zeros((N, 2))
    for i in range(N):
        position[i, :] = [X[6 * i], X[6 * i + 1]]
        vel[i, :] = [X[6 * i + 3] * np.cos(X[6 * i + 2]), X[6 * i + 3] * np.sin(X[6 * i + 2])]

    e = config.waypoints.waypoint_update(position, 1.5)

    # Acting forces
    F0, Fe, ang = HSFM_forces(X, e, N, map_walls, num_walls, r, m, v0)
    FT = F0 + Fe

    # Magnitude of F0
    F_nV=(np.sqrt(np.sum(np.abs(F0)**2, 1)))

    #  desired theta
    thr = np.mod(ang, 2*np.pi).flatten()

    # actual theta
    th = np.mod(X.__getitem__(slice(2, None, 6)), 2*np.pi)

    # angle to rotate
    ang = np.unwrap(th - thr)
    # td = np.vstack((ang, ang+2*np.pi, ang-2*np.pi))
    # I = np.argmin(np.abs(td), 0)

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
        a = ang[i]#td[I[i], i]
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
