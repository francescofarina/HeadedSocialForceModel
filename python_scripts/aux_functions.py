import numpy as np


def map_def():
    segments = {}
    segments[0] = np.array([[0, 0], [0, 25]]).transpose()
    segments[1] = np.array([[5, 0], [5, 25]]).transpose()
    # segments_number = segments.__len__()

    map_walls = np.array([[0,0]])
    for segment in segments:
        map_walls = np.concatenate((map_walls, segments[segment]))
    map_walls = map_walls[1:]

    return map_walls


def parameters_load():
    # SFM Parameters
    tau = 0.5
    A = 2000
    B = 0.08
    Aw = 2000
    Bw = 0.08
    k1 = 1.2 * 10 ** 5
    k2 = 2.4 * 10 ** 5
    # HSFM Parameters
    kd = 500
    ko = 1
    k1g = 200  # forward group cohesion force strength
    k2g = 200  # sideward group cohesion force strength
    d_o = 0.5  # sideward maximum distance from the center of mass
    d_f = 1  # forward maximum distance from the center of mass
    alpha = 3

    return tau, A, B, Aw, Bw, k1, k2, kd, ko, k1g, k2g, d_o, d_f, alpha


def initialization(n_groups, N, rm, rM, mm, mM, v0m, v0M, s, am):
    # Map loading
    map_walls = map_def()
    # Number of walls
    double_num_walls, aux = map_walls.shape
    num_walls = int(double_num_walls / 2)

    v0 = v0m + (v0M - v0m) * np.random.rand(N, 1)  # random desired speed
    v = 0 * np.ones((N, 2))  # initial speed
    th = 2 * np.pi * np.random.rand(N, 1) - np.pi  # initial orientation
    omg = 0  # initial angular velocity

    r = np.empty((N, 1), dtype=float)
    m = np.empty((N, 1), dtype=float)
    group_membership = np.empty((N, 1), dtype=int)
    for i in range(len(n_groups)):  # random radii and masses
        # random radii
        r[sum(n_groups[0: i + 1]) - n_groups[i]: sum(n_groups[0:i + 1])] = np.sort(
            rm + (rM - rm) * np.random.rand(n_groups[i], 1))
        # random masses
        m[sum(n_groups[0: i + 1]) - n_groups[i]: sum(n_groups[0:i + 1])] = np.sort(
            mm + (mM - mm) * np.random.rand(n_groups[i], 1))
        # aux variable
        group_membership[sum(n_groups[0: i + 1]) - n_groups[i]: sum(n_groups[0:i + 1])] = int(i)

    J = 0.5 * r ** 2  # Inertia

    i = 0
    p = {}
    X0 = []
    while i < N:
        gr = int(group_membership[i])
        pos = [s[gr][0] - am + 2 * am * np.random.rand(), s[gr][1] - am + 2 * am * np.random.rand()]
        # minimum distance between pedestrians
        d = []
        for l in range(i):
            d.append(int(np.linalg.norm(pos - np.array(p[l][0:1])) <= r[i] + r[l]))

        # minimum distance from walls
        for l in range(num_walls):
            xp = pos[0]
            yp = pos[1]
            rp = np.array(pos)
            ra = map_walls[2 * l: 2 * l + 2, 0]
            rb = map_walls[2 * l: 2 * l + 2, 1]
            xa = ra[0]
            ya = ra[1]
            xb = rb[0]
            yb = rb[1]
            t = ((xp - xa) * (xb - xa) + (yp - ya) * (yb - ya)) / (((xb - xa) ** 2 + (yb - ya) ** 2))
            t_star = min(max(0, t), 1)
            rh = ra + t_star * (rb - ra)
            d.append(int(np.linalg.norm(rp - rh) <= r[i]))
        if sum(d) == 0:
            p[i] = [pos[0], pos[1], v[i, 0], v[i, 1], r[i], m[i]]
            X0 = np.append(X0, [pos[0], pos[1], th[i], np.linalg.norm(v[i, :]), 0, omg])
            i = i + 1

    return map_walls, num_walls, r, m , J, v0, v, th, omg, group_membership, X0, p


class waypoints_updater():

    def __init__(self,e_seq, e_n, e_ind, e_act, N, n_groups, group_membership):
        self.e_seq = e_seq
        self.e_n = e_n
        self.e_ind = e_ind
        self.e_act =e_act
        self.N = N
        self.group_membership = group_membership
        self.n_groups = n_groups

    def waypoint_update(self, position, coef):
        e = np.zeros((self.N, 2))
        # Determination of the current waypoints
        for i in range(self.N):
            curr_wayp = self.e_act[int(self.group_membership[i])][i - sum(self.n_groups[0:int(self.group_membership[i])])]
            vect = curr_wayp - position[i]
            vect_norm = np.linalg.norm(curr_wayp - position[i])
            e[i, :] = vect / vect_norm
            current_index = self.e_ind[int(self.group_membership[i])][i - sum(self.n_groups[0:int(self.group_membership[i])])]

            if vect_norm <= coef and current_index < self.e_n[int(self.group_membership[i])]-1:
                current_index += 1
                curr_wayp = self.e_seq[int(self.group_membership[i])][:, current_index].transpose()
                vect = curr_wayp - position[i]
                vect_norm = np.linalg.norm(curr_wayp - position[i])
                e[i, :] = vect / vect_norm

                self.e_ind[int(self.group_membership[i])][i - sum(self.n_groups[0:int(self.group_membership[i])])] = current_index
                self.e_act[int(self.group_membership[i])][i - sum(self.n_groups[0:int(self.group_membership[i])])] = curr_wayp

            if vect_norm <= coef and current_index == self.e_n[int(self.group_membership[i])]:
                e[i, :] = ((1 - np.exp(-5 * vect_norm))/(1+np.exp(-5 * vect_norm))) * (vect / vect_norm)

        return e