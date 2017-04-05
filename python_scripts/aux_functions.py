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