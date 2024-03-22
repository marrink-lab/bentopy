import numpy as np


# TODO: There's got to be a nicer way of writing this function.
def sphere_mask(size, r):
    x = np.arange(0, size)
    y = np.arange(0, size)
    z = np.arange(0, size)
    cx = cy = cz = size / 2 - 1
    return (
        (x[:, None, None] - cx) ** 2
        + (y[None, :, None] - cy) ** 2
        + (z[None, None, :] - cz) ** 2
    ) < r**2


def cuboid_mask(width, height, depth, padding):
    np.s_[
        int(padding) : int(width) - int(padding),
        int(padding) : int(height) - int(padding),
        int(padding) : int(depth) - int(padding),
    ]
