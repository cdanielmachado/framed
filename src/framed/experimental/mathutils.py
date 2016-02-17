__author__ = 'daniel'

from scipy.linalg import svd
from scipy import compress
from numpy import ones, concatenate, array


def nullspace(M, eps=1e-12):
    M = array(M)
    u, s, vh = svd(M)
    padding = M.shape[1]-s.shape[0]
    null_mask = concatenate(((s <= eps), ones((padding,), dtype=bool)), axis=0)
    null_space = compress(null_mask, vh, axis=0)
    return null_space