import numpy as np
import scipy as sp



def preimage_of_subspace(TV,Ws):
    # T:V->W is a linear transformation input as a matrix, and TV is the matrix of the image of basis elements of V under T.
    # Vector spaces V, W, and subspace Ws<W are inputted as numpy arrays of the basis elements.
    # Output is the basis for the preimage of Ws under T.
    mat=np.array([v for v in TV]+[w for w in Ws])

    return sp.linalg.null_space(mat.transpose)


