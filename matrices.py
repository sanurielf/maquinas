# -*- coding: utf-8 -*-

from cvxopt import matrix, spmatrix
from cvxopt.lapack import getri, gesv



def inv(A):

    n, _ = A.size
    B = matrix(spmatrix([1.0]*n, range(n), range(n)))
    gesv(A, B)

    return B

