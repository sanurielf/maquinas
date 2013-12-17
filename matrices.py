# -*- coding: utf-8 -*-

from cvxopt import matrix
from cvxopt.lapack import getri, gesv



def inv(A):

    n, _ = A.size

    gesv(A, matrix(1.0, (n, n)))

    return A

