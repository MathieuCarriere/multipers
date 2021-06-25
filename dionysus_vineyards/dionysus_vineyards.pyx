from cython cimport numeric
from cython cimport int
from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "dionysus_vineyards.hpp":
    vector[vector[vector[double]]] vineyards(vector[vector[double]], string, int, int)

def ls_vineyards(filtrations, complex, discard):
    return vineyards(filtrations, complex, discard, 1)
