
from __future__ import annotations

# A basic, beginner-friendly implementation in R^2 as espoused in "Quantum computing for everyone
# https://mitpress.mit.edu/9780262539531/quantum-computing-for-everyone/"
# vectors and matrices and their ops can be easily implemented using a library like numpy, 
# but they are explicitly implemented here, for the sake of clarity and transparency

from lin_alg import Vector, OrderedOrthonormalBases

import math, random
from utils import is_close


class Qubit():
    def __init__(self):
        # using "hidden/private" variables (with double scores) to emphasize these are not observable without the measurement process
        # Also we cannot set the state at qubit creation, it is going to be random
        self.__state = Vector.Ket(random.random()-0.5, random.random()-0.5).normalize()

    # Override copying methods, to make sure it cant be copied. qubits can only be observed, not copied
    def __copy__(self):
        # https://en.wikipedia.org/wiki/No-cloning_theorem
        assert False

    def __deepcopy__(self, _):
        # https://en.wikipedia.org/wiki/No-cloning_theorem
        assert False

    def get_measured(self, msmt: Measurement):
        assert type(msmt) == Measurement
        res = msmt.bases * self.__state
        p0 = res.x * res.x
        p1 = res.y * res.y
        #print(p0, p1)
        #import pdb; pdb.set_trace()
        #if is_close(msmt.bases.k0.x, 1.0):
        #    import pdb; pdb.set_trace()
        #if not is_close(msmt.bases.k0.x, 1.0):
        #    if not (is_close(res.x,1) or is_close(res.x,0)):
        #        import pdb; pdb.set_trace()
        assert is_close(p0 + p1, 1.0)
        rnd = random.random()
        #print(f'{rnd:.3f}, {p0:.3f}, {p1:.3f}', self.__state)
        if p0 >= rnd:
            self.__state = msmt.bases.k0
            return 0
        else:
            self.__state = msmt.bases.k1
            return 1




# given an angle create basis
def create_basis_for_direction(theta):
    k0 = Vector.Ket(math.cos(theta/2), -math.sin(theta/2))
    k1 = Vector.Ket(math.sin(theta/2), math.cos(theta/2))
    return OrderedOrthonormalBases(k0, k1)


class Measurement():
    def __init__(self, bases: OrderedOrthonormalBases = None, theta: float = None):
        assert (bases is None) ^ (theta is None) # exactly one of them must be specified
        # either theta or basis must be specified?
        if theta is not None:
            bases = create_basis_for_direction(theta)
        self.bases = bases
    def __repr__(self):
        return f'msmt[{self.bases.__repr__()}]'


        
