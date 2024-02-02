
from __future__ import annotations

import math
from typing import Union

'''
vectors and matrices and their ops can be easily implemented using a library like numpy, 
but they are explicitly implemented here, for the sake of clarity and transparency
'''

def is_close(x, y, tol=0.00001):
    return abs(x-y) < tol


class Vector():
    '''
    Construct using Vector, Vector.Bra or Vector.Ket
    '''
    def __init__(self, x, y, type) -> Vector:
        self.x = x
        self.y = y
        assert type in ['bra', 'ket']
        self.type = type
    @staticmethod
    def Bra(x, y) -> Vector:
        return Vector(x, y, 'bra')
    @staticmethod
    def Ket(x, y) -> Vector:
        return Vector(x, y, 'ket')

    # Unary ops
    def magnitude(self) -> float:
        return math.sqrt(self.x * self.x + self.y * self.y)
    def is_ket(self):
        return self.type == 'ket'
    def is_bra(self):
        return not self.is_ket()
    def is_unit(self):
        return is_close(1.0, self.magnitude())
    def normalize(self) -> Vector:
        return self * (1/self.magnitude())
    def __repr__(self) -> str:
        return f'{self.type}[x:{self.x}, y:{self.y}]'
    def t(self) -> Vector:
        return Vector(self.x, self.y, ('bra', 'ket')[self.type=='bra'])
    
    # Scalar ops
    def scalar_mult(self, a):
        return Vector(a * self.x, a * self.y, self.type)
    def __rmul__(self, other):
        assert type(other) == type(1.1)
        return self.scalar_mult(other)

    # Binary ops with other Vectors
    def __mul__(self, other: Union[float, Vector]) -> float:
        scalar_mul = type(other) == type(1.1)
        if scalar_mul:
            return self.scalar_mult(other)
        else:
            assert self.is_bra() and other.is_ket()
            return other.x * self.x + other.y * self.y
    def __add__(self, other: Vector) -> Vector:
        assert self.type == other.type
        return Vector(self.x + other.x, self.y + other.y, self.type)
    def __eq__(self, other: Vector) -> bool:
        return self.type == other.type and (is_close(self.x, other.x) and is_close(self.y, other.y))
    



# Like a square, orthonormal matrix.
class OrderedOrthonormalBases():
    def __init__(self, k0 : Vector, k1 : Vector):
        assert k0.is_ket() and k1.is_ket()
        assert k0.is_unit() and k1.is_unit()
        assert is_close(k0.t() * k1, 0.0)
        assert is_close(k1.t() * k0, 0.0)
        self.k0 = k0
        self.k1 = k1
    def __repr__(self):
        return f'[{self.k0.__repr__()}; {self.k1.__repr__()}]'
    # generally matrices support scalar mult. but these are orthonormal, and any scalar other 1.0 will disrupt that property
    # hence not supporting scalar mult
    def __mul__(self, other: Vector):
        assert type(other) == Vector
        return Vector(self.k0.t() * other, self.k1.t() * other, 'ket')
