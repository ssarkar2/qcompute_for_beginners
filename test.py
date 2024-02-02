from lin_alg import Vector, OrderedOrthonormalBases

from basics import Qubit, Measurement
import pytest
from utils import is_close, deg_to_rad



def test_vector():
    v0 = Vector(1.0, 2.0, 'bra')
    v1 = Vector(3.0, 4.0, 'ket')
    assert v1.magnitude() == 5.0
    assert v0.__repr__() == 'bra[x:1.0, y:2.0]'
    assert v0 == v0 and v0 != v1
    assert v0 + v0 == Vector(2.0, 4.0, 'bra')
    assert 2.0*v0 == Vector(2.0, 4.0,'bra') == v0*2.0
    assert is_close(1.0, v0.normalize().magnitude())
    assert is_close(1.0, v0.normalize().magnitude())
    assert v0.normalize() == Vector(0.4472135954999579, 0.8944271909999159, 'bra')

def helper(cls):
    b0 = cls(1.0,2.0)
    b1 = cls(3.0,4.0)
    assert b0.__repr__() == f"{('ket', 'bra')['Bra' in str(cls)]}[x:{b0.x}, y:{b0.y}]"
    assert b0 == b0 and b0 != b1
    assert b0 + b1 == b1 + b0
    assert 2.0 * b0 == b0 * 2.0 == cls(2.0,4.0)


def test_bra():
    helper(Vector.Bra)

def test_ket():
    helper(Vector.Ket)

def test_braket():
    b = Vector.Bra(1.0,2.0)
    k = Vector.Ket(3.0,4.0)
    assert b * k == 11.0
    with pytest.raises(Exception):
        k * b
    with pytest.raises(Exception):
        b * b
    with pytest.raises(Exception):
        k * k
    assert b.t() == Vector.Ket(1.0, 2.0)
    assert b.t().t() == b


def test_mtx():
    k0 = Vector.Ket(1.0,0.0)
    k1 = Vector.Ket(0.0,1.0)
    m = OrderedOrthonormalBases(k0, k1)
    assert m.__repr__() == '[ket[x:1.0, y:0.0]; ket[x:0.0, y:1.0]]'
    with pytest.raises(Exception):
        OrderedOrthonormalBases(k0, k0)
    with pytest.raises(Exception):
        OrderedOrthonormalBases(k0, Vector.Bra(0.0, 1.0))
    assert m * Vector.Ket(2.0,3.0) == Vector.Ket(2.0,3.0)
    with pytest.raises(Exception):
        Vector.Ket(2.0,3.0) * m



def test_qubit():
    q0 = Qubit()
    with pytest.raises(Exception):
        q0.__state

def test_msmt():
    k0 = Vector.Ket(1.0,0.0)
    k1 = Vector.Ket(0.0,1.0)
    m = OrderedOrthonormalBases(k0, k1)
    msmt = Measurement(m)
    q0 = Qubit()
    out = q0.get_measured(msmt)
    assert out in [0, 1]

def test_observe_multi_qubits():
    msmt = Measurement(theta = deg_to_rad(90))  # what if we have other angles?
    num_zero = 0
    num_expts = 10000
    for i in range(num_expts):
        q0 = Qubit()
        out = q0.get_measured(msmt)
        if out == 0:
            num_zero += 1
    assert is_close(0.5, num_zero / num_expts, 0.1)


def test_idempotent_measurement():
    msmt = Measurement(theta = deg_to_rad(90))
    num_expts = 100
    q0 = Qubit()
    for i in range(num_expts):
        out = q0.get_measured(msmt)
        if i == 0:
            res = out
        else:
            assert out == res



@pytest.mark.parametrize("degree,zero_ratio", [(90, 0.5), (60, 0.75)])
def test_serial_measurement(degree,zero_ratio):
    msmt0 = Measurement(theta = deg_to_rad(0))
    msmt90 = Measurement(theta = deg_to_rad(degree))
    num_expts = 1000
    num_zero = 0
    counts = {}
    for i in range(num_expts):
        #print('------------')
        q0 = Qubit()
        out0 = q0.get_measured(msmt0)
        #print('....', q0._Qubit__state)
        out1 = q0.get_measured(msmt90)
        counts[(out0, out1)] = counts.get((out0, out1), 0) + 1
        if out1 == 0:
            num_zero += 1
    # probability that it is 0 after second msmt, given it was 0 after first msmt
    p_m1_0_m0_0 = counts[(0,0)] / (counts[(0,0)] + counts[(0,1)])
    assert is_close(zero_ratio, p_m1_0_m0_0, 0.1)