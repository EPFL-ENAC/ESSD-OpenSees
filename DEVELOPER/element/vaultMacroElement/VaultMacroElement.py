import numpy as np

from Kinematics import *

state = np.array([0])


def update(disp, accel):
    global state

    disp = np.array(disp)
    accel = np.array(accel)

    model = Model_TH() # Test Kinematics import

    state[0] += 1
    print("state: ", state[0]) # Test state update


def get_resisting_force():
    force = np.array([1, 2, 3])
    return force.tolist() # WARNING: must return a list for C++


def get_initial_stiff():
    stiff = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    return stiff.flatten().tolist() # flatten() returns a 1D array, easier to process in C++
