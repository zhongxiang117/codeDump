import math


def func_matrix_rotation_Rx(gamma):
    """
    Rx = [[1,   0,     0],
          [0, cosX, -sinX],         ## gamma -- roll
          [0, sinX,  cosX]]
    """
    c = math.cos(gamma)
    s = math.sin(gamma)
    return [[1,0,0],[0,c,-s],[0,s,c]]


def func_matrix_rotation_Ry(beta):
    """
    Ry = [[cosX,  0,  sinX],
          [0,     1,    0],         ## beta -- pitch
          [-sinX, 0,  cosX]]
    """
    c = math.cos(beta)
    s = math.sin(beta)
    return [[c,0,s],[0,1,0],[-s,0,c]]


def func_matrix_rotation_Rz(alpha):
    """
    Rz = [[cosX,  -sinX, 0],
          [sinX,  cosX,  0],        ## alpha -- yaw
          [0,       0,   1]]
    """
    c = math.cos(alpha)
    s = math.sin(alpha)
    return [[c,-s,0],[s,c,0],[0,0,1]]


def func_matrix_cross_33(X,Y):
    """
    3D*3D cross 3D*3D
    return: X cross Y  -->  3D*3D
    """
    v00 = X[0][0]*Y[0][0] + X[0][1]*Y[1][0] + X[0][2]*Y[2][0]
    v01 = X[0][0]*Y[0][1] + X[0][1]*Y[1][1] + X[0][2]*Y[2][1]
    v02 = X[0][0]*Y[0][2] + X[0][1]*Y[1][2] + X[0][2]*Y[2][2]

    v10 = X[1][0]*Y[0][0] + X[1][1]*Y[1][0] + X[1][2]*Y[2][0]
    v11 = X[1][0]*Y[0][1] + X[1][1]*Y[1][1] + X[1][2]*Y[2][1]
    v12 = X[1][0]*Y[0][2] + X[1][1]*Y[1][2] + X[1][2]*Y[2][2]

    v20 = X[2][0]*Y[0][0] + X[2][1]*Y[1][0] + X[2][2]*Y[2][0]
    v21 = X[2][0]*Y[0][1] + X[2][1]*Y[1][1] + X[2][2]*Y[2][1]
    v22 = X[2][0]*Y[0][2] + X[2][1]*Y[1][2] + X[2][2]*Y[2][2]

    return [[v00,v01,v02],
            [v10,v11,v12],
            [v20,v21,v22]]


def func_matrix_cross_13(X,Y):
    """
    1D*3D cross 3D*3D
    return: X cross Y  -->  1D*3D
    """
    a = X[0]*Y[0][0] + X[1]*Y[1][0] + X[2]*Y[2][0]
    b = X[0]*Y[0][1] + X[1]*Y[1][1] + X[2]*Y[2][1]
    c = X[0]*Y[0][2] + X[1]*Y[1][2] + X[2]*Y[2][2]
    return [a,b,c]


def func_matrix_transpose(X):
    """
    3D transpose
    return: T(X)  -->  3D
    """
    return [[X[0][0],X[1][0],X[2][0]],
            [X[0][1],X[1][1],X[2][1]],
            [X[0][2],X[1][2],X[2][2]]]



