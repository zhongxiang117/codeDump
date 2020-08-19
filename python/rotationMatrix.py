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


def func_matrix_rotation_Rzyx(alpha,beta,gama):
    """
    R = Rz(alpha)Ry(beta)Rx(gama)
    """
    ca = math.cos(alpha)
    sa = math.sin(alpha)
    cb = math.cos(beta)
    sb = math.sin(beta)
    cg = math.cos(gama)
    sg = math.sin(gama)

    R11 = ca*cb
    R12 = ca*sb*sg - sa*cg
    R13 = ca*sb*cg + sa*sg

    R21 = sa*cb
    R22 = sa*sb*sg + ca*cg
    R23 = sa*sb*cg - ca*sg

    R31 = -sb
    R32 = cb*sg
    R33 = cb*cg

    return [[R11, R12, R13],
            [R21, R22, R23],
            [R31, R32, R33]]


def func_matrix_diff_Rzyx(alpha,beta,gama):
    """
    Return: [diff(R(alpha)),  diff(R(beta)),  diff(R(gama))]
    Dimension:     3D            3D               3D
    """
    ca = math.cos(alpha)
    sa = math.sin(alpha)
    cb = math.cos(beta)
    sb = math.sin(beta)
    cg = math.cos(gama)
    sg = math.sin(gama)
    
    dR11_a = -sa*cb
    dR12_a = -sa*sb*sg - ca*cg
    dR13_a = -sa*sb*cg + ca*sg

    dR11_b = -ca*sb
    dR12_b = ca*cb*sg
    dR13_b = ca*cb*cg

    dR11_g = 0
    dR12_g = ca*sb*cg + sa*sg
    dR13_g = -ca*sb*sg + sa*cg

    dR21_a = ca*cb
    dR22_a = ca*sb*sg - sa*cg
    dR23_a = ca*sb*cg + sa*sg

    dR21_b = -sa*sb
    dR22_b = sa*cb*sg
    dR23_b = sa*cb*cg

    dR21_g = 0
    dR22_g = sa*sb*cg - ca*sg
    dR23_g = -sa*sb*sg - ca*cg

    dR31_a = 0
    dR32_a = 0
    dR33_a = 0

    dR31_b = -cb
    dR32_b = -sb*sg
    dR33_b = -sb*cg

    dR31_g = 0
    dR32_g = cb*cg
    dR33_g = -cb*sg

    return [[[dR11_a, dR12_a, dR13_a],[dR21_a, dR22_a, dR23_a],[dR31_a, dR32_a, dR33_a]],
            [[dR11_b, dR12_b, dR13_b],[dR21_b, dR22_b, dR23_b],[dR31_b, dR32_b, dR33_b]],
            [[dR11_g, dR12_g, dR13_g],[dR21_g, dR22_g, dR23_g],[dR31_g, dR32_g, dR33_g]]]


def func_matrix_plus_33(X,Y):
    """
    return 3D + 3D
    """
    R1 = [ X[0][0]+Y[0][0], X[0][1]+Y[0][1], X[0][2]+Y[0][2] ]
    R2 = [ X[1][0]+Y[1][0], X[1][1]+Y[1][1], X[1][2]+Y[1][2] ]
    R3 = [ X[2][0]+Y[2][0], X[2][1]+Y[2][1], X[2][2]+Y[2][2] ]
    
    return [R1,R2,R3]


def func_matrix_minus_33(X,Y):
    """
    return 3D - 3D
    """
    R1 = [ X[0][0]-Y[0][0], X[0][1]-Y[0][1], X[0][2]-Y[0][2] ]
    R2 = [ X[1][0]-Y[1][0], X[1][1]-Y[1][1], X[1][2]-Y[1][2] ]
    R3 = [ X[2][0]-Y[2][0], X[2][1]-Y[2][1], X[2][2]-Y[2][2] ]
    
    return [R1,R2,R3]




