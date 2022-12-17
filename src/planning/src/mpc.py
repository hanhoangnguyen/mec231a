#!/usr/bin/env python
"""
MPC Controller for MEC231A
Author: Han Nguyen
"""

import numpy as np
import pyomo
import pyomo.environ as pyo
from pyomo.opt import SolverStatus, TerminationCondition
import sys

gst0 = [[-0.0004, -0.0017,  1.    ,  1.0162],
        [-0.9846, -0.1751, -0.0006,  0.1606],
        [ 0.1751, -0.9846, -0.0016,  0.3156],
        [ 0.    ,  0.    ,  0.    ,  1.    ]]

twists = [[ 0.0001, -0.3158, -0.0002, -0.3164, -0.0007, -0.316 , -0.0001],
          [ 0.0004,  0.0004,  0.3168,  0.0009,  0.318 ,  0.0022,  0.3173],
          [-0.    ,  0.0812, -0.1927,  0.4819, -0.0233,  0.8816, -0.1613],
          [ 0.0000,  0.0000,  1.    ,  0.0000,  1.    ,  0.0000,  1.    ],
          [ 0.0000,  1.    ,  0.0000,  1.    ,  0.0000,  1.    ,  0.0000],
          [ 1.    ,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000]]


def solve_cftoc(A, B, P, Q, R, N, x0, yref, U_lim, eps):

    model = pyo.ConcreteModel()
    model.N = N
    model.nx = np.size(A, 0)
    model.nu = np.size(B, 1)
    
    # length of finite optimization problem:
    model.tIDX = pyo.Set( initialize= range(model.N+1), ordered=True )  
    model.xIDX = pyo.Set( initialize= range(model.nx), ordered=True )
    model.uIDX = pyo.Set( initialize= range(model.nu), ordered=True )
    
    # these are 2d arrays:
    model.A = A
    model.B = B
    model.Q = Q 
    model.P = P
    model.R = R
    
    # Create state and input variables trajectory:
    model.x = pyo.Var(model.xIDX, model.tIDX)
    model.u = pyo.Var(model.uIDX, model.tIDX)

    def zeros(rows, cols):
        return [[0 for i in range(cols)] for i in range(rows)]

    def eye(N):
        matrix = zeros(N,N)
        for i in range(N):
            matrix[i][i] = 1
        return matrix

    def dot(A, B):
        if len(A) != len(B):
            print("Dimensions don't match in dot product")
        else:
            product = 0
            for i in range(len(A)):
                product += A[i] * B[i]
        return product

    def add_matrix(A,B):
        C = zeros(len(A), len(A[0]))
        for i in range(len(A)):
            for j in range(len(A[0])):
                C[i][j] = A[i][j] + B[i][j]
        return C

    def add_vector(A,B):
        return [A[i] + B[i] for i in range(len(A))]

    def subtract_matrix(A,B):
        C = zeros(len(A), len(A[0]))
        for i in range(len(A)):
            for j in range(len(A[0])):
                C[i][j] = A[i][j] - B[i][j]
        return C

    def scale_matrix(A,c):
        C = zeros(len(A), len(A[0]))
        for i in range(len(A)):
            for j in range(len(A[0])):
                C[i][j] = A[i][j]*c
        return C

    def scale_vector(A,c):
        return [a*c for a in A]

    def outer(a, b):
        M = len(a)
        N = len(b)
        out = zeros(M, N)
        for i in range(M):
            for j in range(N):
                out[i][j] = a[i] * b[j]
        return out

    def norm_squared(v):
        return sum([x**2 for x in v])

    def assign_values(row_start, row_end, col_start, col_end, mat, value, matrix):
        if matrix:
            value_i = 0
            value_j = 0
            for i in range(row_start, row_end):
                value_j = 0
                for j in range(col_start, col_end):
                    mat[i][j] = value[value_i][value_j]
                    value_j += 1
                value_i += 1
        else:
            value_index = 0
            for i in range(row_start, row_end):
                for j in range(col_start, col_end):
                    mat[i][j] = value[value_index]
                    value_index += 1
        return mat

    def matmul(A, B):
        C = zeros(len(A), len(B[0]))
        for i in range(len(A)):
            for j in range(len(B[0])):
                for k in range(len(B)):
                    C[i][j] += A[i][k] * B[k][j]
        return C

    def matvector(A,v):
        return [dot(A[i],v) for i in range(len(A))]

    def skew_3d(omega):
        return [[0, -omega[2], omega[1]],
                [omega[2], 0, -omega[0]],
                [-omega[1], omega[0], 0]]

    def rotation_3d(omega, theta):
        hat_u = skew_3d(omega)
        theta = theta #* np.linalg.norm(omega) Assume norm of omegas are 1
        hat_u = hat_u #/ np.linalg.norm(omega)
        term1 = eye(3)
        term2 = hat_u
        term3 = pyo.sin(theta)
        term4 = matmul(hat_u, hat_u)
        term5 = 1 - pyo.cos(theta)
        term123 = add_matrix(term1, scale_matrix(term2, term3))
        return add_matrix(term123, scale_matrix(term4, term5))
        
    def homog_3d(xi, theta):
        v = xi[:3]
        w = xi[3:]
        I = eye(3)
        # Translation and rotation
        R = rotation_3d(w, theta)
        p1 = (1/np.square(np.linalg.norm(w)))
        p2 = subtract_matrix(I,R)
        p3 = skew_3d(w)
        p23 = matvector(matmul(p2,p3),v)
        p4 = scale_vector(matvector(outer(w,w),v), theta)
        p234 = add_vector(p23, p4)
        p = scale_vector(p234, p1)
        g = eye(4)
        g = assign_values(0,3,0,3,g,R,matrix=True)
        g = assign_values(0,3,3,4,g,p,matrix=False)
        return g

    def prod_exp(xi, theta):
        g = eye(4)
        for i in range(len(xi[0])):
            xi_i = [xi[j][i] for j in range(len(xi))]
            theta_i = theta[i]
            g_i = homog_3d(xi_i, theta_i)
            g = matmul(g, g_i)
        return g

    def forward_kinematics(joint_angles):
        return matmul(prod_exp(twists,joint_angles), gst0)

    def forward_kinematics_position(joint_angles):
        gst = matmul(prod_exp(twists,joint_angles), gst0)
        return [gst[0][3], gst[1][3], gst[2][3]]

    # Objective:
    def objective_rule(model):
        costY = 0.0
        costU = 0.0
        costTerminal = 0.0
        for t in model.tIDX:
            if t < model.N:
                joint_angles = [model.x[i,t] for i in model.xIDX]
                y = forward_kinematics_position(joint_angles)
                for i in range(3):
                    for j in range(3):
                        costY += (y[i]-yref[i]) * model.Q[i,j] * (y[j]-yref[j]) 
        for t in model.tIDX:
            for i in model.uIDX:
                for j in model.uIDX:
                    if t < model.N:
                        costU += model.u[i, t] * model.R[i, j] * model.u[j, t]
        joint_angles_N = [model.x[i,model.N] for i in model.xIDX]
        y_N = forward_kinematics_position(joint_angles_N)
        for i in range(3):
            for j in range(3):               
                costTerminal += (y_N[i]-yref[i]) * model.P[i,j] * (y_N[j]-yref[j])
        return costY + costU + costTerminal
    
    model.cost = pyo.Objective(rule = objective_rule, sense = pyo.minimize)
    
    # Constraints:
    def dynamic_constraint_rule(model, i, t):
        return model.x[i, t+1] - (sum(model.A[i, j] * model.x[j, t] for j in model.xIDX)
                            +  sum(model.B[i, j] * model.u[j, t] for j in model.uIDX) ) == 0.0 if t < model.N else pyo.Constraint.Skip

    model.equality_constraints = pyo.Constraint(model.xIDX, model.tIDX, rule=dynamic_constraint_rule)

    # Initial Constraints
    model.init_const1 = pyo.Constraint(expr = model.x[0, 0] == x0[0])
    model.init_const2 = pyo.Constraint(expr = model.x[1, 0] == x0[1])
    model.init_const3 = pyo.Constraint(expr = model.x[2, 0] == x0[2])
    model.init_const4 = pyo.Constraint(expr = model.x[3, 0] == x0[3])
    model.init_const5 = pyo.Constraint(expr = model.x[4, 0] == x0[4])
    model.init_const6 = pyo.Constraint(expr = model.x[5, 0] == x0[5])
    model.init_const7 = pyo.Constraint(expr = model.x[6, 0] == x0[6])

    # Input Constraints
    model.constraint1 = pyo.Constraint(model.tIDX,  #Maximum velocity constraint
                                    rule=lambda model, t: model.u[0, t] <= U_lim
                                    if t < N else pyo.Constraint.Skip)
    model.constraint2 = pyo.Constraint(model.tIDX,  #Minimum velocity constraint
                                    rule=lambda model, t: model.u[0, t] >= -U_lim
                                    if t < N else pyo.Constraint.Skip)
    model.constraint3 = pyo.Constraint(model.tIDX,  #Maximum velocity constraint
                                    rule=lambda model, t: model.u[1, t] <= U_lim
                                    if t < N else pyo.Constraint.Skip)
    model.constraint4 = pyo.Constraint(model.tIDX,  #Minimum velocity constraint
                                    rule=lambda model, t: model.u[1, t] >= -U_lim
                                    if t < N else pyo.Constraint.Skip)
    model.constraint5 = pyo.Constraint(model.tIDX,  #Maximum velocity constraint
                                    rule=lambda model, t: model.u[2, t] <= U_lim
                                    if t < N else pyo.Constraint.Skip)
    model.constraint6 = pyo.Constraint(model.tIDX,  #Minimum velocity constraint
                                    rule=lambda model, t: model.u[2, t] >= -U_lim
                                    if t < N else pyo.Constraint.Skip)
    model.constraint7 = pyo.Constraint(model.tIDX,  #Maximum velocity constraint
                                    rule=lambda model, t: model.u[3, t] <= U_lim
                                    if t < N else pyo.Constraint.Skip)
    model.constraint8 = pyo.Constraint(model.tIDX,  #Minimum velocity constraint
                                    rule=lambda model, t: model.u[3, t] >= -U_lim
                                    if t < N else pyo.Constraint.Skip)
    model.constraint9 = pyo.Constraint(model.tIDX,  #Maximum velocity constraint
                                    rule=lambda model, t: model.u[4, t] <= U_lim
                                    if t < N else pyo.Constraint.Skip)
    model.constraint10 = pyo.Constraint(model.tIDX,  #Minimum velocity constraint
                                    rule=lambda model, t: model.u[4, t] >= -U_lim
                                    if t < N else pyo.Constraint.Skip)
    model.constraint11 = pyo.Constraint(model.tIDX,  #Maximum velocity constraint
                                    rule=lambda model, t: model.u[5, t] <= U_lim
                                    if t < N else pyo.Constraint.Skip)
    model.constraint12 = pyo.Constraint(model.tIDX,  #Minimum velocity constraint
                                    rule=lambda model, t: model.u[5, t] >= -U_lim
                                    if t < N else pyo.Constraint.Skip)
    model.constraint13 = pyo.Constraint(model.tIDX,  #Maximum velocity constraint
                                    rule=lambda model, t: model.u[6, t] <= U_lim
                                    if t < N else pyo.Constraint.Skip)
    model.constraint14 = pyo.Constraint(model.tIDX,  #Minimum velocity constraint
                                    rule=lambda model, t: model.u[6, t] >= -U_lim
                                    if t < N else pyo.Constraint.Skip)
    #State Constraints
    # model.constraint15 = pyo.Constraint(model.tIDX, #Joint 0 <= 45 degrees
    #                                     rule=lambda model, t: model.x[0, t] <= 0.785
    #                                     if t <= N else pyo.Constraint.Skip)
    # model.constraint16 = pyo.Constraint(model.tIDX, #Joint 0 >= -45 degrees
    #                                     rule=lambda model, t: model.x[0, t] >= -0.785
    #                                     if t <= N else pyo.Constraint.Skip)
    # model.constraint17 = pyo.Constraint(model.tIDX, #Joint 1 <= 30 degrees
    #                                     rule=lambda model, t: model.x[1, t] <= 0.524
    #                                     if t <= N else pyo.Constraint.Skip)
    # model.constraint18 = pyo.Constraint(model.tIDX, #Joint 1 >= -45 degrees
    #                                     rule=lambda model, t: model.x[1, t] >= -0.785
    #                                     if t <= N else pyo.Constraint.Skip)
    # model.constraint19 = pyo.Constraint(model.tIDX, #Joint 2 <= -1 degrees
    #                                     rule=lambda model, t: model.x[2,t] <= -0.017
    #                                     if t <= N else pyo.Constraint.Skip)
    # model.constraint20 = pyo.Constraint(model.tIDX, #Joint 2 >= -120 degrees
    #                                     rule=lambda model, t: model.x[2,t] >= -2.094
    #                                     if t >= N else pyo.Constraint.Skip)
    # model.constraint21 = pyo.Constraint(model.tIDX, #Joint 3 <= 45 degrees
    #                                     rule=lambda model, t: model.x[3,t] <= 0.785
    #                                     if t <= N else pyo.Constraint.Skip)
    # model.constraint22 = pyo.Constraint(model.tIDX, #Joint 3 >= -90 degrees
    #                                     rule=lambda model, t: model.x[3,t] >= -1.571
    #                                     if t <= N else pyo.Constraint.Skip)
    # model.constraint23 = pyo.Constraint(model.tIDX, #Joint 5 <= 60 degrees
    #                                     rule=lambda model, t: model.x[5,t] <= 1.047
    #                                     if t <= N else pyo.Constraint.Skip)
    # model.constraint24 = pyo.Constraint(model.tIDX, #Joint 5 >= -60 degrees
    #                                     rule=lambda model, t: model.x[5,t] >= -1.047
    #                                     if t <= N else pyo.Constraint.Skip)
    # Orientation Constraint
    model.constraint25 = pyo.Constraint(model.tIDX, #gst[0][0] <= -1 + eps
                                        rule=lambda model, t: forward_kinematics([model.x[i,t] for i in model.xIDX])[0][0] <= -1 + eps
                                        if t <= N else pyo.Constraint.Skip)
    model.constraint26 = pyo.Constraint(model.tIDX, #gst[0][0] >= -1 + eps
                                        rule=lambda model, t: forward_kinematics([model.x[i,t] for i in model.xIDX])[0][0] >= -1 - eps
                                        if t <= N else pyo.Constraint.Skip)
    model.constraint27 = pyo.Constraint(model.tIDX, #gst[1][0] <= 0 + eps
                                        rule=lambda model, t: forward_kinematics([model.x[i,t] for i in model.xIDX])[1][0] <= 0 + eps
                                        if t <= N else pyo.Constraint.Skip)
    model.constraint28 = pyo.Constraint(model.tIDX, #gst[1][0] >= 0 - eps
                                        rule=lambda model, t: forward_kinematics([model.x[i,t] for i in model.xIDX])[1][0] >= 0 - eps
                                        if t <= N else pyo.Constraint.Skip)
    model.constraint29 = pyo.Constraint(model.tIDX, #gst[2][0] >= 0 - eps
                                        rule=lambda model, t: forward_kinematics([model.x[i,t] for i in model.xIDX])[2][0] <= 0 + eps
                                        if t <= N else pyo.Constraint.Skip)
    model.constraint30 = pyo.Constraint(model.tIDX, #gst[2][0] >= 0 - eps
                                        rule=lambda model, t: forward_kinematics([model.x[i,t] for i in model.xIDX])[2][0] >= 0 - eps
                                        if t <= N else pyo.Constraint.Skip)
    # model.constraint31 = pyo.Constraint(model.tIDX, 
    #                                     rule=lambda model, t: forward_kinematics([model.x[i,t] for i in model.xIDX])[0][1] <= 0 + eps
    #                                     if t <= N else pyo.Constraint.Skip)                                       
    # model.constraint32 = pyo.Constraint(model.tIDX, 
    #                                     rule=lambda model, t: forward_kinematics([model.x[i,t] for i in model.xIDX])[0][1] >= 0 - eps
    #                                     if t <= N else pyo.Constraint.Skip)
    # model.constraint33 = pyo.Constraint(model.tIDX, 
    #                                     rule=lambda model, t: forward_kinematics([model.x[i,t] for i in model.xIDX])[1][1] <= 1 + eps
    #                                     if t <= N else pyo.Constraint.Skip)
    # model.constraint34 = pyo.Constraint(model.tIDX, 
    #                                     rule=lambda model, t: forward_kinematics([model.x[i,t] for i in model.xIDX])[1][1] >= -1 - eps
    #                                     if t <= N else pyo.Constraint.Skip)
    # model.constraint35 = pyo.Constraint(model.tIDX,
    #                                     rule=lambda model, t: forward_kinematics([model.x[i,t] for i in model.xIDX])[2][1] <= 0 + eps
    #                                     if t <= N else pyo.Constraint.Skip)
    # model.constraint36 = pyo.Constraint(model.tIDX, 
    #                                     rule=lambda model, t: forward_kinematics([model.x[i,t] for i in model.xIDX])[2][1] >= 0 - eps
    #                                     if t <= N else pyo.Constraint.Skip)
    # model.constraint37 = pyo.Constraint(model.tIDX,
    #                                     rule=lambda model, t: forward_kinematics([model.x[i,t] for i in model.xIDX])[0][2] <= 0 + eps
    #                                     if t <= N else pyo.Constraint.Skip)
    # model.constraint38 = pyo.Constraint(model.tIDX, 
    #                                     rule=lambda model, t: forward_kinematics([model.x[i,t] for i in model.xIDX])[0][2] >= 0 - eps
    #                                     if t <= N else pyo.Constraint.Skip)
    # model.constraint39 = pyo.Constraint(model.tIDX, 
    #                                     rule=lambda model, t: forward_kinematics([model.x[i,t] for i in model.xIDX])[1][2] <= 0 + eps
    #                                     if t <= N else pyo.Constraint.Skip)
    # model.constraint40 = pyo.Constraint(model.tIDX, 
    #                                     rule=lambda model, t: forward_kinematics([model.x[i,t] for i in model.xIDX])[1][2] >= 0 - eps
    #                                     if t <= N else pyo.Constraint.Skip)
    # model.constraint41 = pyo.Constraint(model.tIDX, 
    #                                     rule=lambda model, t: forward_kinematics([model.x[i,t] for i in model.xIDX])[2][2] <= -1 + eps
    #                                     if t <= N else pyo.Constraint.Skip)
    # model.constraint42 = pyo.Constraint(model.tIDX, 
    #                                     rule=lambda model, t: forward_kinematics([model.x[i,t] for i in model.xIDX])[2][2] >= -1 - eps
    #                                     if t <= N else pyo.Constraint.Skip)

    # Path Constraint
    model.constraint43 = pyo.Constraint(model.tIDX,
                                        rule=lambda model, t: forward_kinematics_position([model.x[i,t] for i in model.xIDX])[1] <= 0.16701999999999667 + eps
                                        if t <= N else pyo.Constraint.Skip)
    model.constraint44 = pyo.Constraint(model.tIDX, 
                                        rule=lambda model, t: forward_kinematics_position([model.x[i,t] for i in model.xIDX])[1] >= 0.16701999999999667 - eps
                                        if t <= N else pyo.Constraint.Skip)
    model.constraint45 = pyo.Constraint(model.tIDX,
                                        rule=lambda model, t: forward_kinematics_position([model.x[i,t] for i in model.xIDX])[2] <= 0.5184984909676793 + eps
                                        if t <= N else pyo.Constraint.Skip)
    model.constraint46 = pyo.Constraint(model.tIDX, 
                                        rule=lambda model, t: forward_kinematics_position([model.x[i,t] for i in model.xIDX])[2] >= 0.5184984909676793 - eps
                                        if t <= N else pyo.Constraint.Skip)

    solver = pyo.SolverFactory('ipopt')
    results = solver.solve(model)
    
    if str(results.solver.termination_condition) == "optimal":
        feas = True
    else:
        feas = False
            
    xOpt = np.asarray([[model.x[i,t]() for i in model.xIDX] for t in model.tIDX]).T
    uOpt = np.asarray([model.u[:,t]() for t in model.tIDX]).T
    
    JOpt = model.cost()
    
    return [model, feas, xOpt, uOpt, JOpt]

Ts = 1
A = np.eye(7)
B = np.eye(7)*Ts
Q = np.eye(3)*100
R = np.eye(7)
P = Q
N = 10
# x0 = np.array([0.1, 0.2, -0.3, 0.4, 0.5, 0.6, 0.7]) #State Constraints Only
x0 = np.array([0, -1, 0, 1, 0, 1.6, 1.57079632679]) #Orientation Constraints or Path Constraints Only
# yref = [0.5, 0.5, 0.5] #State Constraints Only
# yref = [0.747, 0.312, 0.77] #Orientation Constraints Only
yref = [0.6935683264807795-0.5, 0.16701999999999667, 0.5184984909676793] #Path Constraints Only

U_lim = 0.2
# esp = 0.3 #Orientation Constraints Only
esp = 0.01 #Path Constraints Only

from ttictoc import tic,toc
tic()
[model, feas, xOpt, uOpt, JOpt] = solve_cftoc(A, B, P, Q, R, N, x0, yref, U_lim, esp)
elapsed = toc()
print('Elapsed time:',elapsed)

xOpt = np.array(xOpt)
uOpt = np.array(uOpt)

import kin_func_skeleton as kfs
yOpt = np.zeros((3,xOpt.shape[1]))
for i in range(xOpt.shape[1]):
    joint_angles = xOpt[:,i]
    gst = np.dot(kfs.prod_exp(np.array(twists), joint_angles), gst0)
    yOpt[:,i] = np.array([gst[0][3], gst[1][3], gst[2][3]]).T

# print('JOpt=', JOpt)
# print('xOpt=', xOpt)
# print('uOpt=', uOpt)
print('yOpt=', yOpt)

x0_tuck = [0, -1, 0, 1, 0, 1.6, 1.57079632679]
gst_tuck = np.dot(kfs.prod_exp(np.array(twists), np.array(x0_tuck)), gst0)
y_tuck = [gst_tuck[0][3], gst_tuck[1][3], gst_tuck[2][3]]
# print('y_tuck')
# print(y_tuck)

# Run with filename to take save data on xOpt
if len(sys.argv) == 2:
    import csv
    with open(sys.argv[1] + '_x' + '.csv', 'w') as f:
        mywriter = csv.writer(f, delimiter=',')
        mywriter.writerows(xOpt)
    with open(sys.argv[1] + '_input' + '.csv', 'w') as f:
        mywriter = csv.writer(f, delimiter=',')
        mywriter.writerows(uOpt)
    with open(sys.argv[1] + '_output' + '.csv', 'w') as f:
        mywriter = csv.writer(f, delimiter=',')
        mywriter.writerows(yOpt)
