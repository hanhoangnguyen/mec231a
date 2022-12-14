#!/usr/bin/env python
"""
MPC Controller for MEC231A
Author: Han Nguyen
"""

import numpy as np
import pyomo
import pyomo.environ as pyo
from pyomo.opt import SolverStatus, TerminationCondition

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


def solve_cftoc(A, B, P, Q, R, N, x0, yref):

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

    # Objective:
    def objective_rule(model):
        costY = 0.0
        costU = 0.0
        costTerminal = 0.0
        for t in model.tIDX:
            if t < model.N:
                joint_angles = [model.x[i,t] for i in model.xIDX]
                gst = matmul(prod_exp(twists,joint_angles), gst0)
                y = [gst[3][0], gst[3][0], gst[3][2]]
                for i in range(3):
                    for j in range(3):
                        costY += (y[i]-yref[i]) * model.Q[i,j] * (y[j]-yref[j]) 
        for t in model.tIDX:
            for i in model.uIDX:
                for j in model.uIDX:
                    if t < model.N:
                        costU += model.u[i, t] * model.R[i, j] * model.u[j, t]
        joint_angles_N = [model.x[i,model.N] for i in model.xIDX]
        gst_N = matmul(prod_exp(twists,joint_angles_N), gst0)
        y_N = [gst_N[3][0], gst_N[3][0], gst_N[3][2]]
        for i in range(3):
            for j in range(3):               
                costTerminal += (y_N[i]-yref[i]) * model.P[i,j] * (y_N[j]-yref[j])
        return costY #+ costU #+ costTerminal
    
    model.cost = pyo.Objective(rule = objective_rule, sense = pyo.minimize)
    
    # Constraints:
    def dynamic_constraint_rule(model, i, t):
        return model.x[i, t+1] - (sum(model.A[i, j] * model.x[j, t] for j in model.xIDX)
                            +  sum(model.B[i, j] * model.u[j, t] for j in model.uIDX) ) == 0.0 if t < model.N else pyo.Constraint.Skip

    model.equality_constraints = pyo.Constraint(model.xIDX, model.tIDX, rule=dynamic_constraint_rule)
    model.init_const1 = pyo.Constraint(expr = model.x[0, 0] == x0[0])
    model.init_const2 = pyo.Constraint(expr = model.x[1, 0] == x0[1])
    model.init_const3 = pyo.Constraint(expr = model.x[2, 0] == x0[2])
    model.init_const4 = pyo.Constraint(expr = model.x[3, 0] == x0[3])
    model.init_const5 = pyo.Constraint(expr = model.x[4, 0] == x0[4])
    model.init_const6 = pyo.Constraint(expr = model.x[5, 0] == x0[5])
    model.init_const7 = pyo.Constraint(expr = model.x[6, 0] == x0[6])
        
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

Ts = 0.1
A = np.eye(7)
B = np.eye(7)*Ts
Q = np.eye(3)*1000
R = np.eye(7)
P = Q
N = 1000
x0 = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
yref = [0.5, 0.5, 0.5]

[model, feas, xOpt, uOpt, JOpt] = solve_cftoc(A, B, P, Q, R, N, x0, yref)

# print('JOpt=', JOpt)
print('xOpt=', xOpt)
# print('uOpt=', uOpt)
# print('feas=', feas)

xOpt = np.array(xOpt)
import csv
#Change the filename for data taking!
with open('mpc_data3.csv', 'w') as f:
    mywriter = csv.writer(f, delimiter=',')
    mywriter.writerows(xOpt)