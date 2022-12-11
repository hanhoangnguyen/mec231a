#!/usr/bin/env python
"""
Controller Class for Lab 8
Author: Valmik Prabhu, Chris Correa
"""

import rospy
import sys
import numpy as np
import itertools
from collections import deque
import pyomo
import baxter_interface
import intera_interface
import pyomo.environ as pyo
from pyomo.opt import SolverStatus, TerminationCondition
from moveit_msgs.msg import RobotTrajectory



class mpc_controller(object):
    """
    A Model predictive controller object for EE220B

    Fields:
    


    """

    def __init__(self, Kp, Ki, Kd, Kw, limb):
        """
        Constructor:

        Inputs:
        Kp: 7x' ndarray of proportional constants
        Ki: 7x' ndarray of integral constants
        Kd: 7x' ndarray of derivative constants
        Kw: 7x' ndarray of antiwindup constants
        limb: sawyer_interface.Limb
        """

        # If the node is shutdown, call this function
        rospy.on_shutdown(self.shutdown)

        # self._Kp = Kp
        # self._Ki = Ki
        # self._Kd = Kd
        # self._Kw = Kw

        # self._LastError = np.zeros(len(Kd))
        # self._LastTime = 0;
        # self._IntError = np.zeros(len(Ki))
        # self._ring_buff_capacity = 3
        # self._ring_buff = deque([], self._ring_buff_capacity)

        self._path = RobotTrajectory()
        # self._curIndex = 0;
        # self._maxIndex = 0;

        self._limb = limb

        # For Plotting:
        self._times = list()
        self._actual_positions = list()
        self._actual_velocities = list()
        self._target_positions = list()
        self._target_velocities = list()

    
    def shutdown(self):
        """
        Code to run on shutdown. This is good practice for safety
        """
        rospy.loginfo("Stopping Controller")

        # Set velocities to zero
        dic_vel = {name:np.zeros(len(self._limb.joint_names())) for name in self._limb.joint_names()}

        self._limb.set_joint_velocities(dic_vel)
        rospy.sleep(0.1)

    

    def execute_plan(self, path, timeout=100.0, log=True):
        """
        Execute a given path

        Inputs:
        path: a moveit_msgs/RobotTrajectory message
        timeout: max time the controller will run
        log: should the controller display a plot
        
        """


        '''


        if log:
            import matplotlib.pyplot as plt

            times = np.array(self._times)
            actual_positions = np.array(self._actual_positions)
            actual_velocities = np.array(self._actual_velocities)
            target_positions = np.array(self._target_positions)
            target_velocities = np.array(self._target_velocities)
            plt.figure()
            joint_num = len(self._path.joint_trajectory.joint_names)
            for joint in range(joint_num):
                plt.subplot(joint_num,2,2*joint+1)
                plt.plot(times, actual_positions[:,joint], label='Actual')
                plt.plot(times, target_positions[:,joint], label='Desired')
                plt.xlabel("Time (t)")
                if(joint == 0):
                    plt.ylabel(self._path.joint_trajectory.joint_names[joint] + " Position Error")
                else:
                    plt.ylabel(self._path.joint_trajectory.joint_names[joint])
                plt.legend()

                plt.subplot(joint_num,2,2*joint+2)
                plt.plot(times, actual_velocities[:,joint], label='Actual')
                plt.plot(times, target_velocities[:,joint], label='Desired')
                plt.xlabel("Time (t)")
                if(joint == 0):
                    plt.ylabel(self._path.joint_trajectory.joint_names[joint] + " Velocity Error")
                else:
                    plt.ylabel(self._path.joint_trajectory.joint_names[joint])
                plt.legend()

            print("Close the plot window to continue")
            plt.show()
        '''
        return True

    def step_control(self, t):
        """
        Return the control input given the current controller state at time t

        Inputs:
        t: time from start in seconds

        Output:
        u: 7x' ndarray of velocity commands
        
        """
       
    
    # Setup the CTFOC problem using Pyomo - I don't know if we will acc need this 
    def cftoc(x0,N,input_bounds,vref,Xf={'type':''},p={'type':''}):     
    # x is state vector, u is input, Ts is sampling period.   
        model = pyo.ConcreteModel()
        model.tidx = pyo.Set(initialize=range(0, N+1)) # length of finite optimization problem
        model.tidu = pyo.Set(initialize=range(0, N)) # length of finite optimization problem
        model.xidx = pyo.Set(initialize=range(0, num_states)) #number of states I think this will be 6 but unsure
        model.uidx = pyo.Set(initialize=range(0, num_inputs)) #not sure what our inputs are going to be either

        # Create state and input variables trajectory:
        model.x = pyo.Var(model.xidx, model.tidx)
        model.u = pyo.Var(model.uidx, model.tidu)

        #We are probably going to need to edit this, not sure what our cost function will be 
        model.R=2
        model.Q=1
    
        #This code handles 3 types of terminal constraint: cvxhull, full-dimentional polytope and terminal equalty constraint
        #I think our terminal constraint is going to the be the last one polytope_eq
            if  Xf['type']=='cvxhull':
            model.nf = np.size(Xf['SS'], 0)
            model.nfidx = pyo.Set( initialize= range(model.nf), ordered=True )
            print( model.nf )
            model.SS = Xf['SS']
            model.lambdavar = pyo.Var(model.nfidx)
        
        if Xf['type']=='polytope':
            model.nf = np.size(Xf['Af'], 0)
            model.nfidx= pyo.Set( initialize= range(model.nf), ordered=True )
            model.Af = Xf['Af']
            model.bf = Xf['bf']

        if Xf['type']=='polytope_eq':
            model.nf = np.size(Xf['Aeq'], 0)
            model.nfidx= pyo.Set( initialize= range(model.nf), ordered=True )
            model.Aeq = Xf['Aeq']
            model.beq = Xf['beq']

        #This code handles 2 types of terminal const : cvxhull, and quadratic cost
        #We definitely want to go with quadratic cost
        if  p['type']=='cvxhull':
            model.cvalue=p['patsamples']
        
        if  p['type']=='quadratic':
            model.P=p['P']
        
        # Constraints:
        #Initial condition - I think this is fine 
        model.constraint1 = pyo.Constraint(model.xidx, rule=lambda model, i: model.x[i, 0] == x0[i])

        #State dynamics - These will definitely need to be edited
        # model.constraint2 = pyo.Constraint(model.tidx, rule=lambda model, t: model.x[0, t+1] == model.x[0, t] + Ts*model.x[1, t]
        #                                     if t < N else pyo.Constraint.Skip)
        # model.constraint3 = pyo.Constraint(model.tidx, rule=lambda model, t: model.x[1, t+1] == model.x[1, t] + Ts*model.u[0,t]
        #                                     if t < N else pyo.Constraint.Skip)

        #State bounds
        # model.constraint4 = pyo.Constraint(model.tidx, rule=lambda model, t: model.x[1, t] >= 0) # non-strict inequalities not allowed
        # model.constraint5 = pyo.Constraint(model.tidx, rule=lambda model, t: model.x[1, t] <= v_max) # Speed limit

        #Input bounds
        # model.constraint6 = pyo.Constraint(model.tidu, rule=lambda model, t: model.u[0,t] >= -g*mu)
        # model.constraint7 = pyo.Constraint(model.tidu, rule=lambda model, t: model.u[0,t] <= g*mu)    

        #State constraints
        #model.constraint8 = pyo.Constraint(model.tidx, rule = lambda model, t: model.x[0,t] <= d_bike - w_bike/2-  l_car/2 )              

            
        # if  Xf['type']=='cvxhull':
        #     model.constraint9= pyo.Constraint(model.xidx, rule=lambda model, k: model.x[k,N] == sum( (model.SS[i,k] * model.lambdavar[i]) for i in model.nfidx))
        #     model.constraint10= pyo.Constraint(rule=lambda model: sum(model.lambdavar[i]  for i in model.nfidx)  == 1 )
        #     model.constraint11= pyo.Constraint(model.nfidx, rule=lambda model, i:   model.lambdavar[i]>=0)
        # if  Xf['type']=='polytope':
        #     def final_const_rule(model, i):
        #         return sum(model.Af[i, j] * model.x[j, N] for j in model.xidx) <= model.bf[i] 
        #     model.final_const = pyo.Constraint(model.nfidx, rule=final_const_rule)
        if  Xf['type']=='polytope_eq':
            def final_const_rule(model, i):
                return sum(model.Aeq[i, j] * model.x[j, N] for j in model.xidx) == model.beq[i] 
            model.final_const = pyo.Constraint(model.nfidx, rule=final_const_rule)

        if p['type']=='cvxhull':
            #jerk minimiation
            #model.cost = pyo.Objective(expr = sum((model.u[i, t+1]-model.u[i, t])**2 for i in model.uidx for t in model.tidx if t < N-1) + sum((model.x[1, t]-vref)**2 for t in model.tidx), sense=pyo.minimize)
            #acceleration minimiation
            model.cost = pyo.Objective(expr = (sum(model.R*(model.u[0, t])**2  for t in model.tidu ) + sum(model.Q*(model.x[1, t]-vref)**2 for t in model.tidx) +sum((model.cvalue[i] * model.lambdavar[i]) for i in model.nfidx)), sense=pyo.minimize)
        elif  p['type']=='quadratic':
            model.cost = pyo.Objective(expr = (sum(model.R*(model.u[0, t])**2  for t in model.tidu ) + sum(model.Q*(model.x[1, t]-vref)**2 for t in model.tidu) +model.P*(model.x[1, N]-vref)**2), sense=pyo.minimize)
        else:
            model.cost = pyo.Objective(expr = (sum(model.R*(model.u[0, t])**2  for t in model.tidu ) + sum(model.Q*(model.x[1, t]-vref)**2 for t in model.tidu)), sense=pyo.minimize)
        # Now we can solve:
        results = pyo.SolverFactory('ipopt').solve(model)
        d = pyo.value(model.x[0,:])
        v = pyo.value(model.x[1,:])
        a = pyo.value(model.u[0,:])
        c = pyo.value(model.cost)

        return d, v, a, results.solver.termination_condition,c
    

    #this is Hansung's code for solving the cftoc and I think this will be easier to formulate
    def solve_cftoc(A = 6, B = 6, P, Q, R, N, x0, xL, xU, uL, uU, bf=np.nan, Af=np.nan):
    
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
        model.Q = Q #Not sure what exactly our P Q and R r gonna be 
        model.P = P
        model.R = R
        
        # Create state and input variables trajectory:
        model.x = pyo.Var(model.xIDX, model.tIDX, bounds=(xL,xU))
        model.u = pyo.Var(model.uIDX, model.tIDX, bounds=(uL,uU))
        
        #Objective:
        def objective_rule(model):
            costX = 0.0
            costU = 0.0
            costTerminal = 0.0
            for t in model.tIDX:
                for i in model.xIDX:
                    for j in model.xIDX:
                        if t < model.N:
                            costX += model.x[i, t] * model.Q[i, j] * model.x[j, t] 
            for t in model.tIDX:
                for i in model.uIDX:
                    for j in model.uIDX:
                        if t < model.N:
                            costU += model.u[i, t] * model.R[i, j] * model.u[j, t]
            for i in model.xIDX:
                for j in model.xIDX:               
                    costTerminal += model.x[i, model.N] * model.P[i, j] * model.x[j, model.N]
            return costX + costU + costTerminal
        
        model.cost = pyo.Objective(rule = objective_rule, sense = pyo.minimize)
        
        # Constraints:
        def equality_const_rule(model, i, t):
            return model.x[i, t+1] - (sum(model.A[i, j] * model.x[j, t] for j in model.xIDX)
                                +  sum(model.B[i, j] * model.u[j, t] for j in model.uIDX) ) == 0.0 if t < model.N else pyo.Constraint.Skip

        model.equality_constraints = pyo.Constraint(model.xIDX, model.tIDX, rule=equality_const_rule)
        model.init_const1 = pyo.Constraint(expr = model.x[0, 0] == x0[0])
        model.init_const2 = pyo.Constraint(expr = model.x[1, 0] == x0[1])
            
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


    #main MPC control law generation:
    #might have to re pass in parameters as it is no longer a notebook
    def generate_uOpt():
        nx = np.size(A, 0)
        nu = np.size(B, 1)

        M = 25   # Simulation time

        xOpt = np.zeros((nx, M+1))
        uOpt = np.zeros((nu, M))
        xOpt[:, 0] = x0.reshape(nx, )

        xPred = np.zeros((nx, N+1, M))

        feas = np.zeros((M, ), dtype=bool)

        fig = plt.figure(figsize=(9, 6))
        for t in range(M):
            [model, feas[t], x, u, J] = solve_cftoc(A, B, P, Q, R, N, xOpt[:, t], xL, xU, uL, uU)
            
            if not feas[t]:
                xOpt = []
                uOpt = []
                break
            
            xOpt[:, t+1] = x[:, 1]
            uOpt[:, t] = u[:, 0].reshape(nu, )

    #sorry i had to leave it here becuase I was going to miss the bus. 

if __name__ == '__main__': 
    pass



    





