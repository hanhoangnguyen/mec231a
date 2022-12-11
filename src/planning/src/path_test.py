#!/usr/bin/env python
"""
Path Planning Script for Lab 7
Author: Valmik Prabhu
"""
import sys
from intera_interface import Limb
import rospy
import numpy as np
import traceback

from moveit_msgs.msg import OrientationConstraint
from geometry_msgs.msg import PoseStamped

from path_planner import PathPlanner

try:
    from controller import Controller
except ImportError:
    pass
    
def main():
    """
    Main Script
    """

    # Make sure that you've looked at and understand path_planner.py before starting

    planner = PathPlanner("right_arm")


    Kp = 0.2 * np.array([0.4, 2, 1.7, 1.5, 2, 2, 3])
    Kd = 0.01 * np.array([2, 1, 2, 0.5, 0.8, 0.8, 0.8])
    Ki = 0.01 * np.array([1.4, 1.4, 1.4, 1, 0.6, 0.6, 0.6])
    Kw = np.array([0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9])

    # SOLUTION -- hiden checkpoint 3
    limb = Limb("right")
    # controller = Controller(Kp, Ki, Kd, Kw, limb) # Solution checkpoint 3



    # # Given commented
    # # Add the obstacle to the planning scene here
    # #

    # #Create a path constraint for the arm
    # #UNCOMMENT FOR THE ORIENTATION CONSTRAINTS PART
    # orien_const = OrientationConstraint()
    # orien_const.link_name = "right_hand";
    # orien_const.header.frame_id = "base";
    # orien_const.orientation.y = 1.0;
    # orien_const.absolute_x_axis_tolerance = 0.5;
    # orien_const.absolute_y_axis_tolerance = 0.5;
    # orien_const.absolute_z_axis_tolerance = 0.5;
    # orien_const.weight = 1.0;



    # # SOLUTION ----- hidden checkpoint 2
    # # Add Table obstacle as a Box
    # """
    # Adds a rectangular prism obstacle to the planning scene

    # Inputs:
    # size: 3x' ndarray; (x, y, z) size of the box (in the box's body frame)
    # name: unique name of the obstacle (used for adding and removing)
    # pose: geometry_msgs/PoseStamped object for the CoM of the box in relation to some frame
    # """
    # size = np.array([0.4, 1.2, 0.1])
    # name = "Table1" 
    # pose = PoseStamped()
    # pose.header.frame_id = "base"

    # #x, y, and z position
    # pose.pose.position.x = 0.5
    # pose.pose.position.y = 0.0
    # pose.pose.position.z = 0.0

    # #Orientation as a quaternion
    # pose.pose.orientation.x = 0.0
    # pose.pose.orientation.y = 0.0
    # pose.pose.orientation.z = 0.0
    # pose.pose.orientation.w = 1.0   
    # planner.add_box_obstacle(size, name, pose)
    # # END SOLUTION 

    while not rospy.is_shutdown():

        while not rospy.is_shutdown():
            try:
                x, y, z = 0.8, 0.05, 0.07
                goal_1 = PoseStamped()
                goal_1.header.frame_id = "base"

                #x, y, and z position
                goal_1.pose.position.x = x
                goal_1.pose.position.y = y
                goal_1.pose.position.z = z

                #Orientation as a quaternion
                goal_1.pose.orientation.x = 0.0
                goal_1.pose.orientation.y = 1.0
                goal_1.pose.orientation.z = 0.0
                goal_1.pose.orientation.w = 0.0

                # Might have to edit this . . . 
                plan = planner.plan_to_pose(goal_1, [orien_const]) #solution
                # plan = planner.plan_to_pose(goal_1, [])
                input("Press <Enter> to move the right arm to goal pose 1: ")
                if not planner.execute_plan(plan[1]): 
                # if not controller.execute_plan(plan[1]): # Solution checkpoint 3
                    raise Exception("Execution failed")
            except Exception as e:
                print(e)
                traceback.print_exc()
            else:
                break

if __name__ == '__main__':
    rospy.init_node('moveit_node')
    main()
