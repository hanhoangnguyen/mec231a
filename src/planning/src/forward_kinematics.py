#!/usr/bin/env python
import rospy
import numpy as np
import kin_func_skeleton as kfs
import tf2_ros
import tf
from scipy.spatial.transform import Rotation
from intera_interface import Limb

def main():
    """
    Computes the orientation of the Sawyer's left end-effector given the joint
    angles in radians. Record the initial twists in the zero configuration first
    before moving the arm.
    """
    limb = Limb("right")
    qs = np.ndarray((3,8)) # points on each joint axis in the zero configuration
    ws = np.ndarray((3,7)) # axis vector of each joint axis
    R = None

    tfBuffer = tf2_ros.Buffer()
    tfListener = tf2_ros.TransformListener(tfBuffer)
    transform_joint_names = ['right_l0', 'right_l1', 'right_l2', 'right_l3', 'right_l4', 'right_l5', 'right_l6', 'right_hand']
    input("Go to Zero Configuration Using 'rosrun intera_examples go_to_joint_angles.py -q 0 0 0 0 0 0 0. Press a key when done.")
    recorded_zero_configuration = False
    while not rospy.is_shutdown():
        try:
            if recorded_zero_configuration == False:
                for index,joint in enumerate(transform_joint_names):
                    trans = tfBuffer.lookup_transform('base',joint, rospy.Time())
                    qs[0:3,index] = [trans.transform.translation.x, trans.transform.translation.y, trans.transform.translation.z]
                    r = Rotation.from_quat([trans.transform.rotation.x, trans.transform.rotation.y, trans.transform.rotation.z, trans.transform.rotation.w])
                    if index < 7:
                        ws[0:3, index] = r.as_dcm()[0:3,2] #Assumption: Axis of rotation is the Z vector in rotation matrix
                    elif index == 7:
                        R = np.array(r.as_dcm())

                    vs = np.ndarray((3,7))
                    twists = np.ndarray((6,7))

                    for i in range(vs.shape[1]):
                        vs[0:3,i] = -np.dot(kfs.skew_3d(ws[:, i]), qs[:, i])
                    
                    twists[0:3,0:7] = vs
                    twists[3:6,0:7] = ws

                    gst0 = np.eye(4)
                    gst0[0:3,0:3] = R
                    gst0[0:3,3] = qs[:,7]
                recorded_zero_configuration = True
                print("Finished recording zero configuration.")
            joint_angles = np.array(list(limb.joint_angles().values()))
            gst = np.dot(kfs.prod_exp(twists, joint_angles), gst0)
            print("Gst: ", gst)
            pose = tfBuffer.lookup_transform('base','right_hand', rospy.Time())
            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion([pose.transform.rotation.x, pose.transform.rotation.y,
                                                                           pose.transform.rotation.z, pose.transform.rotation.w])
            # translation = pose.transform.translation
            # print("Euler Transform: ", translation.x, translation.y, translation.z)
            # print('Euler Matrix: ', tf.transformations.euler_matrix(roll, pitch, yaw))
            # print()
        except (tf2_ros.LookupException, tf2_ros.ConnectivityException, tf2_ros.ExtrapolationException) as e:
            pass

if __name__ == '__main__':
    rospy.init_node('sawyer_fk_node')
    main()
