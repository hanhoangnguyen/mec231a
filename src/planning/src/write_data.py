#!/usr/bin/env python
from intera_interface import Limb
import rospy
import time
import tf2_ros
import csv
import numpy as np
import sys
def main():
    """
    Main Script
    """
    tfBuffer = tf2_ros.Buffer()
    tfListener = tf2_ros.TransformListener(tfBuffer)
    limb = Limb("right")
    joint_names = limb.joint_names()
    transform_joint_names = ['right_l0', 'right_l1', 'right_l2', 'right_l3', 'right_l4', 'right_l5', 'right_l6', 'right_hand']
    data = np.array([])
    while not rospy.is_shutdown():
        try:
            joint_angles = limb.joint_angles()
            joint_velocities = limb.joint_velocities()
            endpoint_pose = limb.endpoint_pose()
            endpoint_velocity = limb.endpoint_velocity()
            data_row = np.array([])
            for i in joint_names:
                data_row = np.append(data_row, joint_angles[i])
            for i in joint_names:
                data_row = np.append(data_row, joint_velocities[i])
            data_row = np.append(data_row, endpoint_pose['position'].x)
            data_row = np.append(data_row, endpoint_pose['position'].y)
            data_row = np.append(data_row, endpoint_pose['position'].z)
            data_row = np.append(data_row, endpoint_pose['orientation'].x)
            data_row = np.append(data_row, endpoint_pose['orientation'].y)
            data_row = np.append(data_row, endpoint_pose['orientation'].z)
            data_row = np.append(data_row, endpoint_pose['orientation'].w)
            data_row = np.append(data_row, endpoint_velocity['linear'].x)
            data_row = np.append(data_row, endpoint_velocity['linear'].y)
            data_row = np.append(data_row, endpoint_velocity['linear'].z)
            data_row = np.append(data_row, endpoint_velocity['angular'].x)
            data_row = np.append(data_row, endpoint_velocity['angular'].y)
            data_row = np.append(data_row, endpoint_velocity['angular'].z)
            for joint in transform_joint_names:
                trans = tfBuffer.lookup_transform('base',joint, rospy.Time())
                data_row = np.append(data_row, trans.transform.translation.x)
                data_row = np.append(data_row, trans.transform.translation.y)
                data_row = np.append(data_row, trans.transform.translation.z)
                data_row = np.append(data_row, trans.transform.rotation.x)
                data_row = np.append(data_row, trans.transform.rotation.y)
                data_row = np.append(data_row, trans.transform.rotation.z)
                data_row = np.append(data_row, trans.transform.rotation.w)
            data_row = np.append(data_row, time.time())
            if not data.any():
                data = data_row
            else:
                data = np.vstack((data, data_row))
        except (tf2_ros.LookupException, tf2_ros.ConnectivityException, tf2_ros.ExtrapolationException) as e:
            pass
    with open(sys.argv[1] + '.csv', 'w') as f:
        mywriter = csv.writer(f, delimiter=',')
        mywriter.writerows(data)

if __name__ == '__main__':
    rospy.init_node('writer_node')
    if len(sys.argv) != 2:
        print('Specify file name')
    else:
        main()
