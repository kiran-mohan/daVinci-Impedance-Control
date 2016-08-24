All script files in this package were developed in MATLAB R2012a and have been made compatible in MATLAB R2015a.

The DH parameters used here are taken directly from Intuitive Surgical's documentation. We did find some problems with this DH parameters and have corrected it
in the Impedance Control code package.

The MAIN.m file should be run to see the output of the package. The Inverse kinematics has been done for all 6 joints without decoupling the spherical joint.
The issue with this is that the orientation never converges. Hence we understood that we need to decouple the spherical joint to find Inverse Kinematics for
just three joints and then rotate by the Euler angles. However, we did not have the need nor the time to get this corrected as the Impedance Control model,
which is our primary scope does not need explicit computation of Inverse Kinematics.

The cubicTraj.m file generates a cubic trajectory for the End effector's cartesian position and orientation.

The dh2mat.m file converts a row of the DH table into a 4x4 transformation matrix.

The eul2rotm_ip_rad.m file is a function that converts the Euler angles to a Rotation matrix. The Euler angles need to be input in radians.

The rotm2eul.m file is a function that converts a rotation matrix into euler angles.

The plotarm_rad_6q.m file is a function used to plot the da Vinci PSM given a vector of joint angles.

Developers Contact
Amit Trivedi	    :	  atrivedi@wpi.edu
Terence Carmichael  :     twcarmichael@wpi.edu
Kiran Mohan         :     kmohan@wpi.edu
Aman Rana           :     arana@wpi.edu
Akanksha Devkar     :     acdevkar@wpi.edu