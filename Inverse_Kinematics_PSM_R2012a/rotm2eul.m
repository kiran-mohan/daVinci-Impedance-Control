%% Function to convert Rotation matrix to Euler angles
function eul = rotm2eul(T)
R = T(1:3,1:3);
if (abs(R(3,1)) ~= 1)
    theta1 = -asin(R(3,1));
%     theta2 = pi - theta1;
    si1 = atan2(R(3,2)/cos(theta1),R(3,3)/cos(theta1));
%     si2 = atan2(R(3,2)/cos(theta2),R(3,3)/cos(theta2));
    phi1 = atan2(R(2,1)/cos(theta1),R(1,1)/cos(theta1));
%     phi2 = atan2(R(2,1)/cos(theta2),R(1,1)/cos(theta2));
    theta = theta1;
    si = si1;
    phi = phi1;
else
    phi = 0;
    if (R(3,1) == -1)
        theta = pi/2;
        si = phi + atan2(R(1,2),R(1,3));
    else
        theta = -pi/2;
        si = -phi + atan2(-R(1,2),-R(1,3));
    end
end
eul = [si,theta,phi];
end
