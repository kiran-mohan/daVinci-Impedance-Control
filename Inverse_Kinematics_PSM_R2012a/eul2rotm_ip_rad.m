%% This function converts Euler angles to a Rotation Matrix.
% Assumption: We assume that the Z-Y-X order of rotation is used.
function [rotm] = eul2rotm(s,t,p)
    A = [1 0 0;0 cos(s) -sin(s); 0 sin(s) cos(s)];
    B = [cos(t) 0 sin(t);0 1 0;-sin(t) 0 cos(t)];
    C = [cos(p) -sin(p) 0;sin(p) cos(p) 0; 0 0 1];
    rotm = C*B*A;
end