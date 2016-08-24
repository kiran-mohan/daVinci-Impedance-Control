%% FK for 6 links and IK for 3 links of PSM

%%%%%%% MATLAB's Symbolic Toolbox has been used in this code %%%%%%%

% This code was developed in MATLAB R2012a and has been tested in 
% MATLAB R2015a. Although we do not guarantee it, the code in this file 
% should be backward compatible.

% This file contains MATLAB code for determining the forward and inverse 
% kinematics for all 6 links of the Patient Side Manipulator (PSM) and
% Joint trajectory tracking of all joints.
% However, the inverse kinematics for joints 4,5 and 6 are still
% Work-In-Progress. There are some issues with the computation causing the 
% error to oscillate and never converge.
% Forward kinematics is implemented using DH parameters
% Inverse kinematics is implemented using the Pseudo-inverse-Jacobian
% technique which is based on an iterative solution (numerical methods).
% Change desired position of end effector by changing variables xdx, xdy & 
% xdz to check for different configurations of the robot.
%% Initialization
clc;clear all;close all;
syms q1 q2 q3 q4 q5 q6 real;
syms t1(t) t2(t) t3(t) t4(t) t5(t) t6(t) t real;
xdx = 0.5; % meters
xdy = -0.6; % meters
xdz = 0.2; % meters
xdpsi = 0;
xdtheta = 0;
xdphi = 0;
lRcc = 0.4318; % meters
lTool = 0.4162; % meters
lP2Y = 0.0091; % meters
lY2CP = 0.0102; % meters
tstep = 0.005;
t = 0:tstep:1;
jntLim = 2.6179; % Joint limits
posAcc = 1e-3; % = 1 millimeter position accuracy
oriAcc = 0.1; % = 0.1 radians angular orientation accuracy
% Joint tracking initializations
% All joints are initialized to zeros for variable consistency
jnt1tr = zeros(3,length(t));
jnt2tr = zeros(3,length(t));
jnt3tr = zeros(3,length(t));
jnt4tr = zeros(3,length(t));
jnt5tr = zeros(3,length(t));
jnt6tr = zeros(3,length(t));
jnt7tr = zeros(3,length(t));

%% DH Parameters
DH = [q1+pi/2,     0,         pi/2,      0;...
      q2-pi/2,     0,         -pi/2,     0;...
      0,           q3-lRcc,   pi/2,      0;...
      q4,          lTool,     0,         0;...
      q5-pi/2,     0,         -pi/2,     0;...
      q6-pi/2,     0,         -pi/2,     lP2Y;...
      0,           lY2CP,     -pi/2,     0];

%% Forward Kinematics
% Final transformation matrix (from base to tip) initialized as 
% identity matrix
T_fin = eye(4);

% Loop run for each row in DH table. T is a 3 dimensional tensor which
% stores the transformation matrices T01 T12 T23 and so on. Each matrix is
% stored in an increment in the third dimension of T. The final
% transformation matrix from base to tip is stored in T_fin by multiplying
% the matrices T01*T12*T23*...*T67. The matrix Ti has the position vectors
% P01 P12 P23...P67 stored in columns.
for i = 1:size(DH,1)
    T(:,:,i) = dh2mat(DH(i,1),DH(i,2),DH(i,3),DH(i,4));
    T_fin = T_fin*T(:,:,i);
    Ti(:,i) = T(1:3,4);
end
% Obtaining intermediate transformations with respect to base.
T01 = T(:,:,1);
T02 = T01*T(:,:,2);
T03 = T02*T(:,:,3);
T04 = T03*T(:,:,4);
T05 = T04*T(:,:,5);
T06 = T05*T(:,:,6);
T07 = T06*T(:,:,7);
P_fin = T_fin(1:3,4); %Extracting only the x-y-z position from T matrix

%% Home configuration (0,0,0,0,0,0)
THome = subs(T_fin,[q1 q2 q3 q4 q5 q6],[0 0 0 0 0 0]);
RHome = THome(1:3,1:3);
PHome = THome(1:3,4);
    
%% Calculate symbolic Jacobian
Jsym = simplify(jacobian(P_fin,[q1,q2,q3,q4,q5,q6])); 
omega = [ 0  -1   0   0  -1   0;...
         -1   0   0   0   0  -1;...
          0   0  -1  -1   0   0];
eul = rotm2eul(THome);
phi = eul(1);theta = eul(2);psi = eul(3);
om2eul = [0 -sin(phi) cos(phi)*cos(theta);...
          0 cos(phi) cos(theta)*sin(phi);...    
          1 0 -sin(theta)];
jac_an_ang = om2eul\omega;
Jsym = [Jsym;jac_an_ang];

%% Pseudo-Inverse-Jacobian based Inverse Kinematics
chk=0;
% Desired position and orientation
xdx_tr = cubicTraj(PHome(1),0,xdx,0,2,1);
xdy_tr = cubicTraj(PHome(2),0,xdy,0,2,1);
xdz_tr = cubicTraj(PHome(3),0,xdz,0,2,1);
xdpsi_tr = cubicTraj(psi,0,xdpsi,0,2,1);
xdtheta_tr = cubicTraj(theta,0,xdtheta,0,2,1);
xdphi_tr = cubicTraj(phi,0,xdphi,0,2,1);

xd = [xdx;xdy;xdz;rotm2eul([1 0 0;0 1 0;0 0 1])'];
disp('Desired position:');
disp(xd);
% Robot current position and orientation
x = [PHome;rotm2eul(RHome)'];
% qi has the updated joint values in each iteration
% Initially set to zero
qi = zeros(6,1);
% q_n has the joint data from qi gathered over all iterations
q_n = zeros(6,length(t));
disp('Processing...');
% Set learning rate to 0.1
alpha = 0.1;

time_now = 0:0.01:2;

% Iterate until error is reduced, max iterations is reached or goal state
% is unreachable.
for idx = 1:length(t)
    % Initial plot (Home config)
    if(idx==1)
        plotarm_rad_6q(qi,idx);
        xlabel('X axis (m)');
        ylabel('Y axis (m)');
        zlabel('Z axis (m)');
    end
    % Imposing joint limits
    qchk = double([qi(1:3);qi(5:6)]);
    if(sum(abs(qchk)>jntLim)>0)
        chk=1;
        break; % discontinue loop if joint limit is exceeded
    else
        % Calculate the updated Jacobian
        T_loop = subs(T_fin,[q1 q2 q3 q4 q5 q6],qi');
        rot = double(T_loop(1:3,1:3));
        x = [T_loop(1:3,4);rotm2eul(rot)'];
        
        eul = rotm2eul(T_loop);
        phi = eul(1);theta = eul(2);chi = eul(3);
        om2eul = [0 -sin(phi) cos(phi)*cos(theta);...
          0 cos(phi) cos(theta)*sin(phi);...
          1 0 -sin(theta)];
        jac_an_ang = om2eul\omega;
        J = double(subs(Jsym,[q1,q2,q3,q4,q5,q6],qi'));
        J(4:6,:) = jac_an_ang;
        % Calculate the updated position and orientation of the robot
% % %         x = subs(T_fin,[q1 q2 q3 q4 q5 q6],qi');
        % Slice to separate orientation info into rot
% % %         rot = double(x(1:3,1:3));
        % Slice to separate position info into x
% % %         x = [x(1:3,4);rotm2eul(rot)'];
        % Transformation Matrices for each iteration
        % For joint trajectory tracking, we calculate the intermediate 
        % transformation matrices to obtain the position of each joint in 
        % task space during each iteration
        j1 = subs(T01,[q1,q2,q3,q4,q5,q6],qi');
        j2 = subs(T02,[q1,q2,q3,q4,q5,q6],qi');
        j3 = subs(T03,[q1,q2,q3,q4,q5,q6],qi');
        j4 = subs(T04,[q1,q2,q3,q4,q5,q6],qi');
        j5 = subs(T05,[q1,q2,q3,q4,q5,q6],qi');
        j6 = subs(T06,[q1,q2,q3,q4,q5,q6],qi');
        j7 = subs(T07,[q1,q2,q3,q4,q5,q6],qi');
        % Slice to obtain joint positions (task space) in each iteration
        jnt1tr(:,idx) = j1(1:3,4);
        jnt2tr(:,idx) = j2(1:3,4);
        jnt3tr(:,idx) = j3(1:3,4);
        jnt4tr(:,idx) = j4(1:3,4);
        jnt5tr(:,idx) = j5(1:3,4);
        jnt6tr(:,idx) = j6(1:3,4);
        jnt7tr(:,idx) = j7(1:3,4);
        % Difference between desired pose and current pose
        tim = time_now(idx);
        xdx_tr_now = xdx_tr(1)+xdx_tr(2)*tim+xdx_tr(3)*tim^2+xdx_tr(4)*tim^3;
        xdy_tr_now = xdy_tr(1)+xdy_tr(2)*tim+xdy_tr(3)*tim^2+xdy_tr(4)*tim^3;
        xdz_tr_now = xdz_tr(1)+xdz_tr(2)*tim+xdz_tr(3)*tim^2+xdz_tr(4)*tim^3;
        xdpsi_tr_now = xdpsi_tr(1)+xdpsi_tr(2)*tim+xdpsi_tr(3)*tim^2+xdpsi_tr(4)*tim^3;
        xdtheta_tr_now = xdtheta_tr(1)+xdtheta_tr(2)*tim+xdtheta_tr(3)*tim^2+xdtheta_tr(4)*tim^3;
        xdphi_tr_now = xdphi_tr(1)+xdphi_tr(2)*tim+xdphi_tr(3)*tim^2+xdphi_tr(4)*tim^3;
        
        xd = [xdx_tr_now;xdy_tr_now;xdz_tr_now;xdpsi_tr_now;xdtheta_tr_now;xdphi_tr_now];
        e = double(xd-x);
        err_hist(:,idx) = e;
        % Slice to obtain only positional difference
        ePos = e(1:3);
        % Slice to obtain only orientational difference
        eOri = e(4:6);
        % Combine data from eOri over all iterations in eOriMat
        eOriMat(:,idx) = eOri;
        % Combine data from ePos over all iterations in ePosMat
        ePosMat(:,idx) = ePos;
        % Plotting joint positions of robot
        hold on;
        [~,h,qx,qy,qz] = plotarm_rad_6q(qi,idx);
        xlabel('X axis (m)');
        ylabel('Y axis (m)');
        zlabel('Z axis (m)');
        % Delete previous plots if max_iterations not reached and if 
        % position or orientation accuracy not reached.
        % This is done for the plot to look like an animation
        if(idx<length(t) && (sum(abs(ePos)<posAcc)<3 || sum(abs(eOri)<oriAcc)<3))
            delete(h);
            delete(qx);
            delete(qy);
            delete(qz);
        end
        % Break the loop if positional and orientational accuracy has been 
        % reached.
%         if(sum(abs(ePos)<posAcc) == 3 && sum(abs(eOri)<oriAcc) == 3)
%             break;
%         end
        % Null space projection vector phi is set to reach home
        % configuration if robot reaches singularity.
        phi = [PHome;rotm2eul(RHome)'];
        % Calculate Pseudo inverse of Jacobian
        dagJ = J'/(J*J');
        % Store joint values over all iterations
        q_n(:,idx) = qi;
        % Calculate small change in joint angles (delq) required to produce 
        % a small change in error (e).
        delq = alpha*(dagJ*e + (eye(size(J,2))-dagJ*J)*phi);
        % Update qi by adding delq to the previous value of qi
        qi = qi + delq;
        % This is possible because we assume the system to be linear for
        % small changes in position of the end effector.
    end
end
disp('Final End effector position:');
disp(x);
disp('Joint values (in radians) for this position:');
disp(q_n(:,idx-1));
if(chk==1) % if joint limits violated
    disp('GOAL STATE NOT REACHABLE.!!!');
    disp(['The following joints exceeded joint limits after ',num2str(idx),' iterations']);
    tmp = find(abs(qi)>jntLim==1);
    for i = 1:length(tmp)
        disp(['Joint ',num2str(tmp(i))]);
    end
elseif (idx==length(t)) % if failed to converge within max_iterations
    disp(['FAILED TO CONVERGE WITHIN ',num2str(length(t)),' ITERATIONS']);
    % Plot orientational and positional error
    figure;
    subplot(2,1,1);
    plot(1:length(eOriMat),eOriMat);
    title('Orientation Error vs iterations');
    xlabel('Iterations (number)');
    ylabel('Orientation Error (rad)');
    legend('Joint 4','Joint 5','Joint 6','location','northeastoutside');
    subplot(2,1,2);
    plot(1:length(ePosMat),ePosMat);
    title('Position Error vs time');
    xlabel('Iterations (number)');
    ylabel('Position Error (m)');
    legend('Joint 1','Joint 2','Joint 3','location','northeastoutside');
else
    disp(['REACHED GOAL SUCCESSFULLY IN ',num2str(idx),' ITERATIONS.']);
end
% Truncate joint tracking variables upto converged iterations
jnt1tr = jnt1tr(:,1:idx);
jnt2tr = jnt2tr(:,1:idx);
jnt3tr = jnt3tr(:,1:idx);
jnt4tr = jnt4tr(:,1:idx);
jnt5tr = jnt5tr(:,1:idx);
jnt6tr = jnt6tr(:,1:idx);
jnt7tr = jnt7tr(:,1:idx);

% Plotting tracked joint positions
% only joints 3 to 7 are plotted here because Joints 1 and 2 are at base
% and do not move. 
figure;
title('Joint Position tracking');
plot3(jnt3tr(1,1:idx),jnt3tr(2,1:idx),jnt3tr(3,1:idx),'b');
xlabel('X axis (m)');
ylabel('Y axis (m)');
zlabel('Z axis (m)');
hold on;
plot3(jnt4tr(1,1:idx),jnt4tr(2,1:idx),jnt4tr(3,1:idx),'r');
xlabel('X axis (m)');
ylabel('Y axis (m)');
zlabel('Z axis (m)');
plot3(jnt5tr(1,1:idx),jnt5tr(2,1:idx),jnt5tr(3,1:idx),'g');
xlabel('X axis (m)');
ylabel('Y axis (m)');
zlabel('Z axis (m)');
% origins of joints 4 and 5 are at the same cartesian point and hence the
% trajectories followed by both joints will be the same.
plot3(jnt6tr(1,1:idx),jnt6tr(2,1:idx),jnt6tr(3,1:idx),'m');
xlabel('X axis (m)');
ylabel('Y axis (m)');
zlabel('Z axis (m)');
plot3(jnt7tr(1,1:idx),jnt7tr(2,1:idx),jnt7tr(3,1:idx),'c');
xlabel('X axis (m)');
ylabel('Y axis (m)');
zlabel('Z axis (m)');
legend('Joint 3','Joint 4','Joint 5','Joint 6','Joint 7');
title('Joint Position tracking');
hold off;
% Plot joint values q1 through q6
figure;
plot(1:idx-1,q_n(1:end,1:idx-1));
title('Joint angles');
legend('q1','q2','q3','q4','q5','q6','location','northeastoutside');
xlabel('Iterations (number)');
ylabel('Joint angles (rad)');