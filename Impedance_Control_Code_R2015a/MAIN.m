% This script was written in MATLAB R2012a
% This is the main file in the package. This file should be run to get
% plots and results as in the project report. Forces can be set in the
% interactive UI enabled plots. However, large forces are not recommended
% as the joint limits may be exceeded by the control input. The mass of the
% system is low and hence, it can tolerate only small forces and even
% smaller moments.
%% Initialization of variables
close all;clc;
syms q_vec q1 q2 q3 q4 q5 q6 g real;
syms t t1(t) t2(t) t3(t) t4(t) t5(t) t6(t);
syms q1dot q2dot q3dot q4dot q5dot q6dot real;
syms t1dot(t) t2dot(t) t3dot(t) t4dot(t) t5dot(t) t6dot(t);
syms q1dotdot q2dotdot q3dotdot q4dotdot q5dotdot q6dotdot real;
syms t1dotdot(t) t2dotdot(t) t3dotdot(t) t4dotdot(t) t5dotdot(t) t6dotdot(t);
syms g real;
q_vec = [q1 q2 q3 q4 q5 q6].';
t_vec = [t1 t2 t3 t4 t5 t6].';
tdot_vec = [diff(t1(t),t) diff(t2(t),t) diff(t3(t),t) diff(t4(t),t) diff(t5(t),t) diff(t6(t),t)].';
tdotdot_vec = [diff(t1dot(t),t) diff(t2dot(t),t) diff(t3dot(t),t) diff(t4dot(t),t) diff(t5dot(t),t) diff(t6dot(t),t)].';
tdotsh_vec = [t1dot(t) t2dot(t) t3dot(t) t4dot(t) t5dot(t) t6dot(t)].';
tdotdotsh_vec = [t1dotdot(t) t2dotdot(t) t3dotdot(t) t4dotdot(t) t5dotdot(t) t6dotdot(t)].';
qdot_vec = [q1dot q2dot q3dot q4dot q5dot q6dot].';
qdotdot_vec = [q1dotdot q2dotdot q3dotdot q4dotdot q5dotdot q6dotdot].';

% Load lagrangian matrices M, Cqdot and G
load('MCG');

% Set Joint limits
q1lim = 1.5707;
q2lim = 0.7854;
q3lim = 0.24;
q4lim = 2.2689;
q5lim = 1.5707;
q6lim = 1.3963;
qUlim = [q1lim q2lim q3lim q4lim q5lim q6lim]';
qLlim = [-q1lim -q2lim 0 -q4lim -q5lim -q6lim]';

%% DH Parameters
DH = getPSM_DH();
%% Forward Kinematics
% Calculate intermediate and final transformation matrices, position
% vectors and rotation matrices
T_fin = eye(4);
T = sym('T',[4,4]);
for i = 1:size(DH,1)
    T(:,:,i) = dh2mat(DH(i,1),DH(i,2),DH(i,3),DH(i,4));
    T_fin = T_fin*T(:,:,i);
    switch i
        case 2
            T0q1 = T_fin;
            P0q1 = T0q1(1:3,4);
            R0q1 = T0q1(1:3,1:3);
        case 4
            T0q2 = T_fin;
            P0q2 = T0q2(1:3,4);
            R0q2 = T0q2(1:3,1:3);
        case 7
            T0q3 = T_fin;
            P0q3 = T0q3(1:3,4);
            R0q3 = T0q3(1:3,1:3);
        case 8
            T0com4 = T_fin;
            P0com4 = T0com4(1:3,4);
            R0com4 = T0com4(1:3,1:3);
        case 10
            T0q4 = T_fin;
            P0q4 = T0q4(1:3,4);
            R0q4 = T0q4(1:3,1:3);
        case 14
            T0com5 = T_fin;
            P0com5 = T0com5(1:3,4);
            R0com5 = T0com5(1:3,1:3);
        case 12
            T0q5 = T_fin;
            P0q5 = T0q5(1:3,4);
            R0q5 = T0q5(1:3,1:3);
        case 20
            T0com6 = T_fin;
            P0com6 = T0com6(1:3,4);
            R0com6 = T0com6(1:3,1:3);
        case 17
            T0q6 = T_fin;
            P0q6 = T0q6(1:3,4);
            R0q6 = T0q6(1:3,1:3);
        case 23
            T0tip = T_fin;
            P0tip = T0tip(1:3,4);
            R0tip = T0tip(1:3,1:3);
        case 25
            T0com7 = T_fin;
            P0com7 = T0com7(1:3,4);
            R0com7 = T0com7(1:3,1:3);
    end
            
end
P_fin = T_fin(1:3,4); %Extracting only the x-y-z position from T matrix

P_fint = subs(P_fin,q_vec,t_vec); % substitute time dependence
Pdottemp = diff(P_fint,t); % differentiate w.r.t. time
Pdot = subs(Pdottemp,[t_vec,tdot_vec],[q_vec,qdot_vec]); % substitute time independence

% Initialize the robot at some position as home [0.5 0.5 0.2 0.2 0.2 0.2]';
qHome = [0.5 0.5 0.2 0.2 0.2 0.2]';
% Calculate numerical transformation matrix
THome = subs(T_fin,q_vec,qHome);

%% Geometric and Analytical Jacobian Calculation
% Calculate the first three rows of Jacobian (common to Analytical and
% Geometric Jacobians)
J3 = jacobian(P_fin,q_vec);
% omega = last three rows of geometric Jacobian as defined by axis of
% rotation with respect to the base frame.
omega = [ 0  -1   0   0  -1   0;...
         -1   0   0   0   0  -1;...
          0   0   0  -1   0   0];
% Combine to form the Geometric Jacobian
Jgeom = [J3;omega];

% Calculate Euler Angles from Transformation Matrix using Peter Corke's
% tr2eul function
eul = tr2eul(THome);
phi = eul(1);theta = eul(2);xi = eul(3);
om2eulmat = [0 -sin(phi) cos(phi)*sin(theta);...
             0  cos(phi) sin(phi)*sin(theta);...
             1     0    cos(theta)];
% Calculate (numerically) change in Euler angles with time
euldot = om2eulmat\omega;

% Symbolic Analytical Jacobian
Janasym = [J3;euldot];
% Symbolic First derivative of Analytical Jacobian
Janadot = subs(diff(subs(Janasym,q_vec,t_vec),t),[t_vec,tdot_vec],[q_vec,qdot_vec]);

% Initialize variables for iterations
j=0;
dt = 0.01;time = 6;
iter = 0:dt:time;
errthresh = 0.01;
qnew=qHome;
q=qnew;
qdotnew = [0 0 0 0 0 0]';
euldot_old = euldot;

%% Cubic Trajectory Generation
% Generate Cubic Trajectories for desired position and orientation 
% with the cubicTraj() function 
% Credit: Prof. Jie Fu, Robotics Engineering, WPI.
polytime = 3.5;
xco = cubicTraj(THome(1,4),0,0,0,polytime,1);
yco = cubicTraj(THome(2,4),0,0,0,polytime,1);
zco = cubicTraj(THome(3,4),0,-0.2,0,polytime,1);
phico = cubicTraj(phi,0,deg2rad(0),0,polytime,1);
thco = cubicTraj(theta,0,deg2rad(180),0,polytime,1);
xico = cubicTraj(xi,0,deg2rad(0),0,polytime,1);

%% Initialize parameters for iterations
% Desired Mass of end effector
Md = 0.2*eye(6);
% Desired spring of end effector
Kp = 5*eye(6);
% Desired Damping coefficient of end effector
Kd = 3*eye(6);

% End effector initial force vector
he = [0 0 0 0 0 0]';

% Initial end effector actual cartesian position
xe = [THome(1:3,4);phi;theta;xi];
% Initial end effector desired cartesian acceleration
xddoubled = [0 0 0 0 0 0]';
% Initial end effector desired cartesian velocity
xddot = [0 0 0 0 0 0]';
% Initial end effector actual cartesian velocity
xedot = [subs(Pdot,[q_vec,qdot_vec],[qnew,qdotnew]);0;0;0];
% Initial difference in end effector desired and actual cartesian 
% velocities
xtildedot = xddot - xedot;    

%% Pre-allocating size of matrices for simulation iterations
xori_vec(3,size(iter,2)) = 0;
yori_vec(3,size(iter,2)) = 0;
zori_vec(3,size(iter,2)) = 0;
xdhist(6,size(iter,2)) = 0;
xehist(6,size(iter,2)) = 0;
qhist(6,size(iter,2)) = 0;
EEpos_fin(3,6,size(iter,2)) = 0;
j1pos_fin(3,6,size(iter,2)) = 0;
j2pos_fin(3,6,size(iter,2)) = 0;
j3pos_fin(3,6,size(iter,2)) = 0;
j4pos_fin(3,6,size(iter,2)) = 0;
j5pos_fin(3,6,size(iter,2)) = 0;
j6pos_fin(3,6,size(iter,2)) = 0;

qlimflag = zeros(size(qnew,1),1);

stab_var(6,size(iter,2)) = 0;

EEpos = double(subs(P_fin,q_vec,qHome));
j1pos = double(subs(P0q1,q_vec,qHome));
j2pos = double(subs(P0q2,q_vec,qHome));
j3pos = double(subs(P0q3,q_vec,qHome));
j4pos = double(subs(P0q4,q_vec,qHome));
j5pos = double(subs(P0q5,q_vec,qHome));
j6pos = double(subs(P0q6,q_vec,qHome));

X = [j3pos(1,:);j4pos(1,:);j5pos(1,:);j6pos(1,:); EEpos(1,:)];
Y = [j3pos(2,:);j4pos(2,:);j5pos(2,:);j6pos(2,:); EEpos(2,:)];
Z = [j3pos(3,:);j4pos(3,:);j5pos(3,:);j6pos(3,:); EEpos(3,:)];

scale = 0.2; % scale for quiver plots of line model of PSM
scale2 = 0.005; % scale for quiver plots of zoomed view of end effector
despos = sprintf('Desired end position: \n [0 0 -0.2]');
desori = sprintf('Desired orientation: \n [0 180 0]');

% Save variables that will be laoded into call back function of 
% interactive sliders
save('vars_for_cb','iter','dt','M','Cqdot','G','R0tip','Janasym','Janadot','Jgeom','xco','yco','zco','phico','thco','xico','Md','Kp','Kd','P_fin','P0q1','P0q2','P0q3','P0q4','P0q5','P0q6','q_vec','qdot_vec','xe','xddoubled','xtildedot','j','qnew','qdotnew','g','omega','euldot_old','qUlim','qLlim','Pdot','xddot','qlimflag','scale','scale2','polytime');

% Initial plot with sliders waiting for user input of force
figure;
subplot(1,2,1);
h = plot3(X(:,1),Y(:,1),Z(:,1),'b',X(:,1),Y(:,1),Z(:,1),'mo','LineWidth',5,'MarkerSize',5);
hold on;
h3 = plot3(X(5,1),Y(5,1),Z(5,1),'b','LineWidth',2);
hold on;
plot3(0,0,0,'go','markersize',5,'linewidth',5);
axis([-0.4 0.4 -0.4 0.4 -0.4 0.4]);
xlabel('X axis');
ylabel('Y axis');
zlabel('Z axis');
title('Line model of da Vinci PSM','Fontsize',18,'color','r','fontweight','bold');
txth = text(X(5,1),Y(5,1),Z(5,1),['(' num2str(X(5,1),'%.2f') ',' num2str(Y(5,1),'%.2f') ',' num2str(Z(5,1),'%.2f') ')'],'fontsize',14);
txtf1 = text(-1.5,-0.5,1.55,'Force Applied in X: 0 Newton','fontsize',10);
txtf2 = text(-1.5,-0.5,1.5,'Force Applied in Y: 0 Newton','fontsize',10);
txtf3 = text(-1.5,-0.5,1.45,'Force Applied in Z: 0 Newton','fontsize',10);
txtf7 = text(-1.5,-0.5,1.4,'Moment Applied about X: 0 Newton-Meter','fontsize',10);
txtf8 = text(-1.5,-0.5,1.35,'Moment Applied about Y: 0 Newton-Meter','fontsize',10);
txtf9 = text(-1.5,-0.5,1.3,'Moment Applied about Z: 0 Newton-Meter','fontsize',10);
txtf4 = text(-1.5,-0.5,1.23,despos,'fontsize',10);
txtf10 = text(-1.5,-0.5,1.14,desori,'fontsize',10);
txtf5 = text(-1.5,-0.5,1.05,'Initial position','fontsize',10,'color','b');
txtf6 = text(-1.5,-0.5,1,'Iteration: 1','fontsize',10,'color','b');

% Draw vectors representing orientation of end effector
qx = quiver3(THome(1,4),THome(2,4),THome(3,4),THome(1,1)*scale,THome(2,1)*scale,THome(3,1)*scale,...
    'r','DisplayName','X');
qy = quiver3(THome(1,4),THome(2,4),THome(3,4),THome(1,2)*scale,THome(2,2)*scale,THome(3,2)*scale,...
    'g','DisplayName','Y');
qz = quiver3(THome(1,4),THome(2,4),THome(3,4),THome(1,3)*scale,THome(2,3)*scale,THome(3,3)*scale,...
    'b','DisplayName','Z');

% Plot Zoomed view of end effector
subplot(1,2,2);
h2 = plot3(X(3:5,1),Y(3:5,1),Z(3:5,1),'b',X(3:5,1),Y(3:5,1),Z(3:5,1),'mo','LineWidth',5,'MarkerSize',5);
title('Zoomed View of End effector','Fontsize',18,'color','r','fontweight','bold');
xlabel('X axis');
ylabel('Y axis');
zlabel('Z axis');
hold on
qx2 = quiver3(THome(1,4),THome(2,4),THome(3,4),THome(1,1)*scale2,THome(2,1)*scale2,THome(3,1)*scale2,...
    'r','DisplayName','X');
qy2 = quiver3(THome(1,4),THome(2,4),THome(3,4),THome(1,2)*scale2,THome(2,2)*scale2,THome(3,2)*scale2,...
    'g','DisplayName','Y');
qz2 = quiver3(THome(1,4),THome(2,4),THome(3,4),THome(1,3)*scale2,THome(2,3)*scale2,THome(3,3)*scale2,...
    'b','DisplayName','Z');

%% User Interface
% Sliders definition
hs1 = uicontrol('style','slider','position',[25,5,120,20],'min',-2,'max',2);
hs2 = uicontrol('style','slider','position',[180,5,120,20],'min',-2,'max',2);
hs3 = uicontrol('style','slider','position',[335,5,120,20],'min',-2,'max',2);
hs4 = uicontrol('style','slider','position',[490,5,120,20],'min',-2,'max',2);
hs5 = uicontrol('style','slider','position',[645,5,120,20],'min',-2,'max',2);
hs6 = uicontrol('style','slider','position',[800,5,120,20],'min',-2,'max',2);

% Text for sliders
ht1t = uicontrol('style','text','position',[5,5,20,20],'string','Fx');
ht2t = uicontrol('style','text','position',[160,5,20,20],'string','Fy');
ht3t = uicontrol('style','text','position',[315,5,20,20],'string','Fz');
ht4t = uicontrol('style','text','position',[470,5,20,20],'string','Ux');
ht5t = uicontrol('style','text','position',[625,5,20,20],'string','Uy');
ht6t = uicontrol('style','text','position',[780,5,20,20],'string','Uz');

% Text displaying slider values
ht1v = uicontrol('style','text','position',[65,25,40,20],'string','0');
ht2v = uicontrol('style','text','position',[220,25,40,20],'string','0');
ht3v = uicontrol('style','text','position',[375,25,40,20],'string','0');
ht4v = uicontrol('style','text','position',[530,25,40,20],'string','0');
ht5v = uicontrol('style','text','position',[685,25,40,20],'string','0');
ht6v = uicontrol('style','text','position',[840,25,40,20],'string','0');

% Push buttons for "confirm" and "apply force"
hpb_conf = uicontrol('style','pushbutton','string','Confirm','position',[935,5,100,20]);
hpb_apply = uicontrol('style','pushbutton','string','Apply Force','position',[1050,5,100,20]);

% Set the appropriate callback functions for the corresponding UI elements
set(hs1,'callback',{@cb_txtupdate,ht1v});
set(hs2,'callback',{@cb_txtupdate,ht2v});
set(hs3,'callback',{@cb_txtupdate,ht3v});
set(hs4,'callback',{@cb_txtupdate,ht4v});
set(hs5,'callback',{@cb_txtupdate,ht5v});
set(hs6,'callback',{@cb_txtupdate,ht6v});
set(hpb_conf,'callback',{@cb_conf,hs1,hs2,hs3,hs4,hs5,hs6});
set(hpb_apply,'callback',{@cb_interactive_IMP,hs1,hs2,hs3,hs4,hs5,hs6,h,h2,ht6v,6,txtf1,txtf2,txtf3,txtf4,txtf5,txtf6,txth,qx,qy,qz,qx2,qy2,qz2,txtf7,txtf8,txtf9,h3});