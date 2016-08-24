%% Function to plot the daVinci PSM in task space
function [T,h,qx,qy,qz] = plotarm_rad_6q(qq,idx)
%% Initialization
q1 = qq(1);
q2 = qq(2);
q3 = qq(3);
q4 = qq(4);
q5 = qq(5);
q6 = qq(6);
lRcc = 0.4318;
lTool = 0.4162;
lP2Y = 0.0091;
lY2CP = 0.0102;

%% DH Parameters
DH = [q1+pi/2,     0,         pi/2,      0;...
      q2-pi/2,     0,         -pi/2,     0;...
      0,           q3-lRcc,   pi/2,      0;...
      q4,          lTool,     0,         0;...
      q5-pi/2,     0,         -pi/2,     0;...
      q6-pi/2,     0,         -pi/2,     lP2Y;...
      0,           lY2CP,     -pi/2,     0];
  
%% Transformation matrix (Forward Kinematics)
T = eye(4);
for i = 1:size(DH,1)
    T_temp(:,:,i) = dh2mat(DH(i,1),DH(i,2),DH(i,3),DH(i,4));
    T = T*T_temp(:,:,i);
    Ti(:,i) = T(1:3,4);
end
T
%% Separate X,Y,Z positions
X = [Ti(1,:)];
Y = [Ti(2,:)];
Z = [Ti(3,:)];

%% Plotting
% If first plot (idx=1), plot in black else plot in blue.
if(idx==1)
    h = plot3(X,Y,Z,'k',X,Y,Z,'mo','LineWidth',5,'MarkerSize',5);
    title('Black - Home config; Blue - Desired config');
else
    h = plot3(X,Y,Z,'b',X,Y,Z,'mo','LineWidth',5,'MarkerSize',5);
end
drawnow;
hold on;
% Plot orientation vectors.
qx = quiver3(T(1,4),T(2,4),T(3,4),T(1,1)*0.05,T(2,1)*0.05,T(3,1)*0.05,...
    'r','DisplayName','X');
qy = quiver3(T(1,4),T(2,4),T(3,4),T(1,2)*0.05,T(2,2)*0.05,T(3,2)*0.05,...
    'g','DisplayName','Y');
qz = quiver3(T(1,4),T(2,4),T(3,4),T(1,3)*0.05,T(2,3)*0.05,T(3,3)*0.05,...
    'b','DisplayName','Z');
axis equal;
hold on;
legend([qx qy qz],'X','Y','Z');
end