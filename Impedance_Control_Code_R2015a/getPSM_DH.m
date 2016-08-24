function DH = getPSM_DH()
syms q1 q2 q3 q4 q5 q6 real;
lRcc = 0.4318;      % meters
lTool = 0.4162;     % meters
lP2Y = 0.0091;      % meters
lY2CP = 0.0102;     % meters
com4 = -31.62/1000; %meters % center of mass of link 4
com5 = 3.28/1000;   %meters % center of mass of link 5
com6 = 5.80/1000;  %meters % center of mass of link 6
com7 = -0.52/1000;  %meters % center of mass of link 7
DH = [0 0 pi/2 0;...
      pi/2 0 0 0;... % joint 1
      q1 0 -pi/2 0;...
      -pi/2 0 0 0;... % joint 2 
      q2 0 pi/2 0;...
      0 -lRcc 0 0;...
      0 q3 0 0;... % joint 3
      0 com4 0 0;... % go to center of mass of link 4
      0 -com4 0 0;... % return to joint 3
      0 lTool 0 0;... % joint 4
      q4 0 -pi/2 0;... 
      -pi/2 0 0 0;... % joint 5
      q5 0 0 0;...
      0 0 0 com5;... % go to center of mass of link 5
      0 0 0 -com5;... % return to joint 4
      0 0 -pi/2 lP2Y;...
      -pi/2 0 0 0;... % joint 6
      q6 0 0 0;...
      pi/2 0 0 0;...
      0 0 0 com6;... % go to center of mass of link 6
      0 0 0 -com6;... % return to joint 6
      -pi/2 0 0 0;...
      0 0 -pi/2 0;...
      0 lY2CP 0 0;... % tool tip
      0 com7 0 0;... %go to center of mass of link 7
      0 -com7 0 0]; % return to tool tip
end