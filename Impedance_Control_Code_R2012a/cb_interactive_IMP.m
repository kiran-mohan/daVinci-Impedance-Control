function cb_interactive_IMP(hs,~,hs1,hs2,hs3,hs4,hs5,hs6,h,h2,htv,num,txtf1,txtf2,txtf3,txtf4,txtf5,txtf6,txth,qx,qy,qz,qx2,qy2,qz2,txtf7,txtf8,txtf9,h3)
% This file is the callback function for the "Apply Force" button.
% This function has the impedance control law and the plot updates.
% Load variables stored in MAIN.m that are required for the iterations
load('vars_for_cb');
j=0;

% Get values of forces set by the user using sliders
Fx = get(hs1,'value');
Fy = get(hs2,'value');
Fz = get(hs3,'value');
Ux = get(hs4,'value');
Uy = get(hs5,'value');
Uz = get(hs6,'value');
he_app = [Fx Fy Fz Ux Uy Uz]';

% Start iterating
for i = iter
tic;

    j=j+1;
    disp(['Iteration: ',num2str(j)]);
    if(j>=130 && j<=150)
        he = he_app;
    else
        he = [0 0 0 0 0 0]';
    end
    hehist(:,j) = he;
    M_num = double(subs(M,q_vec,qnew));
    % Cqdot and G are not used in the control law as they get cancelled
    % out. However, they are defined here for consistency but are
    % commented out.
%     Cqdot_num = subs(Cqdot,[q_vec,qdot_vec],[qnew,qdotnew]);
%     G_num = subs(G,[g,q_vec.'],[9.81,qnew.']);
    
    % Calculate numeric rotation matrix
    R0tip_num = double(subs(R0tip,q_vec,qnew));
    % Save orientation vectors over time
    xori_vec(:,j) = R0tip_num(:,1);
    yori_vec(:,j) = R0tip_num(:,2);
    zori_vec(:,j) = R0tip_num(:,3);
    
    % This if condition is to keep the jacobians constant for every 20
    % iterations. We could do this, but we will get plots with a step-like
    % response indicating some level of discontinuity. The iterations also
    % run faster with the if condition. Time per iteration without the if
    % condition (calculating jacobians every iteration) is 0.6 seconds while
    % time per iteration with the if condition (calculating jacobians once
    % every 20 iterations) is 0.35 seconds.
    % To get faster plots, the if condition is used by default.
    % Please comment out the if along with its corresponding end 
    % to get smooth plots.
    if(mod(j,20) == 0 || j==1)
        % Calculate Analytical Jacobian and its 1st derivative
        % numerical updates
        eul = tr2eul(R0tip_num);
        
        phi = eul(1);theta = eul(2);xi = eul(3);
        om2eulmat = [0 -sin(phi) cos(phi)*sin(theta);...
            0  cos(phi) sin(phi)*sin(theta);...
            1     0    cos(theta)];
        euldot = om2eulmat\omega;
        Janasym(4:6,:) = euldot;
        Jana_num = double(subs(Janasym,q_vec,qnew));
        Janadot_num = double(subs(Janadot,[q_vec,qdot_vec],[qnew,qdotnew]));
        Janadot_num(4:6,:) = (euldot-euldot_old);
        
        % Calcualte Geometric Jacobian numeric update
        Jgeom_num = double(subs(Jgeom,q_vec,qnew));
    end    
    
    % Calculate desired position from cubic trajectory generated previously
    % (trajectories generated in the MAIN.m file)
    if(j<=polytime*100)
        xdx = xco(1)+xco(2)*iter(j)+xco(3)*iter(j)^2+xco(4)*iter(j)^3;
        xdy = yco(1)+yco(2)*iter(j)+yco(3)*iter(j)^2+yco(4)*iter(j)^3;
        xdz = zco(1)+zco(2)*iter(j)+zco(3)*iter(j)^2+zco(4)*iter(j)^3;
        xdphi = phico(1)+phico(2)*iter(j)+phico(3)*iter(j)^2+phico(4)*iter(j)^3;
        xdth = thco(1)+thco(2)*iter(j)+thco(3)*iter(j)^2+thco(4)*iter(j)^3;
        xdxi = xico(1)+xico(2)*iter(j)+xico(3)*iter(j)^2+xico(4)*iter(j)^3;
    end
    xd = [xdx xdy xdz xdphi xdth xdxi]';
    
    % Save the desired trajectory
    xdhist(:,j) = xd;
    % Calculate xtilde
    xtilde = xd - xe;
    
    % Control law for the inverse dynamics impedance control
    % Resolved acceleration a_q calculation
    aq = Jana_num\(xddoubled + Md\Kd*xtildedot + Md\Kp*xtilde - Janadot_num*qdotnew);
    aq = double(aq);
    
    % Calculate q_double_dot (Joint accelerations)
    qdotdot = aq - M_num\Jgeom_num'*he;
    
    % Calculate joint velocities by numerically forward integrating 
    % with time.
    qdotnew = qdotnew + qdotdot.*dt;
    % Calculate joint angle values by numerically forward integrating 
    % with time.
    qnew = qnew + qdotnew.*dt;
    
    % Save a history of joint angle values over time
    qhist(:,j) = qnew;
    
    % Check for joint limits and break if joint limits exceeded
    if(sum(qnew>qUlim)>0 || sum(qnew<qLlim)>0)
        qlimflag = 1;
        break;
    end
    
    % Calculate individual joint cartesian positions
    EEpos = double(subs(P_fin,q_vec,qnew));
    j1pos = double(subs(P0q1,q_vec,qnew));
    j2pos = double(subs(P0q2,q_vec,qnew));
    j3pos = double(subs(P0q3,q_vec,qnew));
    j4pos = double(subs(P0q4,q_vec,qnew));
    j5pos = double(subs(P0q5,q_vec,qnew));
    j6pos = double(subs(P0q6,q_vec,qnew));
    
    % Save individual joint cartesian positions over iterations for
    % animating the plot
    EEpos_fin(:,j) = EEpos';
    j1pos_fin(:,j) = j1pos';
    j2pos_fin(:,j) = j2pos';
    j3pos_fin(:,j) = j3pos';
    j4pos_fin(:,j) = j4pos';
    j5pos_fin(:,j) = j5pos';
    j6pos_fin(:,j) = j6pos';
    
    % Calculate current pose (position and orientation) of the end effector
    xe = [EEpos;phi;theta;xi];
    % Save over iterations
    xehist(:,j) = xe;
    % Calculate end effector velocity
    xedot = [subs(Pdot,[q_vec,qdot_vec],[qnew,qdotnew]);euldot*qdotnew];
    % Calculate xtilde
    xtilde = xd-xe;
    % Calculate xtilde_dot
    xtildedot = xddot - xedot;
    
    % Set euldot_old as euldot
    euldot_old = euldot;
    
    disp('Still Processing... Please wait till iteration 601');
    toc
end

% Check if joint limits were exceeded during iterations
if(qlimflag == 1)
        disp('Sorry! Limits exceeded... Try again with lower forces/moments');
else % if joint limits were not exceeded, then plot.!
    % Create a vector of cartesian joint positions for each joint. 
    % Since joints 1 and 2 act about the remote center (0,0,0) and are 
    % always fixed, they are not included in these X,Y,Z vectors. 
    % Joints 3 to 6 and the end effector are included.
    X = [j3pos_fin(1,:);j4pos_fin(1,:);j5pos_fin(1,:);j6pos_fin(1,:); EEpos_fin(1,:)];
    Y = [j3pos_fin(2,:);j4pos_fin(2,:);j5pos_fin(2,:);j6pos_fin(2,:); EEpos_fin(2,:)];
    Z = [j3pos_fin(3,:);j4pos_fin(3,:);j5pos_fin(3,:);j6pos_fin(3,:); EEpos_fin(3,:)];
    
    disp('Finished Processing. Now Plotting the animation...');
    for i = 2:length(j1pos_fin)
        hold on;
        pause(0.001);
        set(h,'xdata',X(:,i),'ydata',Y(:,i),'zdata',Z(:,i));
        set(h2,'xdata',X(3:5,i),'ydata',Y(3:5,i),'zdata',Z(3:5,i));
        set(txth,'position',[X(5,i),Y(5,i),Z(5,i)],'string',['(' num2str(X(5,i),'%.2f') ',' num2str(Y(5,i),'%.2f'),',' num2str(Z(5,i),'%.2f') ')']);
        set(txtf6,'string',['Iteration: ',num2str(i)]);
        set(qx,'xdata',X(5,i),'ydata',Y(5,i),'zdata',Z(5,i),'udata',xori_vec(1,i)*scale,'vdata',xori_vec(2,i)*scale,'wdata',xori_vec(3,i)*scale);
        set(qy,'xdata',X(5,i),'ydata',Y(5,i),'zdata',Z(5,i),'udata',yori_vec(1,i)*scale,'vdata',yori_vec(2,i)*scale,'wdata',yori_vec(3,i)*scale);
        set(qz,'xdata',X(5,i),'ydata',Y(5,i),'zdata',Z(5,i),'udata',zori_vec(1,i)*scale,'vdata',zori_vec(2,i)*scale,'wdata',zori_vec(3,i)*scale);
        set(qx2,'xdata',X(5,i),'ydata',Y(5,i),'zdata',Z(5,i),'udata',xori_vec(1,i)*scale2,'vdata',xori_vec(2,i)*scale2,'wdata',xori_vec(3,i)*scale2);
        set(qy2,'xdata',X(5,i),'ydata',Y(5,i),'zdata',Z(5,i),'udata',yori_vec(1,i)*scale2,'vdata',yori_vec(2,i)*scale2,'wdata',yori_vec(3,i)*scale2);
        set(qz2,'xdata',X(5,i),'ydata',Y(5,i),'zdata',Z(5,i),'udata',zori_vec(1,i)*scale2,'vdata',zori_vec(2,i)*scale2,'wdata',zori_vec(3,i)*scale2);
        if(i>130 && i<150)
            set(h3,'xdata',X(5,i),'ydata',Y(5,i),'zdata',Z(5,i),'erasemode','none');
        else
            set(h3,'xdata',X(5,i),'ydata',Y(5,i),'zdata',Z(5,i),'erasemode','none');
        end
        if(i==130)
            set(txtf1,'string',['Force Applied in X: ',num2str(Fx), ' Newton'],'color','b','fontweight','bold');
            set(txtf2,'string',['Force Applied in Y: ',num2str(Fy), ' Newton'],'color','b','fontweight','bold');
            set(txtf3,'string',['Force Applied in Z: ',num2str(Fz), ' Newton'],'color','b','fontweight','bold');
            set(txtf5,'string','Force being applied on end effector');
            set(txtf7,'string',['Moment Applied about X: ',num2str(Ux), ' Newton-Meter'],'color','b','fontweight','bold');
            set(txtf8,'string',['Moment Applied about Y: ',num2str(Uy), ' Newton-Meter'],'color','b','fontweight','bold');
            set(txtf9,'string',['Moment Applied about Z: ',num2str(Uz), ' Newton-Meter'],'color','b','fontweight','bold');
        elseif(i==150)
            set(txtf1,'string','Force Applied in X: 0 Newton','color','k');
            set(txtf2,'string','Force Applied in Y: 0 Newton','color','k');
            set(txtf3,'string','Force Applied in Z: 0 Newton','color','k');
            set(txtf7,'string','Moment Applied about X: 0 Newton-Meter','color','k');
            set(txtf8,'string','Moment Applied about Y: 0 Newton-Meter','color','k');
            set(txtf9,'string','Moment Applied about Z: 0 Newton-Meter','color','k');
            set(txtf5,'string','Force released');
        end
    end
    
    % Plot X,Y,Z positions of end effector over time.
    figure;
    subplot(4,1,1);
    plot(iter,hehist(1,:),iter,hehist(2,:),iter,hehist(3,:),iter,hehist(4,:),iter,hehist(5,:),iter,hehist(6,:),'linewidth',2);
    title('Force on end effector');
    legend('Fx','Fy','Fz','Ux','Uy','Uz');
    xlabel('time');
    ylabel('Force (N)');
    subplot(4,1,2);
    plot(iter,xdhist(1,:),'b',iter,xehist(1,:),'g','linewidth',2);
    title('End effector X component desired trajectory vs actual trajectory');
    legend('desired','actual');
    xlabel('time');
    ylabel('EE X motion (m)');
    subplot(4,1,3);
    plot(iter,xdhist(2,:),'b',iter,xehist(2,:),'g','linewidth',2);
    title('End effector Y component desired trajectory vs actual trajectory');
    legend('desired','actual');
    xlabel('time');
    ylabel('EE Y motion (m)');
    subplot(4,1,4);
    plot(iter,xdhist(3,:),'b',iter,xehist(3,:),'g','linewidth',2);
    title('End effector Z component desired trajectory vs actual trajectory');
    legend('desired','actual');
    xlabel('time');
    ylabel('EE Z motion (m)');
    
    % Plot phi,theta,xi orientation of end effector over time.
    figure;
    subplot(4,1,1);
    plot(iter,hehist(1,:),iter,hehist(2,:),iter,hehist(3,:),iter,hehist(4,:),iter,hehist(5,:),iter,hehist(6,:),'linewidth',2);
    title('Force on end effector');
    legend('Fx','Fy','Fz','Ux','Uy','Uz');
    xlabel('time');
    ylabel('Force (N)');
    subplot(4,1,2);
    plot(iter,xdhist(4,:),'b',iter,xehist(4,:),'g','linewidth',2);
    title('End effector phi component desired trajectory vs actual trajectory');
    legend('desired','actual');
    xlabel('time');
    ylabel('EE phi (rad)');
    subplot(4,1,3);
    plot(iter,xdhist(5,:),'b',iter,xehist(5,:),'g','linewidth',2);
    title('End effector theta component desired trajectory vs actual trajectory');
    legend('desired','actual');
    xlabel('time');
    ylabel('EE theta (rad)');
    subplot(4,1,4);
    plot(iter,xdhist(6,:),'b',iter,xehist(6,:),'g','linewidth',2);
    title('End effector xi component desired trajectory vs actual trajectory');
    legend('desired','actual');
    xlabel('time');
    ylabel('EE xi (rad)');
    
    % Plot joints 1,2,3 joint space positions over time.
    figure;
    subplot(4,1,1);
    plot(iter,hehist(1,:),iter,hehist(2,:),iter,hehist(3,:),iter,hehist(4,:),iter,hehist(5,:),iter,hehist(6,:),'linewidth',2);
    title('Force on end effector');
    legend('Fx','Fy','Fz','Ux','Uy','Uz');
    xlabel('time');
    ylabel('Force (N)');
    subplot(4,1,2);
    plot(iter,qhist(1,:),'r');
    title('Joints 1 joint-space trajectory');
    legend('actual');
    xlabel('time');
    ylabel('Joint motion (rad)');
    subplot(4,1,3);
    plot(iter,qhist(2,:),'r');
    title('Joints 2 joint-space trajectory');
    legend('actual');
    xlabel('time');
    ylabel('Joint motion (rad)');
    subplot(4,1,4);
    plot(iter,qhist(3,:),'r');
    title('Joints 3 joint-space trajectory');
    legend('actual');
    xlabel('time');
    ylabel('Joint motion (rad)');
    
    % Plot joints 4,5,6 joint space positions over time.
    figure;
    subplot(4,1,1);
    plot(iter,hehist(1,:),iter,hehist(2,:),iter,hehist(3,:),iter,hehist(4,:),iter,hehist(5,:),iter,hehist(6,:),'linewidth',2);
    title('Force on end effector');
    legend('Fx','Fy','Fz','Ux','Uy','Uz');
    xlabel('time');
    ylabel('Force (N)');
    subplot(4,1,2);
    plot(iter,qhist(4,:),'r');
    title('Joints 4 joint-space trajectory');
    legend('actual');
    xlabel('time');
    ylabel('Joint motion (rad)');
    subplot(4,1,3);
    plot(iter,qhist(5,:),'r');
    title('Joints 5 joint-space trajectory');
    legend('actual');
    xlabel('time');
    ylabel('Joint motion (rad)');
    subplot(4,1,4);
    plot(iter,qhist(6,:),'r');
    title('Joints 6 joint-space trajectory');
    legend('actual');
    xlabel('time');
    ylabel('Joint motion (rad)');
    
end
end