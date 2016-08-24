% This script was written in MATLAB R2012a
% The output of this file is the MCG.mat file which contains the symbolic
% variables M, Cqdot and G which are representative of the dynamics of the
% system. This is later loaded into the MAIN.m file in this package
%% Initialization of variables
clear all;close all;clc;
% Load the symbolic torque values calculated and saved by the
% proper_lagrangian_PSM.m file.
load('taus');
syms q_vec q1 q2 q3 q4 q5 q6 g real;
syms t t1(t) t2(t) t3(t) t4(t) t5(t) t6(t) real;
syms q1dot q2dot q3dot q4dot q5dot q6dot real;
syms t1dot(t) t2dot(t) t3dot(t) t4dot(t) t5dot(t) t6dot(t) real;
syms q1dotdot q2dotdot q3dotdot q4dotdot q5dotdot q6dotdot real;
syms t1dotdot(t) t2dotdot(t) t3dotdot(t) t4dotdot(t) t5dotdot(t) t6dotdot(t) real;
q_vec = [q1 q2 q3 q4 q5 q6].';
t_vec = [t1 t2 t3 t4 t5 t6].';
tdot_vec = [diff(t1(t),t) diff(t2(t),t) diff(t3(t),t) diff(t4(t),t) diff(t5(t),t) diff(t6(t),t)].';
tdotdot_vec = [diff(t1dot(t),t) diff(t2dot(t),t) diff(t3dot(t),t) diff(t4dot(t),t) diff(t5dot(t),t) diff(t6dot(t),t)].';
tdotsh_vec = [t1dot(t) t2dot(t) t3dot(t) t4dot(t) t5dot(t) t6dot(t)].';
tdotdotsh_vec = [t1dotdot(t) t2dotdot(t) t3dotdot(t) t4dotdot(t) t5dotdot(t) t6dotdot(t)].';
qdot_vec = [q1dot q2dot q3dot q4dot q5dot q6dot].';
qdotdot_vec = [q1dotdot q2dotdot q3dotdot q4dotdot q5dotdot q6dotdot].';

% Take the coefficients of q_double_dot and g
[t1_val,c1] = coeffs(tau1,[qdotdot_vec.',g]);
[t2_val,c2] = coeffs(tau2,[qdotdot_vec.',g]);
[t3_val,c3] = coeffs(tau3,[qdotdot_vec.',g]);
[t4_val,c4] = coeffs(tau4,[qdotdot_vec.',g]);
[t5_val,c5] = coeffs(tau5,[qdotdot_vec.',g]);
[t6_val,c6] = coeffs(tau6,[qdotdot_vec.',g]);

% Separate out M(q) matrix according to the coefficients
M = [t1_val(1) t1_val(2) t1_val(3) t1_val(4) t1_val(5) t1_val(6);...
     t2_val(1) t2_val(2) t2_val(3) t2_val(4) t2_val(5) t2_val(6);...
     t3_val(1) t3_val(2) t3_val(3)    0      t3_val(4) t3_val(5);...
     t4_val(1) t4_val(2)    0      t4_val(3) t3_val(4) t3_val(5);...
     t5_val(1) t5_val(2) t5_val(3) t5_val(4) t5_val(5) t5_val(6);...
     t6_val(1) t6_val(2) t6_val(3) t6_val(4) t6_val(5) t6_val(6)];
% Separate out G(q) vector according to the coefficients
G = [t1_val(7);t2_val(7);t3_val(6);t4_val(6);t5_val(7);t6_val(7)].*g;
% Everything else is C(q,qdot)*qdot
Cqdot = [0;t2_val(8);t3_val(7);t4_val(7);t5_val(8);t6_val(8)];

% Save M, Cqdot and G
save('MCG','M','Cqdot','G');