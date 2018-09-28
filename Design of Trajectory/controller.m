function [F, M] = controller(t, state, des_state, params)
%CONTROLLER  Controller for the quadrotor
%
%   state: The current state of the robot with the following fields:
%   state.pos = [x; y; z], state.vel = [x_dot; y_dot; z_dot],
%   state.rot = [phi; theta; psi], state.omega = [p; q; r]
%
%   des_state: The desired states are:
%   des_state.pos = [x; y; z], des_state.vel = [x_dot; y_dot; z_dot],
%   des_state.acc = [x_ddot; y_ddot; z_ddot], des_state.yaw,
%   des_state.yawdot
%
%   params: robot parameters

%   Using these current and desired states, you have to compute the desired
%   controls


% =================== Your code goes here ===================

% Variable initialization
m = params.mass;
g = params.gravity;
I = params.I;
inv_I = params.invI;
L = params.arm_length;

x = state.pos(1);
y = state.pos(2);
z = state.pos(3);
x_dot = state.vel(1);
y_dot = state.vel(2);
z_dot = state.vel(3);
phi = state.rot(1);
theta = state.rot(2);
psi = state.rot(3);
p = state.omega(1);
q = state.omega(2);
r = state.omega(3);

des_x = des_state.pos(1);
des_y = des_state.pos(2);
des_z = des_state.pos(3);
des_x_dot = des_state.vel(1);
des_y_dot = des_state.vel(2);
des_z_dot = des_state.vel(3);
des_x_ddot = des_state.acc(1);
des_y_ddot = des_state.acc(2);
des_z_ddot = des_state.acc(3);
des_psi = des_state.yaw;
des_psi_dot = des_state.yawdot;

% Constants
k_d_x = 30;
k_p_x = 3;
k_d_y = 30;
k_p_y = 3;
k_p_z = 800;
k_d_z = 30;
k_p_phi = 160;
k_d_phi = 3;
k_p_theta = 160;
k_d_theta = 3;
k_p_psi = 160;
k_d_psi = 3;

% Commanded accelerations
commanded_r_ddot_x = des_x_ddot + k_d_x * (des_x_dot - x_dot) + k_p_x * (des_x - x);
commanded_r_ddot_y = des_y_ddot + k_d_y * (des_y_dot - y_dot) + k_p_y * (des_y - y);
commanded_r_ddot_z = des_z_ddot + k_d_z * (des_z_dot - z_dot) + k_p_z * (des_z - z);

% Thrust
F = m * (g + commanded_r_ddot_z);

% Moment
p_des = 0;
q_des = 0;
r_des = des_psi_dot;
des_phi = 1/g * (commanded_r_ddot_x * sin(des_psi) - commanded_r_ddot_y * cos(des_psi));
des_theta = 1/g * (commanded_r_ddot_x * cos(des_psi) + commanded_r_ddot_y * sin(des_psi));
M = [k_p_phi * (des_phi - phi) + k_d_phi * (p_des - p);
     k_p_theta * (des_theta - theta) + k_d_theta * (q_des - q);
     k_p_psi * (des_psi - psi) + k_d_psi * (r_des - r)];

% =================== Your code ends here ===================

end
