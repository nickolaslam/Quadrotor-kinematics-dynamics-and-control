function [ desired_state ] = traj_generator(t, state, waypoints)
% TRAJ_GENERATOR: Generate the trajectory passing through all
% positions listed in the waypoints list
%
% NOTE: This function would be called with variable number of input arguments.
% During initialization, it will be called with arguments
% trajectory_generator([], [], waypoints) and later, while testing, it will be
% called with only t and state as arguments, so your code should be able to
% handle that. This can be done by checking the number of arguments to the
% function using the "nargin" variable, check the MATLAB documentation for more
% information.
%
% t,state: time and current state (same variable as "state" in controller)
% that you may use for computing desired_state
%
% waypoints: The 3xP matrix listing all the points you much visited in order
% along the generated trajectory
%
% desired_state: Contains all the information that is passed to the
% controller for generating inputs for the quadrotor
%
% It is suggested to use "persistent" variables to store the waypoints during
% the initialization call of trajectory_generator.


persistent waypoints0 traj_time d0
persistent x_coefficients y_coefficients z_coefficients
if nargin > 2
    d = waypoints(:,2:end) - waypoints(:,1:end-1);
    d0 = 2 * sqrt(d(1,:).^2 + d(2,:).^2 + d(3,:).^2);
    traj_time = [0, cumsum(d0)];
    waypoints0 = waypoints;
    
    x_coefficients = getTrajCoefficients(waypoints0(1,1:end))';
    y_coefficients = getTrajCoefficients(waypoints0(2,1:end))';
    z_coefficients = getTrajCoefficients(waypoints0(3,1:end))';
    
else
    if(t > traj_time(end))
        t = traj_time(end) - 0.0001;
    end
    
    t_index = find(traj_time >= t,1)-1; %between 1:n
    
    if (t_index == 0)
        t_index = 1;
    end
    if(t == 0)
        desired_state.pos = waypoints0(:,1);
        desired_state.vel = 0*waypoints0(:,1);
        desired_state.acc = 0*waypoints0(:,1);
    else
        scale = (t-traj_time(t_index))/d0(t_index);
        index = (t_index - 1) * 8 + 1 : t_index * 8;
        
        p = getPolyCoefficients(8, 0, scale)';
        desired_state.pos = [x_coefficients(index) * p; y_coefficients(index) * p; z_coefficients(index) * p];
        
        pdot = getPolyCoefficients(8, 1, scale)';
        desired_state.vel = [x_coefficients(index) * pdot; y_coefficients(index) * pdot; z_coefficients(index) * pdot].*(1/d0(t_index));
        
        pddot = getPolyCoefficients(8, 2, scale)';
        desired_state.acc = [x_coefficients(index) * pddot; y_coefficients(index) * pddot; z_coefficients(index) * pddot].*(1/d0(t_index)^2);
        
    end
    desired_state.yaw = 0;
    desired_state.yawdot = 0;
end
end

function [coefficients] = getTrajCoefficients(waypoints)
    n = size(waypoints,2) - 1;
    A = zeros(8*n, 8*n);
    b = zeros(8*n, 1);
    
    % constraint 1: p_i(Si-1) = wi-1
    % pattern: 1 0 0 0 0 0 0 0
    for i = 1:n
        x_shift = 8*(i-1);
        A(i,x_shift+1) = 1;
        b(i) = waypoints(i);
    end
    
    % constraint 2: p_i(Si) = wi
    % pattern: 1 1 1 1 1 1 1 1
    y_shift = n;
    for i = 1:n
        x_shift = 8*(i-1);
        A(y_shift+i,x_shift+1:x_shift+8) = 1;
        b(y_shift+i) = waypoints(i+1);
    end
    
    % constraint 3: p_1^k(S0) = 0
    % values:
    % 0 1 0 0 0 0 0 0
    % 0 0 2 0 0 0 0 0
    % 0 0 0 6 0 0 0 0
    y_shift = 2*n;
    for k = 1:3
        A(y_shift+k,1:8) = getPolyCoefficients(8, k, 0);
    end
    
    % constraint 4: p_n^k(Sn) = 0
    % values:
    % 0 1 0 0 0 0 0 0
    % 0 0 2 0 0 0 0 0
    % 0 0 0 6 0 0 0 0
    y_shift = 2*n+3;
    for k = 1:3
        A(y_shift+k, end-7:end) = getPolyCoefficients(8, k, 1);
    end
    
    % constraint 5: p_i^k(Si) = p^k_i+1(Si)
    for i = 1:n-1
        x_shift = 8*(i-1);
        for k = 1:6
            y_shift = 2*n+6 + 6*(i-1);
            A(y_shift+k, x_shift+1:x_shift+16) = [getPolyCoefficients(8, k, 1), -1 * getPolyCoefficients(8, k, 0)];
        end
    end
    
    coefficients = inv(A) * b;
    
end

function [elements] = getPolyCoefficients(length, derivative, value)
    coefficients = ones(1, length);
    exponents = zeros(1, length);
    for i = 1:length
        exponents(1,i) = i-1;
    end
    for d = 1:derivative
        for i = 1:length
            coefficients(1,i) = coefficients(1,i) * exponents(1,i);
            exponents(1,i) = max(0, exponents(1,i) - 1);
        end
    end
    elements = coefficients .* (value .^ exponents);
end