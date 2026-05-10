function ctrl = controller_config()
% CONTROLLER_CONFIG  Returns gains and parameters for inner/outer loop controllers.
%
% Inner loop:  Geometric SO(3)  AND  INDI  (selectable at runtime via ctrl.inner.type)
% Outer loop:  Nonlinear MPC with fmincon
% Allocator:   Pseudoinverse over [T_1..T_4, delta_e, delta_a, delta_r]^T

%% ----- Inner attitude loop -----
ctrl.inner.type = 'SO3';   % 'SO3' | 'INDI'  (toggleable in main.m)

% Geometric SO(3) gains scaled by tailsitter inertia
% Tuning rule (Lee-Leok-McClamroch 2010 closed-loop natural frequency):
%   kR    ~ J * omega_n^2     omega_n: desired bandwidth
%   kOm   ~ 2 * zeta * omega_n * J
% Tailsitter J = diag([2200, 3800, 5500]) kg*m^2 (with J_xz cross term).
% Selection: omega_n = 5 rad/s, zeta = 0.7 (well-damped, stays inside 4*4500*5 N*m envelope).
omega_n  = 5.0;
zeta_att = 0.7;
J_diag = [2200; 3800; 5500];                    % must match aircraft_config J diagonal
ctrl.so3.kR    = diag(J_diag * omega_n^2);              % [55000, 95000, 137500]
ctrl.so3.kOm   = diag(2 * zeta_att * omega_n * J_diag); % [15400, 26600, 38500]
ctrl.so3.kI    = diag(0.05 * J_diag);                   % integral matched to inertia
ctrl.so3.I_max = deg2rad(15.0);                          % anti-windup saturation

% INDI parameters (Smeur 2016, Mancinelli 2024)
ctrl.indi.G       = [];      % effectiveness Jacobian (filled at init from aircraft model)
ctrl.indi.lp_omega = 50.0;   % [rad/s] low-pass on omega_dot measurement
ctrl.indi.lp_u    = 50.0;    % [rad/s] low-pass on actuator feedback
ctrl.indi.K_rate  = diag([20, 25, 12]); % rate-loop proportional gain
ctrl.indi.K_att   = diag([10, 12,  6]); % attitude-loop proportional gain (cascaded)

%% ----- Outer position/velocity loop -----
ctrl.outer.type = 'NMPC';   % nonlinear MPC

% NMPC parameters
ctrl.nmpc.N      = 20;            % horizon steps
ctrl.nmpc.dt     = 0.05;          % [s] horizon step (20 Hz)
ctrl.nmpc.Q_pos  = diag([10, 10, 25]);    % position tracking weight (Z heavier)
ctrl.nmpc.Q_vel  = diag([2, 2, 4]);
ctrl.nmpc.R_acc  = 0.05 * eye(3);         % accel command effort
ctrl.nmpc.Qf     = 50 * eye(6);           % terminal weight on [pos; vel]
ctrl.nmpc.acc_max  = 2.5 * 9.81;          % [m/s^2] max commanded specific force
ctrl.nmpc.tilt_max = deg2rad(85);         % max tilt of thrust vector from local-vertical
% Solver tolerances (fmincon)
ctrl.nmpc.opts = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display','off', ...
    'MaxIterations', 50, ...
    'OptimalityTolerance', 1e-4, ...
    'StepTolerance', 1e-6, ...
    'SpecifyObjectiveGradient', false);

%% ----- Control allocator -----
% Maps virtual controls [F_x; F_y; F_z; M_x; M_y; M_z] (body) -> u = [T_1..T_4; d_e; d_a; d_r]
ctrl.alloc.method = 'wpinv';   % weighted pseudoinverse with cant compensation
ctrl.alloc.W      = diag([1, 1, 1, 1, 5, 5, 5]); % effort weight (rotors cheaper than surfaces)
ctrl.alloc.T_min  = 0.0;
ctrl.alloc.T_max  = 4500.0;
ctrl.alloc.dT_max = 50000.0;   % [N/s] thrust rate limit (battery+ESC dynamics)

%% ----- Loop frequencies -----
ctrl.f_inner = 500.0;   % [Hz] inner loop (matches IMU bandwidth for INDI)
ctrl.f_outer = 50.0;    % [Hz] outer loop (NMPC at 20 Hz horizon)

end
