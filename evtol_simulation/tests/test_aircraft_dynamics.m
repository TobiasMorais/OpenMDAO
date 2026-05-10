function tr = test_aircraft_dynamics()
% TEST_AIRCRAFT_DYNAMICS  Stage 1c — 6-DOF rigid-body dynamics validation.
tr = test_helpers('init', 'test_aircraft_dynamics');

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'core'));
addpath(fullfile(here, '..', 'config'));
addpath(fullfile(here, '..', 'simulation'));

ac_cfg = aircraft_config();
ac     = Aircraft(ac_cfg);

% D1: Free fall — no F, no M, identity orientation
x0 = [0;0;0;  0;0;0;  1;0;0;0;  0;0;0];
xdot = ac.dynamics(x0, zeros(3,1), zeros(3,1), zeros(3,1), zeros(3,1), zeros(3,1));
tr = test_helpers('assert_near', tr, xdot(4:6), [0;0;ac.gravity], 1e-12, 'D1 free fall: a_NED = g');

% D2: Quaternion norm preserved through 1000 RK4 steps in hover
x = ac.initial_state('hover');
% Apply just enough thrust to balance gravity (vertical tailsitter)
F_prop_B = [ac_cfg.mass * ac.gravity; 0; 0];   % thrust along +X_B (which points up)
M_prop_B = zeros(3,1);
F_aero = zeros(3,1); M_aero = zeros(3,1); M_gyro = zeros(3,1);
odefun = @(t,xx) ac.dynamics(xx, F_aero, M_aero, F_prop_B, M_prop_B, M_gyro);
dt = 0.002;
for k = 1:1000
    x = rk4_step(odefun, k*dt, x, dt);
end
qnorm_err = abs(norm(x(7:10)) - 1);
tr = test_helpers('assert_lt', tr, qnorm_err, 1e-6, 'D2 quaternion norm preserved over 1000 RK4 steps');

% D3: Hover trim — verify acceleration is near zero
x_h = ac.initial_state('hover');
xdot_h = ac.dynamics(x_h, zeros(3,1), zeros(3,1), F_prop_B, zeros(3,1), zeros(3,1));
% Note: at hover state q points x_B up, so R*F_prop_B should yield -m*g vector (NED up)
tr = test_helpers('assert_lt', tr, norm(xdot_h(4:6)), 1e-9, 'D3 hover trim: a_NED ~ 0');

% D4: Pure moment on M_x => angular accel = J^-1 [1;0;0]
M_test = [1; 0; 0];
xdot4 = ac.dynamics(x_h, zeros(3,1), M_test, zeros(3,1), zeros(3,1), zeros(3,1));
expected_wdot = ac_cfg.J_inv * M_test;
tr = test_helpers('assert_near', tr, xdot4(11:13), expected_wdot, 1e-9, 'D4 J^-1 * M_x rate');

% D5: With omega only, no torque, position must remain fixed (only attitude rotates)
x5 = x_h;
x5(11:13) = [0; 0.5; 0];
xdot5 = ac.dynamics(x5, zeros(3,1), zeros(3,1), F_prop_B, zeros(3,1), zeros(3,1));
% velocity NED should still be ~0
tr = test_helpers('assert_lt', tr, norm(xdot5(4:6)), 1e-6, 'D5 no force => v_dot ~ 0 in hover');

% D6: J_xz coupling — omega_x produces nonzero alpha_z component
x6 = x_h;
x6(11:13) = [1.0; 0; 0];
% Inertial torque: -omega x J*omega = -[1;0;0] x J*[1;0;0]
% J*[1;0;0] = [J_xx; 0; -J_xz]
% [1;0;0] x [J_xx; 0; -J_xz] = [0*(-J_xz) - 0*0; 0*J_xx - 1*(-J_xz); 1*0 - 0*J_xx] = [0; J_xz; 0]
% so the inertial cross term contributes -[0; J_xz; 0] to RHS, then J^-1 of that...
xdot6 = ac.dynamics(x6, zeros(3,1), zeros(3,1), zeros(3,1), zeros(3,1), zeros(3,1));
% wdot expected = -J^-1 * (omega x J*omega)
inertial_term = cross([1;0;0], ac_cfg.J*[1;0;0]);   % [0; -J_xz; 0]... wait
% The dynamics returns J^-1 * (M_total - omega x J*omega) where M_total=0
expected = -ac_cfg.J_inv * inertial_term;
tr = test_helpers('assert_near', tr, xdot6(11:13), expected, 1e-9, 'D6 J_xz inertial coupling');

% D7: Rotor gyroscopic torque — rotors spinning, body angular rate produces precession
Omega_signed = [200; -200; -200; 200];   % CCW/CW pairs
w_B = [0; 1; 0];
M_gyro = ac.compute_rotor_gyro(Omega_signed, w_B);
% Sum of angular momenta H_rotor = J_r * sum(Omega_signed_i * e_i)
% e_i are nominally +X_B with small cants. Sum ~ J_r * (200-200-200+200) * X_B = 0!
% This means with paired CCW/CW the net rotor H is zero -> M_gyro = 0
% Verify zero or near-zero
tr = test_helpers('assert_lt', tr, norm(M_gyro), 1e-3, 'D7 paired rotors cancel net H => M_gyro ~ 0');

tr = test_helpers('report', tr);
end
