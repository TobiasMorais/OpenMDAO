function tr = test_rk4()
% TEST_RK4  Stage 5a — fixed-step RK4 integrator with quaternion renorm.
tr = test_helpers('init', 'test_rk4');

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'simulation'));
addpath(fullfile(here, '..', 'core'));

% R1: linear scalar dx/dt = -x converges with order 4
% Analytic: x(t) = x0 * exp(-t)
% RK4 should match closely for small dt
x0 = 1.0;
dt = 0.01;
T = 1.0;
N = round(T/dt);
x = x0;
for k = 1:N
    x = rk4_step(@(t, xx) -xx, k*dt, x, dt);
end
x_exact = x0 * exp(-T);
err1 = abs(x - x_exact);
fprintf('       R1 RK4 error on exp decay: %.3e\n', err1);
tr = test_helpers('assert_lt', tr, err1, 1e-8, 'R1 RK4 4th order accuracy');

% R2: quaternion renormalization via state with q in 7:10
% Use a fake odefun that drives a small drift
fake_ode = @(t, x) [zeros(6,1); 0.01; 0; 0; 0; zeros(3,1)];   % adds to q components
x = [zeros(6,1); 1; 0; 0; 0; zeros(3,1)];
for k = 1:1000
    x = rk4_step(fake_ode, k*0.01, x, 0.01);
end
qnorm_err = abs(norm(x(7:10)) - 1);
tr = test_helpers('assert_lt', tr, qnorm_err, 1e-9, 'R2 quaternion renormalized post-RK4');

% R3: harmonic oscillator energy approximately conserved
% dx/dt = v, dv/dt = -x. E = 0.5*(x^2+v^2)
x = 1.0; v = 0.0;
state = [x; v];
ode_h = @(t, s) [s(2); -s(1)];
dt = 0.001;
N = round(2*pi / dt);   % one period
E0 = 0.5 * (x^2 + v^2);
for k = 1:N
    state = rk4_step(ode_h, k*dt, state, dt);
end
E1 = 0.5 * (state(1)^2 + state(2)^2);
energy_err = abs(E1 - E0) / E0;
fprintf('       R3 harmonic oscillator energy err over 1 period: %.3e\n', energy_err);
tr = test_helpers('assert_lt', tr, energy_err, 1e-3, 'R3 RK4 energy preservation 1 period');

% R4: stability with bounded inputs — sin forced osc stays bounded
x = 0; v = 0;
state = [x; v];
ode_f = @(t, s) [s(2); -s(1) + sin(t)];
for k = 1:5000
    state = rk4_step(ode_f, k*0.01, state, 0.01);
end
tr = test_helpers('assert_lt', tr, max(abs(state)), 1e2, 'R4 forced osc stays bounded');

tr = test_helpers('report', tr);
end
