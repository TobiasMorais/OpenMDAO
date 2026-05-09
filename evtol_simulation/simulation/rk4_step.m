function x_next = rk4_step(odefun, t, x, dt)
% RK4_STEP  Classical fixed-step Runge-Kutta 4 integrator.
%
%   x_{n+1} = x_n + dt/6 * (k1 + 2 k2 + 2 k3 + k4)
%   k1 = f(t,         x)
%   k2 = f(t + dt/2,  x + dt/2 * k1)
%   k3 = f(t + dt/2,  x + dt/2 * k2)
%   k4 = f(t + dt,    x + dt   * k3)
%
% Quaternion components are renormalized after the step to maintain
% unit norm (drift control); columns 7:10 are assumed to be the quaternion.

    k1 = odefun(t,         x);
    k2 = odefun(t + dt/2,  x + dt/2 * k1);
    k3 = odefun(t + dt/2,  x + dt/2 * k2);
    k4 = odefun(t + dt,    x + dt   * k3);
    x_next = x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);

    if numel(x_next) >= 10
        q = x_next(7:10);
        nq = norm(q);
        if nq > 1e-9
            x_next(7:10) = q / nq;
        else
            x_next(7:10) = [1;0;0;0];
        end
    end
end
