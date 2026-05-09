classdef DrydenWind < handle
% DRYDENWIND  MIL-F-8785C Dryden turbulence model + steady mean wind.
%
% Linear filter formulation: drives unit-variance white noise n(t) through
% three first-order shaping filters (one per wind axis: u, v, w):
%
%   H_u(s) = sigma_u * sqrt(2*L_u/(pi*V)) / (1 + L_u/V * s)
%   H_v(s) = sigma_v * sqrt(L_v/(pi*V))  * (1 + sqrt(3)*L_v/V * s) / (1 + L_v/V * s)^2
%   H_w(s) = sigma_w * sqrt(L_w/(pi*V))  * (1 + sqrt(3)*L_w/V * s) / (1 + L_w/V * s)^2
%
% with low-altitude (h<1000ft) length scales and intensities:
%   L_w = h, L_u = L_v = h/(0.177+0.000823*h)^1.2
%   sigma_w = 0.1 * W20, sigma_u = sigma_v = sigma_w / (0.177+0.000823*h)^0.4
%
% W20 is wind speed at 20 ft altitude (default 7.7 m/s for light turbulence).
% This implementation uses Tustin-discretized 1st/2nd-order filters.

    properties
        cfg
        rng
        % filter states
        x_u  = 0;
        x_v  = [0; 0];
        x_w  = [0; 0];
        constant_NED
    end

    methods
        function obj = DrydenWind(cfg)
            obj.cfg = cfg;
            obj.rng = RandStream('mt19937ar', 'Seed', cfg.dryden.seed);
            obj.constant_NED = cfg.constant_NED;
        end

        function v_wind_NED = step(obj, V_airspeed, h_m, dt)
            % Returns 3x1 wind velocity in NED.
            if ~obj.cfg.dryden.enable
                v_wind_NED = obj.constant_NED;
                return;
            end
            h_ft = max(50, h_m * 3.281);
            V    = max(2.0, V_airspeed);
            W20  = obj.cfg.dryden.W20;

            denom = (0.177 + 0.000823 * h_ft).^1.2;
            L_w = h_ft;
            L_u = h_ft / denom;
            L_v = L_u;

            sigma_w = 0.1 * W20;
            sigma_u = sigma_w / (0.177 + 0.000823*h_ft)^0.4;
            sigma_v = sigma_u;

            % --- u channel (1st-order) ---
            tau_u = L_u / V;
            a_u = exp(-dt/tau_u);
            n1 = randn(obj.rng);
            obj.x_u = a_u * obj.x_u + sigma_u * sqrt(1 - a_u^2) * n1;
            u_g = obj.x_u;

            % --- v channel (2nd-order) ---
            tau_v = L_v / V;
            % Discrete state-space (continuous: 1/(1+tau s)^2 with zero (1+sqrt(3)tau s))
            a_v = exp(-dt/tau_v);
            n2 = randn(obj.rng);
            new_v1 = a_v * obj.x_v(1) + sigma_v * sqrt(1 - a_v^2) * n2;
            new_v2 = a_v * obj.x_v(2) + (1 - a_v) * (new_v1 + sqrt(3) * tau_v * (new_v1 - obj.x_v(1))/dt);
            obj.x_v = [new_v1; new_v2];
            v_g = obj.x_v(2);

            % --- w channel (2nd-order) ---
            tau_w = L_w / V;
            a_w = exp(-dt/tau_w);
            n3 = randn(obj.rng);
            new_w1 = a_w * obj.x_w(1) + sigma_w * sqrt(1 - a_w^2) * n3;
            new_w2 = a_w * obj.x_w(2) + (1 - a_w) * (new_w1 + sqrt(3) * tau_w * (new_w1 - obj.x_w(1))/dt);
            obj.x_w = [new_w1; new_w2];
            w_g = obj.x_w(2);

            v_wind_NED = obj.constant_NED + [u_g; v_g; w_g];
        end

        function reset(obj)
            obj.x_u = 0;
            obj.x_v = [0; 0];
            obj.x_w = [0; 0];
            obj.rng.reset();
        end
    end
end
