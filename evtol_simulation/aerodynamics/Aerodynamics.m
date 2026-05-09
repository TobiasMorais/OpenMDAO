classdef Aerodynamics < handle
% AERODYNAMICS  Computes F_aero^B and M_aero^B with 360-deg AoA support.
%
% Strip-theory decomposition of fuselage into surfaces: right wing, left wing,
% horizontal tail, vertical tail. Each strip computes local airspeed including
% (a) free-stream velocity in body frame V_inf,
% (b) propeller slipstream contribution V_slip from any rotor that bathes it,
% (c) angular-rate-induced velocity (omega x r),
% then evaluates lift and drag using a piecewise model:
%
%   - Linear regime |alpha| < alpha_stall:  C_L = C_Lalpha * (alpha - alpha0)
%   - Post-stall via Viterna-Corrigan extrapolation up to +/- 90 deg:
%
%     C_L(alpha) = A1 sin(2 alpha) + A2 cos(alpha)^2 / sin(alpha)
%     C_D(alpha) = B1 sin(alpha)^2 + B2 cos(alpha)
%
%     where A1 = CD_max/2, A2 = (CL_s - CD_max sin(as) cos(as)) sin(as)/cos(as)^2
%           B1 = CD_max,  B2 = (CD_s - CD_max sin(as)^2)/cos(as)
%     (Viterna & Corrigan, 1981)
%
% Strip force per surface in stability axes:
%   L = q_local * S_strip * C_L(alpha)
%   D = q_local * S_strip * C_D(alpha)
% rotated back to body frame, then summed; moment = r_strip x F_strip.

    properties
        cfg
        viterna_const   % cached Viterna coefficients per surface
    end

    methods
        function obj = Aerodynamics(cfg)
            obj.cfg = cfg;
            obj.viterna_const = obj.precompute_viterna();
        end

        function vc = precompute_viterna(obj)
            % Compute Viterna A1,A2,B1,B2 for each lifting surface.
            CD_max = 1.11 + 0.018 * obj.cfg.wing.AR;   % flat-plate at 90 deg
            as     = obj.cfg.wing.alpha_stall;
            CL_s   = obj.cfg.wing.CL_alpha * (as - obj.cfg.wing.alpha0);
            CD_s   = obj.cfg.wing.CD0 + (CL_s^2)/(pi * obj.cfg.wing.AR * obj.cfg.wing.e);

            A1 = CD_max / 2;
            A2 = (CL_s - CD_max*sin(as)*cos(as)) * sin(as)/cos(as)^2;
            B1 = CD_max;
            B2 = (CD_s - CD_max*sin(as)^2)/cos(as);

            vc.A1 = A1;  vc.A2 = A2;  vc.B1 = B1;  vc.B2 = B2;
            vc.CD_max = CD_max;  vc.alpha_stall = as;
        end

        function [CL, CD] = airfoil_360(obj, alpha)
            % Piecewise CL,CD with smooth transitions; alpha in [-pi, pi].
            alpha = atan2(sin(alpha), cos(alpha));   % wrap
            as = obj.viterna_const.alpha_stall;

            if abs(alpha) <= as
                CL = obj.cfg.wing.CL_alpha * (alpha - obj.cfg.wing.alpha0);
                CD = obj.cfg.wing.CD0 + (CL^2)/(pi*obj.cfg.wing.AR*obj.cfg.wing.e);
            else
                a_abs = abs(alpha);
                if a_abs > pi/2
                    a_abs = pi - a_abs;   % mirror for back-side flow
                end
                vc = obj.viterna_const;
                CL_v = vc.A1 * sin(2*a_abs) + vc.A2 * cos(a_abs)^2 / max(sin(a_abs), 1e-3);
                CD_v = vc.B1 * sin(a_abs)^2 + vc.B2 * cos(a_abs);
                CL = sign(alpha) * CL_v;
                CD = max(0.01, CD_v);
            end
        end

        function [F_B, M_B] = compute_forces_moments(obj, V_B, w_B, surf_def, V_slip_per_surf, rho)
            % V_B: 3x1 body-frame airspeed (relative to wind)
            % w_B: 3x1 body angular rates
            % surf_def: struct array with .position [3x1], .area, .role  per strip
            % V_slip_per_surf: 3 x N_surfaces additional velocity from slipstream
            % rho: air density
            F_B = zeros(3,1);
            M_B = zeros(3,1);

            for k = 1:numel(surf_def)
                r  = surf_def(k).position;
                S  = surf_def(k).area;
                role = surf_def(k).role;

                V_local = V_B + cross(w_B, r) + V_slip_per_surf(:, k);
                V_mag   = norm(V_local);
                if V_mag < 0.5
                    continue;
                end

                % AoA in surface plane; for a wing strip, alpha = atan2(-w, u)
                u = V_local(1); w = V_local(3);
                alpha = atan2(-w, u);                 % rad
                beta  = asin(max(-1,min(1, V_local(2)/V_mag)));   % side-slip

                [CL, CD] = obj.airfoil_360(alpha);
                qbar = 0.5 * rho * V_mag^2;
                L = qbar * S * CL;
                D = qbar * S * CD;

                % Surface contribution: lift normal to V_local in vertical plane,
                % drag along -V_local. For vertical tail rotate to xy-plane.
                switch role
                    case {'wing_R','wing_L','htail'}
                        e_drag = -V_local / V_mag;
                        % lift vector: rotate -V projection in x-z plane by +90 (perp to V, up)
                        Vxz = [V_local(1); 0; V_local(3)];
                        nVxz = max(norm(Vxz), 1e-6);
                        e_lift = [-Vxz(3); 0; Vxz(1)] / nVxz;   % perp to V_local in x-z, +z_lift
                        F_strip = L * e_lift + D * e_drag;
                    case 'vtail'
                        % Side force from beta primarily
                        Y = qbar * S * obj.cfg.vtail.CY_beta * beta;
                        F_strip = [-D; Y; 0];
                end

                F_B = F_B + F_strip;
                M_B = M_B + cross(r, F_strip);
            end

            % Static pitching moment about CG (wing camber + tail)
            V_inf_mag = norm(V_B);
            if V_inf_mag > 0.5
                qbar = 0.5 * rho * V_inf_mag^2;
                u = V_B(1); w = V_B(3);
                alpha = atan2(-w, u);
                Cm = obj.cfg.wing.Cm0 + obj.cfg.wing.Cm_alpha * alpha;
                M_B(2) = M_B(2) + qbar * obj.cfg.wing.area * obj.cfg.wing.chord * Cm;
            end
        end

        function defs = build_surface_strips(obj)
            % Returns struct array of surface strips for compute_forces_moments.
            ac = obj.cfg;
            defs(1) = struct('position', [0.05; +ac.wing.span/4; 0], ...
                             'area', ac.wing.area/2, 'role','wing_R');
            defs(2) = struct('position', [0.05; -ac.wing.span/4; 0], ...
                             'area', ac.wing.area/2, 'role','wing_L');
            defs(3) = struct('position', [-ac.htail.arm; 0; 0], ...
                             'area', ac.htail.area, 'role','htail');
            defs(4) = struct('position', [-ac.vtail.arm; 0; -0.5], ...
                             'area', ac.vtail.area, 'role','vtail');
        end
    end
end
