function [qd, T_total, Wd, Wd_dot] = force_to_attitude(F_des_NED, psi_d, dt, prev)
% FORCE_TO_ATTITUDE  Map desired NED specific force to (T_total, qd, Wd, Wd_dot).
%
% Construction (SE(3)-consistent for tailsitter):
%   x_B = F_des_NED / |F_des_NED|             body x = thrust direction
%   c_psi = [cos(psi); sin(psi); 0]            heading reference vector in NED ground plane
%   y_B = (c_psi x x_B) / |c_psi x x_B|        body y in plane perp to x_B and to c_psi
%   z_B = x_B x y_B                            body z completes right-handed frame
%
% Sanity checks:
%   - Hover (F = -m g e_z, psi=0): x_B=[0;0;-1], c_psi=[1;0;0]
%       y_B = cross([1;0;0],[0;0;-1]) / |.| = [0;1;0]   (east)
%       z_B = cross([0;0;-1],[0;1;0]) = [1;0;0]         (north)
%     This corresponds to tailsitter pointing up with belly facing north.
%   - Hover psi=90deg: c_psi=[0;1;0]
%       y_B = cross([0;1;0],[0;0;-1]) = [-1;0;0]
%       z_B = cross([0;0;-1],[-1;0;0]) = [0;1;0]
%     Belly now faces east — heading rotates body about thrust axis. PASS F6b.
%
% Singularity: when c_psi parallel to x_B (e.g., level cruise heading north
% with thrust horizontal-northward), use altitude reference fallback.
%
% Wd, Wd_dot estimated by finite difference between successive calls.

    Fmag = norm(F_des_NED);
    if Fmag < 1e-3
        qd = quat_utils('fromEuler', 0, deg2rad(90), psi_d);
        T_total = 0;
        Wd = zeros(3,1);
        Wd_dot = zeros(3,1);
        return;
    end

    x_B = F_des_NED / Fmag;
    c_psi = [cos(psi_d); sin(psi_d); 0];

    % Build y_B as cross(c_psi, x_B)
    y_B_unnorm = cross(c_psi, x_B);
    if norm(y_B_unnorm) < 1e-3
        % Singular: thrust direction parallel to heading vector. Use down-fallback.
        % Pick reference perpendicular to x_B.
        if abs(x_B(3)) < 0.99
            y_B_unnorm = cross([0;0;1], x_B);
        else
            y_B_unnorm = cross([0;1;0], x_B);
        end
    end
    y_B = y_B_unnorm / norm(y_B_unnorm);
    z_B = cross(x_B, y_B);
    z_B = z_B / norm(z_B);   % numerical safety

    R = [x_B, y_B, z_B];
    % Right-handed safety
    if det(R) < 0
        z_B = -z_B;
        R = [x_B, y_B, z_B];
    end

    qd = quat_utils('fromR', R);
    T_total = Fmag;

    % Estimate Wd, Wd_dot
    if isempty(prev) || ~isfield(prev,'qd') || isempty(prev.qd) || dt <= 0
        Wd = zeros(3,1);
        Wd_dot = zeros(3,1);
    else
        qrel = quat_utils('mul', quat_utils('conj', prev.qd), qd);
        if qrel(1) < 0, qrel = -qrel; end
        phi = quat_utils('logMap', qrel);
        Wd = phi / dt;

        if isfield(prev,'Wd') && ~isempty(prev.Wd)
            Wd_dot = (Wd - prev.Wd) / dt;
            % Saturate Wd_dot to avoid spikes from FD noise
            Wd_dot = max(-50, min(50, Wd_dot));
        else
            Wd_dot = zeros(3,1);
        end
    end
end
