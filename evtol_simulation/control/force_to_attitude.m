function [qd, T_total, Wd, Wd_dot] = force_to_attitude(F_des_NED, psi_d, dt, prev)
% FORCE_TO_ATTITUDE  Map desired NED force vector to (T_total, qd, Wd, Wd_dot).
%
% Differential-flatness inspired mapping (Mellinger & Kumar 2011 generalized):
%   - T points along -z_B in conventional quad; for a tailsitter we choose
%     thrust along +x_B (since rotor axes are nominally +X_B), so:
%
%       x_B = F_des_NED / |F_des_NED|        (desired body-X = thrust axis)
%       z_B chosen to maximize alignment with NED -Z (or with cruise direction
%            via psi_d when in level flight)
%       y_B = z_B x x_B    then re-orthogonalize z_B = x_B x y_B
%
% Heading is parameterized by psi_d (a flat output); during hover psi_d
% rotates the wing about thrust axis; during cruise it acts as conventional yaw.
%
% Wd, Wd_dot are estimated by finite difference between successive calls.
%
% Inputs
%   F_des_NED: 3x1 desired specific force in NED (gravity-compensated, * mass = thrust vector)
%   psi_d:     scalar desired heading (rad)
%   dt:        timestep since last call (for Wd estimation)
%   prev:      struct with prev.qd, prev.Wd  (use [] on first call)
%
% Outputs
%   qd:      desired quaternion (body->world)
%   T_total: total thrust magnitude commanded (scalar, sum of 4 rotors nominal)
%   Wd:      desired body angular rate (body-frame)
%   Wd_dot:  desired body angular accel (body-frame)

    Fmag = norm(F_des_NED);
    if Fmag < 1e-3
        % Degenerate: fall back to upright
        qd = quat_utils('fromEuler', 0, deg2rad(90), psi_d);
        T_total = 0;
        Wd = zeros(3,1);
        Wd_dot = zeros(3,1);
        return;
    end

    x_B = F_des_NED / Fmag;

    % Heading reference vector in NED ground plane
    h_ref = [cos(psi_d); sin(psi_d); 0];

    % y_B perpendicular to x_B and to h_ref
    y_B = cross([0;0;-1] - dot([0;0;-1], x_B)*x_B, x_B);
    if norm(y_B) < 1e-3
        % x_B nearly parallel to vertical -> use heading reference instead
        y_B = cross(h_ref, x_B);
    end
    y_B = y_B / max(norm(y_B), 1e-6);

    % Apply heading rotation about x_B (yaw in tailsitter coordinate)
    R0 = [x_B, y_B, cross(x_B, y_B)];
    Rpsi = [1, 0, 0; 0, cos(psi_d), -sin(psi_d); 0, sin(psi_d), cos(psi_d)];
    R = R0 * Rpsi;

    % Re-orthogonalize via SVD-free correction
    [U,~,V] = svd(R);
    R = U * V';
    if det(R) < 0
        R(:,3) = -R(:,3);
    end

    qd = quat_utils('fromR', R);
    T_total = Fmag;     % * mass  is applied in the caller; here we return specific force magnitude

    % Estimate Wd, Wd_dot by finite differences on prev.qd
    if isempty(prev) || ~isfield(prev,'qd') || isempty(prev.qd) || dt <= 0
        Wd = zeros(3,1);
        Wd_dot = zeros(3,1);
    else
        % Relative rotation Rprev -> R, then log map gives angle*axis
        qrel = quat_utils('mul', quat_utils('conj', prev.qd), qd);
        if qrel(1) < 0, qrel = -qrel; end
        phi = quat_utils('logMap', qrel);
        Wd = phi / dt;

        if isfield(prev,'Wd') && ~isempty(prev.Wd)
            Wd_dot = (Wd - prev.Wd) / dt;
            % Bound Wd_dot to avoid spikes
            Wd_dot = max(-50, min(50, Wd_dot));
        else
            Wd_dot = zeros(3,1);
        end
    end
end
