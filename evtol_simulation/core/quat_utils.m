function out = quat_utils(action, varargin)
% QUAT_UTILS  Hamilton-convention unit-quaternion utilities (q = [q0; qv]).
%
% Singularity-free attitude representation used throughout the simulator.
% Hamilton convention (right-handed, scalar-first):  q = [q0; q1; q2; q3]
% with rotation matrix R(q) = (q0^2 - qv'*qv)*I + 2*qv*qv' - 2*q0*[qv]_x
% mapping body-to-world (R * v_body = v_world) -- matches PDF NED frame.
%
% Operations:
%   q  = quat_utils('mul', q1, q2)         % Hamilton product
%   qc = quat_utils('conj', q)
%   q  = quat_utils('norm', q)             % normalize
%   R  = quat_utils('toR', q)              % body-to-world rotation matrix
%   q  = quat_utils('fromR', R)
%   q  = quat_utils('fromEuler', phi,th,psi)   % Z-Y-X Tait-Bryan
%   eu = quat_utils('toEuler', q)              % [phi; theta; psi]
%   qd = quat_utils('kinematic', q, omega_body) % q-dot from body rates
%   qe = quat_utils('errMul', q, qd)            % multiplicative error qe = q^-1 * qd
%   v  = quat_utils('rotate', q, v)             % rotate vector by q
%   q  = quat_utils('expMap', phi)               % so(3) -> SU(2): exp(phi/2)
%   phi = quat_utils('logMap', q)                % SU(2) -> so(3): 2*log(q)

switch action
    case 'mul'
        a = varargin{1}; b = varargin{2};
        out = [a(1)*b(1) - a(2:4).'*b(2:4);
               a(1)*b(2:4) + b(1)*a(2:4) + cross(a(2:4), b(2:4))];
    case 'conj'
        q = varargin{1};
        out = [q(1); -q(2:4)];
    case 'norm'
        q = varargin{1};
        n = norm(q);
        if n < 1e-12
            out = [1;0;0;0];
        else
            out = q / n;
        end
    case 'toR'
        q = varargin{1};
        q0 = q(1); qx = q(2); qy = q(3); qz = q(4);
        out = [1-2*(qy^2+qz^2),  2*(qx*qy - qz*q0), 2*(qx*qz + qy*q0);
               2*(qx*qy + qz*q0), 1-2*(qx^2+qz^2),  2*(qy*qz - qx*q0);
               2*(qx*qz - qy*q0), 2*(qy*qz + qx*q0), 1-2*(qx^2+qy^2)];
    case 'fromR'
        R = varargin{1};
        tr = trace(R);
        if tr > 0
            S = 2*sqrt(tr + 1);
            q0 = 0.25*S;
            qx = (R(3,2)-R(2,3))/S;
            qy = (R(1,3)-R(3,1))/S;
            qz = (R(2,1)-R(1,2))/S;
        elseif (R(1,1) > R(2,2)) && (R(1,1) > R(3,3))
            S = 2*sqrt(1 + R(1,1) - R(2,2) - R(3,3));
            q0 = (R(3,2)-R(2,3))/S;
            qx = 0.25*S;
            qy = (R(1,2)+R(2,1))/S;
            qz = (R(1,3)+R(3,1))/S;
        elseif R(2,2) > R(3,3)
            S = 2*sqrt(1 + R(2,2) - R(1,1) - R(3,3));
            q0 = (R(1,3)-R(3,1))/S;
            qx = (R(1,2)+R(2,1))/S;
            qy = 0.25*S;
            qz = (R(2,3)+R(3,2))/S;
        else
            S = 2*sqrt(1 + R(3,3) - R(1,1) - R(2,2));
            q0 = (R(2,1)-R(1,2))/S;
            qx = (R(1,3)+R(3,1))/S;
            qy = (R(2,3)+R(3,2))/S;
            qz = 0.25*S;
        end
        out = quat_utils('norm', [q0; qx; qy; qz]);
    case 'fromEuler'
        phi = varargin{1}; th = varargin{2}; psi = varargin{3};
        cphi = cos(phi/2); sphi = sin(phi/2);
        cth  = cos(th/2);  sth  = sin(th/2);
        cpsi = cos(psi/2); spsi = sin(psi/2);
        out = [ cphi*cth*cpsi + sphi*sth*spsi;
                sphi*cth*cpsi - cphi*sth*spsi;
                cphi*sth*cpsi + sphi*cth*spsi;
                cphi*cth*spsi - sphi*sth*cpsi];
    case 'toEuler'
        q = varargin{1};
        q0=q(1); qx=q(2); qy=q(3); qz=q(4);
        phi   = atan2(2*(q0*qx + qy*qz), 1 - 2*(qx^2 + qy^2));
        sinp  = max(-1, min(1, 2*(q0*qy - qz*qx)));
        theta = asin(sinp);
        psi   = atan2(2*(q0*qz + qx*qy), 1 - 2*(qy^2 + qz^2));
        out = [phi; theta; psi];
    case 'kinematic'
        % q_dot = 0.5 * Omega(omega) * q
        q = varargin{1}; w = varargin{2};
        Omega = [ 0,    -w(1), -w(2), -w(3);
                  w(1),   0,    w(3), -w(2);
                  w(2), -w(3),   0,    w(1);
                  w(3),  w(2), -w(1),   0  ];
        out = 0.5 * Omega * q;
    case 'errMul'
        % Multiplicative error: qe = qd^-1 (x) q  (so qe -> [1;0;0;0] when q=qd)
        q  = varargin{1};
        qd = varargin{2};
        out = quat_utils('mul', quat_utils('conj', qd), q);
        if out(1) < 0,  out = -out;  end   % shortest-path / unwinding fix
    case 'rotate'
        q = varargin{1}; v = varargin{2};
        out = quat_utils('toR', q) * v;
    case 'expMap'
        phi = varargin{1};
        a = norm(phi);
        if a < 1e-9
            out = [1; 0.5*phi];
        else
            out = [cos(a/2); sin(a/2)/a * phi];
        end
    case 'logMap'
        q = varargin{1};
        if q(1) < 0, q = -q; end
        v = q(2:4); s = q(1);
        nv = norm(v);
        if nv < 1e-9
            out = 2*v;
        else
            out = 2*atan2(nv, s) * (v / nv);
        end
    otherwise
        error('quat_utils: unknown action "%s"', action);
end
end
