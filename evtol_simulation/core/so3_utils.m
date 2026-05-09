function out = so3_utils(action, varargin)
% SO3_UTILS  Lie-algebra helpers on SO(3).
%
% Operations:
%   X    = so3_utils('hat', v)          % skew-symmetric [v]_x
%   v    = so3_utils('vee', X)          % inverse hat
%   eR   = so3_utils('errMat', R, Rd)   % geometric attitude error eR = 0.5*vee(Rd^T R - R^T Rd)
%   R    = so3_utils('exp', phi)        % so(3) -> SO(3) Rodrigues
%   phi  = so3_utils('log', R)          % SO(3) -> so(3)

switch action
    case 'hat'
        v = varargin{1};
        out = [   0, -v(3),  v(2);
              v(3),    0, -v(1);
             -v(2),  v(1),    0];
    case 'vee'
        X = varargin{1};
        out = [X(3,2); X(1,3); X(2,1)];
    case 'errMat'
        % Lee, Leok, McClamroch (2010) Eq. (10): e_R = 0.5 vee(Rd^T R - R^T Rd)
        R  = varargin{1};
        Rd = varargin{2};
        S  = 0.5 * (Rd' * R - R' * Rd);
        out = so3_utils('vee', S);
    case 'exp'
        phi = varargin{1};
        a = norm(phi);
        if a < 1e-9
            out = eye(3) + so3_utils('hat', phi);
        else
            ax = phi/a;
            K  = so3_utils('hat', ax);
            out = eye(3) + sin(a)*K + (1-cos(a))*(K*K);
        end
    case 'log'
        R = varargin{1};
        c = max(-1, min(1, (trace(R)-1)/2));
        a = acos(c);
        if a < 1e-9
            out = 0.5 * so3_utils('vee', R - R');
        else
            out = (a / (2*sin(a))) * so3_utils('vee', R - R');
        end
    otherwise
        error('so3_utils: unknown action "%s"', action);
end
end
