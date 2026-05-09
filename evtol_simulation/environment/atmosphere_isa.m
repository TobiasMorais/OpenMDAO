function rho = atmosphere_isa(h_m, sim_env)
% ATMOSPHERE_ISA  Density via ISA standard atmosphere (troposphere only).
%
%   T(h) = T0 - L * h        (L = 0.0065 K/m)
%   p(h) = p0 * (T/T0)^(g/(R*L))
%   rho  = p / (R * T)

    if nargin < 2
        T0 = 288.15;
        rho0 = 1.225;
    else
        T0 = sim_env.T0_K;
        rho0 = sim_env.rho0;
    end
    L  = 0.0065;
    g  = 9.80665;
    R  = 287.058;
    h  = max(0, min(11000, h_m));
    T  = T0 - L*h;
    rho = rho0 * (T/T0)^(g/(R*L) - 1);
end
