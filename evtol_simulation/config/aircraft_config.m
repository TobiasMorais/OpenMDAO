function ac = aircraft_config()
% AIRCRAFT_CONFIG  Returns parameter struct for the 1500 kg biplace tailsitter eVTOL.
%
% Reference configuration: traditional fixed-wing layout (fuselage, main wing,
% conventional empennage) operating as a tailsitter (vertical takeoff/landing).
% Four fixed rotors are mounted in pairs on each wingtip with configurable
% pitch/yaw cant angles (geometric design parameters, swept across runs).
%
% All values follow Joby S4 / Beta Alia scaling, refs [2,17,21,46,55] of PDF.
% Units: SI (kg, m, s, rad, N).

%% --- Mass and inertia (Section "Modelagem de Massa, Inércia e Flexibilidade") ---
ac.mass            = 1500.0;                  % [kg] total mass (battery included, constant)
ac.cg_body         = [0.05; 0.0; 0.0];        % [m] CG offset from geometric ref in body axes

% Inertia tensor [kg*m^2], plane of symmetry X-Z (so J_xy = J_yz = 0)
% Non-zero J_xz captures longitudinal asymmetry of tailsitter (wing forward, tail aft, cockpit).
ac.J = [ 2200.0   0.0    180.0;
            0.0  3800.0    0.0;
          180.0   0.0   5500.0 ];   % [kg*m^2]
ac.J_inv = inv(ac.J);

%% --- Wing (main lifting surface) ---
ac.wing.span       = 11.0;                    % [m] wingspan b
ac.wing.chord      = 1.4;                     % [m] mean aerodynamic chord
ac.wing.area       = 15.0;                    % [m^2] reference area S
ac.wing.AR         = ac.wing.span^2 / ac.wing.area; % aspect ratio
ac.wing.e          = 0.85;                    % Oswald efficiency
ac.wing.airfoil    = 'NACA-0015';             % symmetric for bidirectional ops
ac.wing.alpha0     = 0.0;                     % [rad] zero-lift AoA
ac.wing.CL_alpha   = 5.7;                     % [1/rad] lift slope (linear regime)
ac.wing.CL_max     = 1.35;                    % stall CL
ac.wing.alpha_stall = deg2rad(15.0);          % [rad] stall onset
ac.wing.CD0        = 0.022;                   % parasitic drag
ac.wing.Cm0        = 0.0;                     % zero-AoA pitching moment coef
ac.wing.Cm_alpha   = -0.45;                   % static stability derivative

%% --- Horizontal stabilizer (elevator) ---
ac.htail.area      = 2.5;                     % [m^2]
ac.htail.arm       = 4.5;                     % [m] tail moment arm (CG to AC)
ac.htail.CL_alpha  = 4.5;                     % [1/rad]
ac.htail.eta       = 0.9;                     % dynamic pressure ratio
ac.htail.elev_eff  = 0.55;                    % elevator effectiveness (dCL/d_delta_e)

%% --- Vertical stabilizer (rudder) ---
ac.vtail.area      = 1.6;                     % [m^2]
ac.vtail.arm       = 4.5;                     % [m]
ac.vtail.CY_beta   = -0.55;                   % side force per sideslip
ac.vtail.rud_eff   = 0.45;                    % rudder effectiveness

%% --- Ailerons / elevons ---
ac.ail.eff         = 0.18;                    % rolling moment per unit deflection
ac.ail.delta_max   = deg2rad(25.0);

%% --- Control surface limits ---
ac.surf.elev_max   = deg2rad(25.0);
ac.surf.rud_max    = deg2rad(25.0);
ac.surf.tau        = 0.04;                    % [s] servo first-order time constant

%% --- Rotors (4 fixed thrusters, 2 per wingtip) ---
% Position vectors r_i [m] from CG to each rotor hub, body frame.
% Layout: wingtip-mounted, 2 stacked per side (one slightly forward and one aft of wing LE).
ac.n_rotors = 4;

% Rotor positions [m] (body frame: X fwd, Y right, Z down)
ac.rotor.position = [ ...
     0.30,  +5.0, -0.15;   % rotor 1: right wingtip, upper
     0.30,  +5.0, +0.15;   % rotor 2: right wingtip, lower
     0.30,  -5.0, -0.15;   % rotor 3: left  wingtip, upper
     0.30,  -5.0, +0.15];  % rotor 4: left  wingtip, lower

% Rotation direction (+1 CCW viewed from front, -1 CW). Pairs cancel reaction torque.
ac.rotor.spin_dir = [+1; -1; -1; +1];

% --- Cant angles (PARAMETRIC SWEEP — set per run) ---
% Each rotor thrust axis is rotated from +X_B by:
%   first by yaw_cant about Z_B, then by pitch_cant about Y_B
% Positive pitch_cant tilts thrust slightly downward (toward +Z_B) in hover.
% Positive yaw_cant tilts thrust outboard (toward +Y_B for right rotors).
ac.rotor.pitch_cant = deg2rad([+2.0; -2.0; +2.0; -2.0]);  % [rad]
ac.rotor.yaw_cant   = deg2rad([+1.0; +1.0; -1.0; -1.0]);  % [rad]

% Rotor performance
ac.rotor.radius        = 0.90;                % [m]
ac.rotor.disk_area     = pi * ac.rotor.radius^2;
ac.rotor.n_blades      = 5;
ac.rotor.blade_chord   = 0.10;                % [m]
ac.rotor.blade_pitch   = deg2rad(12.0);       % collective
ac.rotor.solidity      = ac.rotor.n_blades * ac.rotor.blade_chord ...
                         / (pi * ac.rotor.radius);
ac.rotor.cl_alpha_blade = 5.7;                % [1/rad] blade airfoil
ac.rotor.cd0_blade     = 0.012;
ac.rotor.J_rotor       = 0.32;                % [kg*m^2] each rotor moment of inertia
ac.rotor.omega_max     = 230.0;               % [rad/s] (~2200 rpm)
ac.rotor.thrust_max    = 4500.0;              % [N] per rotor (4*4500=18000 N > 1500*9.81)
ac.rotor.tau_motor     = 0.06;                % [s] ESC/motor first-order time constant
ac.rotor.kT            = ac.rotor.thrust_max / ac.rotor.omega_max^2; % static thrust coef
ac.rotor.kQ            = 0.025 * ac.rotor.kT * ac.rotor.radius;       % torque coef

% Slipstream geometry: which surfaces are bathed by which rotors
% 1 = bathed, 0 = not. Surfaces order: [wing_R, wing_L, htail, vtail]
ac.slipstream.bath_matrix = [ ...
    1, 0, 1, 0;   % rotor 1 (right upper) bathes right wing + htail
    1, 0, 1, 1;   % rotor 2 (right lower) bathes right wing + htail + vtail
    0, 1, 1, 0;   % rotor 3 (left upper)
    0, 1, 1, 1];  % rotor 4 (left lower)
ac.slipstream.k_s   = 1.8;                    % wake development factor (Patterson 2014)
ac.slipstream.x_off = 0.8;                    % [m] longitudinal distance rotor->surface

%% --- Energy storage (constant mass for eVTOL) ---
ac.battery.capacity_kWh = 280.0;              % nominal pack
ac.battery.mass_kg      = 420.0;              % included in ac.mass

end
