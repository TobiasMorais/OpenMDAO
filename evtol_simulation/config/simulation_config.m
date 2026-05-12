function sim = simulation_config()
% SIMULATION_CONFIG  Top-level simulation settings (solver, time, scenario, environment).

%% --- Numerical integration ---
sim.solver   = 'rk4';     % 'rk4' (fixed-step, HIL-ready) | 'ode45' (variable, accuracy)
sim.dt       = 0.002;     % [s] fixed-step (500 Hz, matches inner loop)
sim.t_final  = 600.0;     % [s] full mission with conservative trajectory (rest-to-rest only)

%% --- Default mission scenario ---
% Selectable: 'full_mission', 'transition_only', 'hover_disturbed'
sim.scenario = 'full_mission';

%% --- Environment ---
sim.env.gravity   = 9.80665;                  % [m/s^2]
sim.env.rho0      = 1.225;                    % [kg/m^3] sea-level density
sim.env.T0_K      = 288.15;                   % ISA standard
sim.env.h0        = 0.0;                      % [m] reference altitude

% Wind: constant + Dryden turbulence (MIL-F-8785C, low-altitude category)
sim.wind.constant_NED = [3.0; -1.0; 0.0];     % [m/s] mean wind in NED
sim.wind.dryden.enable = true;
sim.wind.dryden.W20    = 7.7;                 % [m/s] wind speed at 20 ft (light turbulence)
sim.wind.dryden.h_ft   = 100.0;               % altitude for L_w lookup
sim.wind.dryden.seed   = 42;

%% --- Sensors (for HIL realism) ---
sim.sensors.imu.bias_gyro  = deg2rad([0.05; -0.03; 0.04]);  % [rad/s]
sim.sensors.imu.noise_gyro = deg2rad(0.01);                  % stddev
sim.sensors.imu.bias_acc   = [0.02; -0.01; 0.03];            % [m/s^2]
sim.sensors.imu.noise_acc  = 0.05;
sim.sensors.imu.dt         = 0.002;                          % 500 Hz
sim.sensors.imu.enable_noise = true;

%% --- Logging ---
sim.log.fields = {'t','pos_NED','vel_NED','quat','omega_body', ...
                  'thrust_cmd','thrust_actual','elev','ail','rud', ...
                  'pitch_deg','alpha','beta','tracking_err'};
sim.log.downsample = 1;   % 1 = log every step

end
