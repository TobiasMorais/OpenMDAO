% MAIN  Entry point for the eVTOL Tailsitter Simulation Framework.
%
% Plataforma de Simulação Virtual para Veículos Aéreos Não Convencionais
% Tailsitter biplace 1500 kg | 4 rotores fixos com cant configurável
%
% MISSION (realistic, 4-phase):
%   1. Vertical takeoff and climb to 20 m altitude
%   2. Hover hold at 20 m for 20 s
%   3. Climb-and-transition to 1000 m altitude (rest -> cruise)
%   4. Steady cruise at 1000 m, north direction, V at max L/D efficiency
%
% Pipeline:
%   1) Carrega configurações
%   2) Constrói missão realista (MissionTrajectory) com cálculos físicos
%   3) Executa simulação fechada (RK4 fixed-step, 500 Hz, Dryden)
%   4) Plota telemetria de Grau Aeroespacial
%
% Toggles:
%   ctrl.inner.type   = 'SO3' ou 'INDI'
%   MISSION           = 'realistic' (default) | 'legacy_full' | 'hover'

clear; clc; close all;

%% --- Add paths ---
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'config'));
addpath(fullfile(here, 'core'));
addpath(fullfile(here, 'aerodynamics'));
addpath(fullfile(here, 'control'));
addpath(fullfile(here, 'trajectory'));
addpath(fullfile(here, 'environment'));
addpath(fullfile(here, 'simulation'));
addpath(fullfile(here, 'visualization'));

%% --- Toggles ---
INNER_LAW = 'SO3';        % 'SO3' or 'INDI'
MISSION   = 'realistic';  % 'realistic' | 'legacy_full' | 'hover'

%% --- Load configurations ---
ac   = aircraft_config();
ctrl = controller_config();
sim  = simulation_config();
ctrl.inner.type = INNER_LAW;

fprintf('=== eVTOL Tailsitter Simulation Framework ===\n');
fprintf('Mass: %.0f kg | Rotors: %d | Inner: %s | Mission: %s\n', ...
    ac.mass, ac.n_rotors, ctrl.inner.type, MISSION);

%% --- Mission and trajectory ---
switch lower(MISSION)
    case 'realistic'
        % Conservative 4-phase mission (DifferentialFlatness engine).
        % All rest-to-rest, 30% envelope, 2x safety on time.
        [traj, mission_info] = build_realistic_mission(ac);
        sim.t_final = traj.total_time() + 10;
        sim.init_mode = mission_info.init_mode;   % start at -1m

    case 'legacy_full'
        waypoints = build_mission('full_mission');
        opt = TrajectoryOptimizer(waypoints, ac);
        traj = opt.optimize('method', 'flat');
        sim.t_final = traj.total_time() + 5;

    case 'hover'
        waypoints = build_mission('hover_only');
        opt = TrajectoryOptimizer(waypoints, ac);
        traj = opt.optimize('method', 'flat');
        sim.t_final = 30;
        sim.init_mode = 'hover_at_altitude';

    otherwise
        error('Unknown mission: %s', MISSION);
end
fprintf('Trajectory total time: %.2f s\n', traj.total_time());
fprintf('Simulation t_final:   %.2f s\n', sim.t_final);

%% --- Simulate ---
log = run_simulation(ac, ctrl, sim, traj);

%% --- Visualize ---
plot_telemetry(log, ac);

%% --- Physical validation against analytical predictions ---
addpath(fullfile(here, 'tests'));
if exist('mission_info', 'var')
    validate_physics(log, ac, mission_info);
else
    validate_physics(log, ac, struct());
end

%% --- Final report ---
err_norm = vecnorm(log.tracking_err, 2, 1);
fprintf('\n========================================================\n');
fprintf('   Mission Summary\n');
fprintf('========================================================\n');
fprintf('   Final position error    : %.2f m\n', err_norm(end));
fprintf('   Max position error      : %.2f m\n', max(err_norm));
fprintf('   Mean position error     : %.2f m\n', mean(err_norm));
fprintf('   Total energy proxy J    : %.3e (int T^1.5 dt)\n', ...
    trapz(log.t, log.thrust_total .^ 1.5));
fprintf('   Final altitude          : %.1f m\n', -log.pos_NED(3, end));
fprintf('   Final position N/E      : %.1f / %.1f m\n', ...
    log.pos_NED(1, end), log.pos_NED(2, end));
fprintf('========================================================\n');

fprintf('\nDone. Inspect figures for validation.\n');
