% MAIN  Entry point for the eVTOL Tailsitter Simulation Framework.
%
% Plataforma de Simulação Virtual para Veículos Aéreos Não Convencionais
% Tailsitter biplace 1500 kg | 4 rotores fixos com cant configurável
%
% Pipeline:
%   1) Carrega configurações
%   2) Constrói missão (waypoints) e otimiza trajetória (Flatness + GA opcional)
%   3) Executa simulação fechada (RK4 fixed-step, 500 Hz, Dryden)
%   4) Plota telemetria de Grau Aeroespacial
%
% Uso:
%   >> cd evtol_simulation
%   >> main
%
% Toggles principais (edite as linhas abaixo):
%   ctrl.inner.type    = 'SO3' ou 'INDI'
%   sim.scenario       = 'full_mission' | 'transition_only' | 'hover_disturbed'
%   USE_GA_REFINER     = true/false (otimização energética via GA)

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
USE_GA_REFINER = false;        % set true to enable GA energy refinement (slower)
INNER_LAW      = 'SO3';        % 'SO3' or 'INDI'

%% --- Load configurations ---
ac   = aircraft_config();
ctrl = controller_config();
sim  = simulation_config();
ctrl.inner.type = INNER_LAW;

fprintf('=== eVTOL Tailsitter Simulation Framework ===\n');
fprintf('Mass: %.0f kg | Rotors: %d | Inner: %s | Scenario: %s\n', ...
    ac.mass, ac.n_rotors, ctrl.inner.type, sim.scenario);

%% --- Mission and trajectory ---
waypoints = build_mission(sim.scenario);
opt = TrajectoryOptimizer(waypoints, ac);
if USE_GA_REFINER
    fprintf('Optimizing trajectory via GA energy refiner...\n');
    traj = opt.optimize('method','ga','ga_pop',24,'ga_gens',15);
else
    traj = opt.optimize('method','flat');
end
fprintf('Trajectory total time: %.2f s\n', traj.total_time());

%% --- Simulate ---
sim.t_final = min(sim.t_final, traj.total_time() + 5.0);
log = run_simulation(ac, ctrl, sim, traj);

%% --- Visualize ---
plot_telemetry(log, ac);

fprintf('\nDone. Inspect figures for validation.\n');
