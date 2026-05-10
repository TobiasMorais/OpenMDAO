function summary = run_all_tests()
% RUN_ALL_TESTS  Master runner for the eVTOL Tailsitter unit test suite.
%
% Executes all 15 test suites and prints aggregated PASS/FAIL summary.
% Each suite is independent; failures in one do not block others.
%
% Usage:
%   cd evtol_simulation
%   tests/run_all_tests
%
% Returns: struct with overall stats.

% Add all needed paths
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, '..', 'config'));
addpath(fullfile(here, '..', 'core'));
addpath(fullfile(here, '..', 'aerodynamics'));
addpath(fullfile(here, '..', 'control'));
addpath(fullfile(here, '..', 'trajectory'));
addpath(fullfile(here, '..', 'environment'));
addpath(fullfile(here, '..', 'simulation'));

print_header();

stages = {
    {'STAGE 1 — Cinemática/Dinâmica', {'test_quat_utils','test_so3_utils','test_aircraft_dynamics'}}, ...
    {'STAGE 2 — Aero-Propulsão',       {'test_aerodynamics','test_propulsion'}}, ...
    {'STAGE 3 — Controle',             {'test_so3_controller','test_indi_controller','test_nmpc','test_allocator','test_force_to_attitude'}}, ...
    {'STAGE 4 — Trajetória',           {'test_flatness','test_trajectory_opt'}}, ...
    {'STAGE 5 — Ambiente',             {'test_rk4','test_dryden','test_atmosphere'}}, ...
    {'STAGE 6 — Integração (cascata)', {'test_integration_hover'}}
};

total_suites  = 0;
total_pass_s  = 0;
total_tests   = 0;
total_pass_t  = 0;
fail_log      = {};

for s = 1:numel(stages)
    fprintf('\n[%s]\n', stages{s}{1});
    suites = stages{s}{2};
    for k = 1:numel(suites)
        suite_name = suites{k};
        total_suites = total_suites + 1;
        try
            tr = feval(suite_name);
            total_tests  = total_tests + tr.ntotal;
            total_pass_t = total_pass_t + tr.npass;
            if tr.passed
                total_pass_s = total_pass_s + 1;
            else
                fail_log{end+1} = sprintf('%s: %d/%d', suite_name, tr.npass, tr.ntotal);
            end
        catch ME
            fprintf('  %-30s ERROR: %s\n', suite_name, ME.message);
            fail_log{end+1} = sprintf('%s: ERROR (%s)', suite_name, ME.message);
        end
    end
end

print_footer(total_suites, total_pass_s, total_tests, total_pass_t, fail_log);

summary.suites_total  = total_suites;
summary.suites_passed = total_pass_s;
summary.tests_total   = total_tests;
summary.tests_passed  = total_pass_t;
summary.fails         = fail_log;
summary.ready_for_integration = (total_pass_s == total_suites);

end

function print_header()
fprintf('\n');
fprintf('===================================================================\n');
fprintf('   eVTOL Tailsitter Simulation Framework — Test Suite Report\n');
fprintf('   Date: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('===================================================================\n');
end

function print_footer(ns, np_s, nt, np_t, fails)
fprintf('\n===================================================================\n');
fprintf('   SUMMARY\n');
fprintf('   Suites: %d/%d passed (%.1f%%)\n', np_s, ns, 100*np_s/max(ns,1));
fprintf('   Tests:  %d/%d passed (%.1f%%)\n', np_t, nt, 100*np_t/max(nt,1));
if isempty(fails)
    fprintf('   STATUS: READY FOR INTEGRATION (run main.m)\n');
else
    fprintf('   STATUS: NOT READY — fix failures below before integration:\n');
    for k = 1:numel(fails)
        fprintf('     - %s\n', fails{k});
    end
end
fprintf('===================================================================\n\n');
end
