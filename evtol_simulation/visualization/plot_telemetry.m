function plot_telemetry(log, ac_cfg)
% PLOT_TELEMETRY  Aerospace-grade validation plots from a simulation log.
%
% Plots created:
%   1) 3D trajectory NED  (and projections)
%   2) Position vs time (NED)
%   3) Velocity NED vs time
%   4) Body angular rates p,q,r
%   5) Pitch angle theta(t) showing the 90-deg transition
%   6) Per-rotor thrust (actual vs commanded)
%   7) Differential motor effort (rotor 1-4 spread)
%   8) Tracking error (||e_pos||)

    t = log.t;

    figure('Name','3D Trajectory','Position',[100 100 800 600]);
    plot3(log.pos_NED(2,:), log.pos_NED(1,:), -log.pos_NED(3,:), 'b-', 'LineWidth', 1.5);
    hold on; grid on;
    plot3(log.pos_NED(2,1), log.pos_NED(1,1), -log.pos_NED(3,1), 'go', 'MarkerSize', 10, 'LineWidth', 2);
    plot3(log.pos_NED(2,end), log.pos_NED(1,end), -log.pos_NED(3,end), 'rs', 'MarkerSize', 10, 'LineWidth', 2);
    xlabel('East [m]'); ylabel('North [m]'); zlabel('Altitude [m]');
    title('3D Trajectory  (start=green, end=red)');
    axis equal; view(45, 25);

    figure('Name','Position NED','Position',[100 100 800 600]);
    subplot(3,1,1); plot(t, log.pos_NED(1,:),'b'); grid on; ylabel('North [m]');
    title('Position (NED frame)');
    subplot(3,1,2); plot(t, log.pos_NED(2,:),'b'); grid on; ylabel('East [m]');
    subplot(3,1,3); plot(t, -log.pos_NED(3,:),'b'); grid on; ylabel('Altitude [m]'); xlabel('t [s]');

    figure('Name','Velocity NED','Position',[100 100 800 600]);
    subplot(3,1,1); plot(t, log.vel_NED(1,:),'b'); grid on; ylabel('V_N [m/s]');
    subplot(3,1,2); plot(t, log.vel_NED(2,:),'b'); grid on; ylabel('V_E [m/s]');
    subplot(3,1,3); plot(t, log.vel_NED(3,:),'b'); grid on; ylabel('V_D [m/s]'); xlabel('t [s]');

    figure('Name','Body Rates p,q,r','Position',[100 100 800 600]);
    plot(t, rad2deg(log.omega_body(1,:)),'r', t, rad2deg(log.omega_body(2,:)),'g', ...
         t, rad2deg(log.omega_body(3,:)),'b');
    grid on; legend('p','q','r'); ylabel('Body rates [deg/s]'); xlabel('t [s]');
    title('Angular rates in body frame');

    figure('Name','Pitch transition','Position',[100 100 800 600]);
    plot(t, log.pitch_deg, 'b-', 'LineWidth', 1.5); grid on;
    yline(90, 'k--', 'Hover (\theta=90\circ)');
    yline(0,  'k--', 'Cruise (\theta=0\circ)');
    xlabel('t [s]'); ylabel('\theta [deg]');
    title('Pitch angle vs time  (90\circ -> 0\circ -> 90\circ for full mission)');

    figure('Name','Per-rotor thrust','Position',[100 100 800 600]);
    plot(t, log.thrust_actual(1,:),'r', ...
         t, log.thrust_actual(2,:),'g', ...
         t, log.thrust_actual(3,:),'b', ...
         t, log.thrust_actual(4,:),'m'); grid on;
    ylabel('Thrust [N]'); xlabel('t [s]');
    legend('Rotor 1 (R, upper)','Rotor 2 (R, lower)','Rotor 3 (L, upper)','Rotor 4 (L, lower)');
    title('Per-rotor actual thrust');

    figure('Name','Differential motor effort','Position',[100 100 800 600]);
    avg = mean(log.thrust_actual, 1);
    diff_eff = log.thrust_actual - avg;
    plot(t, diff_eff', 'LineWidth', 1.0); grid on;
    ylabel('\DeltaT_i [N]'); xlabel('t [s]');
    title('Differential thrust per rotor (T_i - mean)');
    legend('R1','R2','R3','R4');

    figure('Name','Tracking error','Position',[100 100 800 600]);
    err = vecnorm(log.tracking_err, 2, 1);
    plot(t, err, 'b-', 'LineWidth', 1.5); grid on;
    ylabel('||p - p_{ref}|| [m]'); xlabel('t [s]');
    title('Position tracking error');

    fprintf('Maximum tracking error: %.2f m\n', max(err));
    fprintf('Total energy proxy J = int T^1.5 dt = %.3e\n', ...
        trapz(t, log.thrust_total .^ 1.5));
end
