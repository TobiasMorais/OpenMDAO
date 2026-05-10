function waypoints = build_mission(scenario)
% BUILD_MISSION  Returns waypoint matrix [x y z psi] for canonical scenarios.
%
% NED frame: x=North, y=East, z=Down (negative z = above ground).
% psi = heading (rad). All scenarios start at origin on ground.

    switch lower(scenario)
        case 'full_mission'
            % Vertical takeoff, climb, transition, cruise, decel, vertical land
            waypoints = [
                  0,    0,     0,    0;        % takeoff
                  0,    0,   -10,    0;        % rise to 10 m
                  0,    0,   -50,    0;        % climb to 50 m (still hover-like)
                 60,    0,   -50,    0;        % begin cruise
                400,    0,   -80,    0;        % cruise outbound 1 km
                900,    0,  -100,    0;        % continue cruise
               1200,    0,   -50,    0;        % decelerate descent
               1280,    0,   -10,    0;        % approach
               1280,    0,     0,    0];       % vertical landing

        case 'transition_only'
            waypoints = [
                  0,    0,   -50,   0;
                 50,    0,   -50,   0;
                200,    0,   -60,   0];

        case 'hover_disturbed'
            waypoints = [
                  0,    0,    0,    0;
                  0,    0,   -30,   0;
                  0,    0,   -30,   0];

        case 'hover_only'
            % Pure stationary hover hold (no climb). Used to test cascade
            % stability without trajectory tracking dynamics.
            % Tiny epsilon between waypoints to avoid degenerate polynomial.
            waypoints = [
                  0,    0,   -10,   0;
                  0,    0,   -10.001, 0];

        otherwise
            error('Unknown scenario: %s', scenario);
    end
end
