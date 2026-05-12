function tr = test_so3_utils()
% TEST_SO3_UTILS  Validation suite for SO(3) Lie-algebra helpers.
tr = test_helpers('init', 'test_so3_utils');
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'core'));

% S1: vee(hat(v)) = v
v = [0.7; -0.2; 0.5];
H = so3_utils('hat', v);
v_back = so3_utils('vee', H);
tr = test_helpers('assert_near', tr, v_back, v, 1e-15, 'S1 vee(hat(v)) = v');

% S2: hat antisymmetric
tr = test_helpers('assert_lt', tr, norm(H + H', 'fro'), 1e-15, 'S2 hat is antisymmetric');

% S3: exp(0) = I
R0 = so3_utils('exp', [0;0;0]);
tr = test_helpers('assert_near', tr, R0, eye(3), 1e-15, 'S3 exp(0) = I');

% S4: exp(log(R)) = R for random R
rng(11);
for k = 1:5
    phi = randn(3,1) * 0.7;
    R = so3_utils('exp', phi);
    R_back = so3_utils('exp', so3_utils('log', R));
    err = norm(R - R_back, 'fro');
end
tr = test_helpers('assert_lt', tr, err, 1e-9, 'S4 exp(log(R)) = R');

% S5: errMat(R, R) = 0
R = so3_utils('exp', [0.3; 0.1; -0.2]);
e = so3_utils('errMat', R, R);
tr = test_helpers('assert_lt', tr, norm(e), 1e-12, 'S5 errMat(R,R) = 0');

% S6: small angle errMat ~ small rotation
small_phi = [0.0; 0.005; 0.0];
Rd = so3_utils('exp', small_phi);
R0 = eye(3);
e_R = so3_utils('errMat', R0, Rd);
expected = -small_phi;   % errMat(R, Rd) gives error of R wrt Rd (negative for forward)
err_rel = norm(e_R - expected) / norm(expected);
tr = test_helpers('assert_lt', tr, err_rel, 1e-2, 'S6 small angle errMat linearization');

% S7: det(exp(phi)) = 1
phi_t = [0.4; -0.3; 0.6];
R_t = so3_utils('exp', phi_t);
tr = test_helpers('assert_lt', tr, abs(det(R_t)-1), 1e-12, 'S7 det(exp(phi)) = 1');

% S8: Rodrigues consistent with axis-angle
ax = [0;0;1];
ang = pi/3;
Rrod = so3_utils('exp', ang*ax);
Rexpected = [cos(ang), -sin(ang), 0; sin(ang), cos(ang), 0; 0, 0, 1];
tr = test_helpers('assert_near', tr, Rrod, Rexpected, 1e-12, 'S8 Rodrigues axis-angle');

tr = test_helpers('report', tr);
end
