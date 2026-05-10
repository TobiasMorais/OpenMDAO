function tr = test_quat_utils()
% TEST_QUAT_UTILS  Validation suite for Hamilton quaternion utilities.
% See docs/TEST_PROCEDURE.md (Stage 1a) for criteria.

tr = test_helpers('init', 'test_quat_utils');
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'core'));

q_id = [1; 0; 0; 0];

% Q1: Hamilton identity
q  = quat_utils('norm', [0.5; 0.3; -0.7; 0.2]);
qm = quat_utils('mul', q, q_id);
tr = test_helpers('assert_near', tr, qm, q, 1e-12, 'Q1 Hamilton identity multiplication');

% Q2: Conjugate inverse
qc = quat_utils('conj', q);
qe = quat_utils('mul', q, qc);
tr = test_helpers('assert_near', tr, qe, q_id, 1e-12, 'Q2 q (x) q* = identity');

% Q3: R(I) = I3
R_id = quat_utils('toR', q_id);
tr = test_helpers('assert_near', tr, R_id, eye(3), 1e-12, 'Q3 R(identity) = I_3');

% Q4: R orthogonal
R = quat_utils('toR', q);
tr = test_helpers('assert_near', tr, R'*R, eye(3), 1e-10, 'Q4 R^T R = I');

% Q5: Euler roundtrip (random sample)
rng(42);
errs = zeros(20,1);
for k = 1:20
    phi = -pi + 2*pi*rand();
    th  = -pi/2 + 0.95*pi*rand() - 0.45*pi;   % avoid gimbal lock |th|>=85deg
    psi = -pi + 2*pi*rand();
    qq  = quat_utils('fromEuler', phi, th, psi);
    eu  = quat_utils('toEuler', qq);
    errs(k) = norm([phi;th;psi] - eu);
end
tr = test_helpers('assert_lt', tr, max(errs), 1e-9, 'Q5 Euler->quat->Euler roundtrip');

% Q6: kinematic with omega=0 -> qdot=0
qd0 = quat_utils('kinematic', q, [0;0;0]);
tr = test_helpers('assert_lt', tr, norm(qd0), 1e-15, 'Q6 zero omega gives zero qdot');

% Q7: kinematic for omega=[0;0;1], q=I -> qdot = 0.5*[0;0;0;1]
qd1 = quat_utils('kinematic', q_id, [0;0;1]);
tr = test_helpers('assert_near', tr, qd1, [0;0;0;0.5], 1e-12, 'Q7 kinematic axis-z');

% Q8: ExpMap(0) = identity
qe0 = quat_utils('expMap', [0;0;0]);
tr = test_helpers('assert_near', tr, qe0, q_id, 1e-15, 'Q8 ExpMap(0) = identity');

% Q9: LogMap(ExpMap(phi)) = phi
phi_test = [0.3; -0.5; 0.2];
qe1 = quat_utils('expMap', phi_test);
phi_back = quat_utils('logMap', qe1);
tr = test_helpers('assert_near', tr, phi_back, phi_test, 1e-9, 'Q9 LogMap(ExpMap)');

% Q10: rotation about z by 90deg: [1;0;0] -> [0;1;0]
qz90 = quat_utils('fromEuler', 0, 0, pi/2);
v_rot = quat_utils('rotate', qz90, [1;0;0]);
tr = test_helpers('assert_near', tr, v_rot, [0;1;0], 1e-12, 'Q10 Rz(90) [1;0;0] = [0;1;0]');

% Q11: errMul with same quaternion gives identity
qe_self = quat_utils('errMul', q, q);
tr = test_helpers('assert_near', tr, qe_self, q_id, 1e-12, 'Q11 errMul(q,q) = identity');

% Q12: norm preservation through 100 random multiplications
qaccum = q_id;
rng(7);
for k = 1:100
    qrand = quat_utils('norm', randn(4,1));
    qaccum = quat_utils('norm', quat_utils('mul', qaccum, qrand));
end
tr = test_helpers('assert_lt', tr, abs(norm(qaccum)-1), 1e-9, 'Q12 norm preservation 100 muls');

tr = test_helpers('report', tr);
end
