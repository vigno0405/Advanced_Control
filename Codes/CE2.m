%% 2.1.1. Gf and G creation

clear
close all
clc
load('data_position.mat')

% Number of experiments and structure initialization
N_exp = 4;
Gf_struct(N_exp) = struct('Gf', []);
G_struct(N_exp) = struct('G', []);
G_der_struct(N_exp) = struct('G_der', []);

% Frequency
freqs = (pi/4096:pi/4096:pi)/Ts;

for ii=1:N_exp

    % Data extraction
    experiment_name = data{ii}.name;
    y = data{ii}.y;
    u = data{ii}.u;
    Z = detrend(iddata(y, u, Ts, "Period", 8191));

    % Frequency response of the system
    Gf_struct(ii).Gf = spa(Z, 8191, freqs);

    % Compute output derivative
    y_derivative = lsim(1 - tf('z',Ts)^-1, y);
    Z_der = detrend(iddata(y_derivative, u, Ts, "Period", 8191));

    % Use Output Error method to compute the model
    z = tf('z', Ts);
    G_der_struct(ii).G_der = oe(Z_der, [10, 10, 1]);    % order (x2), shift

    % Add back the integrator
    G_struct(ii).G = G_der_struct(ii).G_der / (1 - z^-1);

end

save('Model.mat')

%% 2.1.2. Non-parametric models

clear
close all
load('Model.mat')

% Nyquist options
opts = nyquistoptions;
opts.ConfidenceRegionDisplaySpacing = 3;
opts.ShowFullContour = 'off';

for ii=1:N_exp
    figure()
    nyquistplot(Gf_struct(ii).Gf, freqs, opts, 'sd', 2.45)
    hold on
    axis equal
    title(['Non-par, exp.', num2str(ii)])
end

% Nonparametric
figure()
bodemag(Gf_struct(1).Gf)
hold on
bodemag(Gf_struct(2).Gf)
bodemag(Gf_struct(3).Gf)
bodemag(Gf_struct(4).Gf)
legend('1', '2', '3', '4')
title('Nonparametric')

%% 2.1.3. Parametric models

clear
close all
load('Model.mat')

% Nyquist options
opts = nyquistoptions;
opts.ConfidenceRegionDisplaySpacing = 3;
opts.ShowFullContour = 'off';

for ii=1:N_exp
    figure()
    nyquistplot(G_der_struct(ii).G_der, freqs, opts, 'sd', 2.45)
    axis equal
    hold on
    title(['Param, exp.', num2str(ii)])
end

% Parametric
figure()
bodemag(G_struct(1).G)
hold on
bodemag(G_struct(2).G)
bodemag(G_struct(3).G)
bodemag(G_struct(4).G)
legend('1', '2', '3', '4')
title('Parametric')

%% 2.1.4. Multimodel set from parametric models

clear
close all
load('Model.mat')

Gmm = stack(1, G_struct(1).G, G_struct(2).G, G_struct(3).G, ...
    G_struct(4).G);

%% 2.1.5. W2 filter

clear
close all
load('Model.mat')

Gmm = stack(1, G_struct(1).G, G_struct(2).G, G_struct(3).G, ...
    G_struct(4).G);
Gmm_der = stack(1, G_der_struct(1).G_der, G_der_struct(2).G_der,...
    G_der_struct(3).G_der, G_der_struct(4).G_der);

% Order 20 of filter
order = 20;

% Computation for every one
info = struct('info', [], 'W2', []);
for ii=1:4
    Gnom = G_der_struct(ii).G_der;
    [Gu, info(ii).info] = ucover(Gmm_der, Gnom, order, 'InputMult');
    info(ii).W2 = info(ii).info.W1opt;
end

% Plotting on bodemag
figure()
for ii = 1:4
    bodemag(info(ii).W2)
    hold on
end
legend('1', '2', '3', '4')
title('Uncertainty filters for the original system')

% We choose the second one because at low frequencies they are more or less
% the same (1 dB), while 1 and 3 are too big at high frequencies and the
% fourth one has too big bandwidth in high frequencies.

% Plot the filter
figure()
bodemag(info(2).info.W1, info(2).info.W1opt)
legend('W2', 'W2opt')

% Saving information
W2 = info(2).info.W1;
G_nom = G_struct(2).G;
save('W2.mat', 'W2')
save('G_functions.mat', 'G_nom', 'G_struct', 'Gf_struct')

%% 2.1.6. Comparison with original

clc
close all
clear

load('W2.mat')
load('G_functions.mat')

Gmm = stack(1, G_struct(1).G, G_struct(2).G, ...
    G_struct(3).G, G_struct(4).G);

figure()
hold on
Gf = frd(Gmm, logspace(-1,3,60));
bodemag((G_nom - Gf)/G_nom, 'b--', W2, 'r')
title('Relative gaps vs magnitude of W2')
grid on

%% 2.2.1. W1

% 1. Zero steady-state tracking error:
% - integrator needed, for numerical stability use a quasi-integrator
% - modulus margin of at least 0.5, so since it is the distance between the
% Nyquist plot and (-1, 0), and the filter is increasing in magnitude, at
% high frequencies the filter should have a gain lower than 0.5
% - minimum settling time for the nominal model

clear
close all
clc
load('Model.mat')

num = [1 4];        % settling time
den = [1 1e-5];     % quasi-integrator
W1 = tf(num, den)*1/2;
W1 = c2d(W1, Ts, 'zoh');

disp(W1)

% Checking behavior in high frequencies
figure()
bodemag(W1^-1, tf(2))
legend('W1^{-1}', '6dB')

save('W1.mat', 'W1')

%% 2.2.2. Mixed sensitivity approach

clear
close all
clc

load('Model.mat')
load('W2.mat')
load('W1.mat')
load('G_functions.mat')

Gmm = stack(1, G_struct(1).G, G_struct(2).G, ...
    G_struct(3).G, G_struct(4).G);
Kinf = mixsyn(G_nom, W1, [], W2);

%% 2.2.3. Validation plots and condition for robust performance

clear
close all
clc

load('Model.mat')
load('W2.mat')
load('W1.mat')
load('G_functions.mat')

Gmm = stack(1, G_struct(1).G, G_struct(2).G, ...
    G_struct(3).G, G_struct(4).G);
Kinf = mixsyn(G_nom, W1, [], W2);

Tnom = feedback(G_nom*Kinf, 1);
Snom = feedback(1, G_nom*Kinf);
Unom = feedback(Kinf, G_nom);
T = feedback(Gmm*Kinf, 1);
S = feedback(1, Gmm*Kinf);
U = feedback(Kinf, Gmm);

% Plot step response, control signal
figure()
subplot(1,2,1)
step(T,Tnom)
title('Step response')
legend('Multimodel','Nominal')
subplot(1,2,2)
step(U,Unom)
title('Control signal')
legend('Multimodel','Nominal')

figure()
subplot(1,3,1)
bodemag(U,Unom)
title('Sensitivity function U')
legend('Multimodel','Nominal')
subplot(1,3,2)
bodemag(S,Snom, W1^-1)
title('Sensitivity function S')
legend('Multimodel','Nominal', 'W1 inv')
subplot(1,3,3)
bodemag(T,Tnom,W2^-1)
title('Sensitivity function T')
legend('Multimodel','Nominal','W2 inv')

nominal_performance = norm(W1*S, inf);
disp('Nominal performance')
disp(nominal_performance)

robust_stability = norm(W2*Tnom, inf);
disp('Robust stability')
disp(robust_stability)

robust_performance = norm([W1*Snom, W2*Tnom], inf);
disp('Robust performance')
disp(robust_performance)

%% 2.2.4. W3

clear
clc
close all

load('Model.mat')
load('W2.mat')
load('W1.mat')
load('G_functions.mat')

Gmm = stack(1, G_struct(1).G, G_struct(2).G, ...
    G_struct(3).G, G_struct(4).G);

W3 = tf(0.05);     % trial and error
save('W3.mat', 'W3')

Kinf_3 = mixsyn(G_nom, W1, W3, W2);

Tnom = feedback(G_nom*Kinf_3, 1);
Snom = feedback(1, G_nom*Kinf_3);
Unom = feedback(Kinf_3, G_nom);
T = feedback(Gmm*Kinf_3, 1);
S = feedback(1, Gmm*Kinf_3);
U = feedback(Kinf_3, Gmm);

figure()
subplot(1,2,1)
step(T,Tnom)
title('Step response')
legend('Multimodel','Nominal')
subplot(1,2,2)
step(U,Unom)
title('Control signal')
legend('Multimodel','Nominal')

figure()
subplot(1,3,1)
bodemag(U,Unom,W3^-1)
title('Sensitivity function U')
legend('Multimodel','Nominal')
subplot(1,3,2)
bodemag(S,Snom, W1^-1)
title('Sensitivity function S')
legend('Multimodel','Nominal', 'W1 inv')
subplot(1,3,3)
bodemag(T,Tnom,W2^-1)
title('Sensitivity function T')
legend('Multimodel','Nominal','W2 inv')

nominal_performance = norm(W1*S, inf);
disp('Nominal performance')
disp(nominal_performance)

robust_stability = norm(W2*Tnom, inf);
disp('Robust stability')
disp(robust_stability)

robust_performance = norm([W1*Snom, W2*Tnom], inf);
disp('Robust performance')
disp(robust_performance)

save('H_inf_controller.mat')

%% 2.2.5. The order of the controller may be too large: pzmap.

% Let's consider the last controller on position
clear
clc
close all
load('H_inf_controller.mat')

Gmm = stack(1, G_struct(1).G, G_struct(2).G, ...
    G_struct(3).G, G_struct(4).G);
Kinf = Kinf_3;

% Order of Kinf
orderK = size(Kinf.A,1);
disp(['Order of the controller: ', num2str(orderK)])

% High order of the controller, especially due to the high order of W2
figure()
pzmap(Kinf)
title('pzmap - before')

% Lots of poles and zeros are really close, this means that there is
% zero/pole cancellation in the controller.

% Reduction
Kreduced = reduce(Kinf, 25);
orderK = size(Kreduced.A, 1);
disp(['Order of the controller: ', num2str(orderK)])
figure()
pzmap(Kreduced)
title('pzmap - reduced')

Tnom = feedback(G_nom*Kreduced, 1);
Snom = feedback(1, G_nom*Kreduced);
Unom = feedback(Kreduced, G_nom);
T = feedback(Gmm*Kreduced, 1);
S = feedback(1, Gmm*Kreduced);
U = feedback(Kreduced, Gmm);

figure()
subplot(1,2,1)
step(T,Tnom)
title('Step response')
legend('Multimodel','Nominal')
subplot(1,2,2)
step(U,Unom)
title('Control signal')
legend('Multimodel','Nominal')

figure()
subplot(1,3,1)
bodemag(U,Unom,W3^-1)
title('Sensitivity function U')
legend('Multimodel','Nominal')
subplot(1,3,2)
bodemag(S,Snom, W1^-1)
title('Sensitivity function S')
legend('Multimodel','Nominal', 'W1 inv')
subplot(1,3,3)
bodemag(T,Tnom,W2^-1)
title('Sensitivity function T')
legend('Multimodel','Nominal','W2 inv')

% Plot of the original controller and of the reduced one
figure()
bodemag(Kinf, Kreduced)
legend('Kinf', 'Kreduced')

nominal_performance = norm(W1*S, inf);
disp('Nominal performance')
disp(nominal_performance)

robust_stability = norm(W2*Tnom, inf);
disp('Robust stability')
disp(robust_stability)

robust_performance = norm([W1*Snom, W2*Tnom], inf);
disp('Robust performance')
disp(robust_performance)

save('Norm_validation.mat')

%% 2.3.1. Conversion to continuous time model

clc, clear, close all
load("Model.mat")

Gct = d2c(G_struct(2).G);
[A, B, C, D] = ssdata(Gct);
save('continuous_model.mat')

%% 3.3.4. H2 controller

clear, clc, close all
load('continuous_model.mat')
n = size(A, 1);
m = size(B, 2);

% Decision variables
L = sdpvar(n, n, 'symmetric');
X = sdpvar(m, n);
M = sdpvar(m, m, 'symmetric');

% Objective
obj = trace(C*L*C') + trace(M);

% LMIs
lmi1 = A*L + L*A' - B*X - X'*B' + B*B' <= 0;
lmi2 = [M X; X' L] >= 0;
lmi3 = L >= eps;
lmi = [lmi1, lmi2, lmi3];

% MOSEK optimization
options = sdpsettings('solver', 'mosek');
optimize(lmi, obj, options);

% H2 controller
K_H2 = value(X) * inv(value(L));
save('H_2_controller.mat', 'A', 'B', 'C', 'D', 'K_H2', 'm', 'n')

%% 3.3.5. Step response of the closed-loop system

clear, clc, close all
load('H_2_controller.mat')

% State-space representation
Acl = A - B*K_H2;
Bcl = B;
Ccl = C;
Dcl = D;
sys_closedloop = ss(Acl,Bcl,Ccl,Dcl);

figure()
step(sys_closedloop)
title('Step response using H2-controller')
save('H_2_system.mat')

%% 3.3.6. Comparison with LQR

clear, clc, close all
load('H_2_system.mat')

% LQR
Q = C'*C;           % Q (min. errors on state)
R = eye(m);         % R (min. errors on input)
[K_lqr, ~, ~] = lqr(A, B, Q, R);
Acl_lqr = A - B*K_lqr;
sys_closedloop_lqr = ss(Acl_lqr, Bcl, Ccl, Dcl);

figure()
step(sys_closedloop)
hold on
step(sys_closedloop_lqr)
title('Comparison between LQR and H2 controller')
legend('H-2', 'LQR')

%% 4.4.1. Data-driven controller: multiplicative uncertainty

clear, clc, close all
rng(22)
load("Model.mat")
load("G_functions.mat")
load('W1.mat')
load('W2.mat')
load('W3.mat')

Gmm = stack(1, G_struct(1).G, G_struct(2).G, ...
    G_struct(3).G, G_struct(4).G);
omegas = unique([ ...
    logspace(log10(0.4), log10(pi/Ts), 200), ...
    linspace(300, 400, 101), ...
    linspace(100, 1000, 101), ...
    ]);

% 1. Proportional initial controller
order = 5;
ny = 1;
nu = 1;
K_c = c2d(tf(1), Ts, 'zoh');

% Check that K_c is a stabilizing controller
CL_nom = feedback(1, K_c*G_nom);
isstable(CL_nom)

% LFR model
P = augw(G_nom, W1, W3, W2);
P = mktito(P, ny, nu);

% Problem formulation
K = datadriven.Controller.SS(order, ny, nu, Ts);
K.setinitial(K_c);

% Objective and ensuring stability (not needed in bounded real lemma)
synthesiser = Synthesiser(K);
synthesiser.add_Hinf_objective(P, omegas);
synthesiser.ensure_controller_stability(omegas);

% Output
output = synthesiser.synthesise();
K_final = output.Controller;

% Plotting
Tnom = feedback(G_nom*K_final, 1);
Snom = feedback(1, G_nom*K_final);
Unom = feedback(K_final, G_nom);
T = feedback(Gmm*K_final, 1);
S = feedback(1, Gmm*K_final);
U = feedback(K_final, Gmm);

figure()
subplot(1,2,1)
step(T,Tnom)
title('Step response')
legend('Multimodel','Nominal')
subplot(1,2,2)
step(U,Unom)
title('Control signal')
legend('Multimodel','Nominal')

figure()
subplot(1,3,1)
bodemag(U,Unom,W3^-1)
title('Sensitivity function U')
legend('Multimodel','Nominal')
subplot(1,3,2)
bodemag(S,Snom, W1^-1)
title('Sensitivity function S')
legend('Multimodel','Nominal', 'W1 inv')
subplot(1,3,3)
bodemag(T,Tnom,W2^-1)
title('Sensitivity function T')
legend('Multimodel','Nominal','W2 inv')

% 2. Comparison with model-based controller
load('H_inf_controller.mat')
Kinf = Kinf_3;

% Order of Kinf
orderK = size(Kinf.A, 1);
disp(['Order of the controller: ', num2str(orderK)])

% Reduction to order 5
Kreduced = reduce(Kinf, 5);
orderK = size(Kreduced.A, 1);
disp(['Order of the controller: ', num2str(orderK)])

Tnom = feedback(G_nom*Kreduced, 1);
Snom = feedback(1, G_nom*Kreduced);
Unom = feedback(Kreduced, G_nom);
T = feedback(Gmm*Kreduced, 1);
S = feedback(1, Gmm*Kreduced);
U = feedback(Kreduced, Gmm);

figure()
subplot(1,2,1)
step(T,Tnom)
title('Step response')
legend('Multimodel','Nominal')
subplot(1,2,2)
step(U,Unom)
title('Control signal')
legend('Multimodel','Nominal')

figure()
subplot(1,3,1)
bodemag(U,Unom,W3^-1)
title('Sensitivity function U')
legend('Multimodel','Nominal')
subplot(1,3,2)
bodemag(S,Snom, W1^-1)
title('Sensitivity function S')
legend('Multimodel','Nominal', 'W1 inv')
subplot(1,3,3)
bodemag(T,Tnom,W2^-1)
title('Sensitivity function T')
legend('Multimodel','Nominal','W2 inv')

% 3. Robust stability constraint
P_constraint = augw(G_nom, [], [], W2);
P_constraint = mktito(P_constraint, ny, nu);
P = augw(G_nom, W1, W3, []);
P = mktito(P, ny, nu);
K = datadriven.Controller.SS(order, ny, nu, Ts);
K.setinitial(K_c);
synthesiser = Synthesiser(K);
synthesiser.add_Hinf_objective(P, omegas);
synthesiser.add_Hinf_constraint(P_constraint, omegas);
synthesiser.ensure_controller_stability(omegas);

% Output
output = synthesiser.synthesise();
K_final_constrained = output.Controller;

% Plotting
Tnom = feedback(G_nom*K_final_constrained, 1);
Snom = feedback(1, G_nom*K_final_constrained);
Unom = feedback(K_final_constrained, G_nom);
T = feedback(Gmm*K_final_constrained, 1);
S = feedback(1, Gmm*K_final_constrained);
U = feedback(K_final_constrained, Gmm);

figure()
subplot(1,2,1)
step(T,Tnom)
title('Step response')
legend('Multimodel','Nominal')
subplot(1,2,2)
step(U,Unom)
title('Control signal')
legend('Multimodel','Nominal')

figure()
subplot(1,3,1)
bodemag(U,Unom,W3^-1)
title('Sensitivity function U')
legend('Multimodel','Nominal')
subplot(1,3,2)
bodemag(S,Snom, W1^-1)
title('Sensitivity function S')
legend('Multimodel','Nominal', 'W1 inv')
subplot(1,3,3)
bodemag(T,Tnom,W2^-1)
title('Sensitivity function T')
legend('Multimodel','Nominal','W2 inv')

nominal_performance = norm(W1*S, inf);
disp('Nominal performance')
disp(nominal_performance)

robust_stability = norm(W2*Tnom, inf);
disp('Robust stability')
disp(robust_stability)

%% 4.4.2 Data-driven controller: multimodel uncertainty

clear, clc, close all
rng(22)
load("Model.mat")
load("G_functions.mat")
load('W3.mat')
load('W2.mat')
load('W1.mat')

omegas = unique([ ...
    logspace(log10(0.4), log10(pi/Ts), 200), ...
    linspace(300, 400, 101), ...
    linspace(100, 1000, 101), ...
    ]);

% 1. Proportional initial controller
order = 5;
ny = 1;
nu = 1;
K_c = c2d(tf(1), Ts, 'zoh');

% Check that K_c is a stabilizing controller
Gmm = stack(1, G_struct(1).G, G_struct(2).G, ...
    G_struct(3).G, G_struct(4).G);
isstable(feedback(1, Gmm))

% LFR model
P = augw(Gmm, W1, W3, []);
P = mktito(P, ny, nu);

% Problem formulation
K = datadriven.Controller.SS(order, ny, nu, Ts);
K.setinitial(K_c);

% Objective and ensuring stability (not needed in bounded real lemma)
synthesiser = Synthesiser(K);
synthesiser.add_Hinf_objective(P, omegas);
synthesiser.ensure_controller_stability(omegas);

% Output
output = synthesiser.synthesise();
K_final = output.Controller;

% Plotting
T = feedback(Gmm*K_final, 1);
S = feedback(1, Gmm*K_final);
U = feedback(K_final, Gmm);

figure()
subplot(1,2,1)
step(T)
title('Step response')
legend('Multimodel')
subplot(1,2,2)
step(U)
title('Control signal')
legend('Multimodel')

figure()
subplot(1,2,1)
bodemag(U,W3^-1)
title('Sensitivity function U')
legend('Multimodel','W3 inv')
subplot(1,2,2)
bodemag(S,W1^-1)
title('Sensitivity function S')
legend('Multimodel','W1 inv')

% 2. Improved settling time (comparison)
num = [1, 20];
den = [1, 1e-5];
W1 = tf(num, den)*1/2;
W1 = c2d(W1, Ts, 'zoh');

P = augw(G_nom, W1, W3, W2);
P = mktito(P, ny, nu);

% Problem formulation
K = datadriven.Controller.SS(order, ny, nu, Ts);
K.setinitial(K_c);

% Objective and ensuring stability (not needed in bounded real lemma)
synthesiser = Synthesiser(K);
synthesiser.add_Hinf_objective(P, omegas);
synthesiser.ensure_controller_stability(omegas);

% Output
output = synthesiser.synthesise();
K_final = output.Controller;

% Plotting
Tnom = feedback(G_nom*K_final, 1);
Snom = feedback(1, G_nom*K_final);
Unom = feedback(K_final, G_nom);
T = feedback(Gmm*K_final, 1);
S = feedback(1, Gmm*K_final);
U = feedback(K_final, Gmm);

figure()
subplot(1,2,1)
step(T,Tnom)
title('Step response')
legend('Multimodel','Nominal')
subplot(1,2,2)
step(U,Unom)
title('Control signal')
legend('Multimodel','Nominal')

figure()
subplot(1,3,1)
bodemag(U,Unom,W3^-1)
title('Sensitivity function U')
legend('Multimodel','Nominal')
subplot(1,3,2)
bodemag(S,Snom, W1^-1)
title('Sensitivity function S')
legend('Multimodel','Nominal', 'W1 inv')
subplot(1,3,3)
bodemag(T,Tnom,W2^-1)
title('Sensitivity function T')
legend('Multimodel','Nominal','W2 inv')

% Multimodel
P = augw(Gmm, W1, W3, []);
P = mktito(P, ny, nu);

% Problem formulation
K = datadriven.Controller.SS(order, ny, nu, Ts);
K.setinitial(K_c);

% Objective and ensuring stability (not needed in bounded real lemma)
synthesiser = Synthesiser(K);
synthesiser.add_Hinf_objective(P, omegas);
synthesiser.ensure_controller_stability(omegas);

% Output
output = synthesiser.synthesise();
K_final = output.Controller;

% Plotting
T = feedback(Gmm*K_final, 1);
S = feedback(1, Gmm*K_final);
U = feedback(K_final, Gmm);

figure()
subplot(1,2,1)
step(T)
title('Step response')
legend('Multimodel')
subplot(1,2,2)
step(U)
title('Control signal')
legend('Multimodel')

figure()
subplot(1,2,1)
bodemag(U,W3^-1)
title('Sensitivity function U')
legend('Multimodel','W3 inv')
subplot(1,2,2)
bodemag(S,W1^-1)
title('Sensitivity function S')
legend('Multimodel','W1 inv')