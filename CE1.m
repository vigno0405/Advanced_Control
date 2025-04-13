%% 1.1 Norm of SISO systems

clc
close all
clear all
num = [1 -10];
den = [1 5 26.5];
G = tf(num, den);   % second-order model
figure()
bode(G);

% 1.1.1 2-NORM
% 1. Residue theorem (sum of two residuals)
poles = roots(den);
p1 = poles(1);
p2 = poles(2);
Res1 = (p1-10)*(-p1-10)/((p1-p2)*(-p1-p2)*(-2*p1));
Res2 = (p2-10)*(-p2-10)/((p2-p1)*(p2+p1)*(p2*2));
Norm1 = sqrt(Res1+Res2);
disp(['2-norm = ', num2str(Norm1)]);

% 2. Frequency response of G (integral)
omega = linspace(0,5e2,1e5);    % relevant interval
G_jomega = (j*omega-10)./((j*omega).^2+5.*j*omega+26.5);
Norm2 = sqrt(1/pi*trapz(omega, abs(G_jomega).^2));
disp(['2-norm = ', num2str(Norm2)]);

% 3. Impulse response of G
t = linspace(0, 600, 1e5);  % time
g = impulse(G, t);
Norm3 = sqrt(trapz(t, abs(g).^2));
disp(['2-norm = ', num2str(Norm3)]);

% 4. State-space method
[A,B,C,D] = ssdata(G);
L = are(A', zeros(2), B*B');
Norm4 = sqrt(trace(C*L*C'));
disp(['2-norm = ', num2str(Norm4)]);

% 5. Validate results with MATLAB command norm()
disp(['2-norm = ', num2str(norm(G))]);

% 1.1.2. INFINITY-NORM
% 1. Frequency response of G
Norm5 = max(abs(G_jomega));
disp(['Infinity-norm = ', num2str(Norm5)]);

% 2. Bounded real lemma (iterative bisection algorithm)
gamma_u = 1;
gamma_l = 1e-5;     % boundaries
eps = 1e-8;
flag = true;
niter = 0;

while flag
    if (gamma_u-gamma_l)/gamma_l < eps
        flag = false;
        Norm6 = mean([gamma_u, gamma_l]);
    else
    gamma = mean([gamma_u, gamma_l]);
    H = [A, gamma^(-2)*B*B'; -C'*C, -A'];
    if any(abs(real(eig(H)))<1e-3)
        gamma_l = gamma;
    else
        gamma_u = gamma;
    end
    end
    niter = niter + 1;
end
disp(['Infinity-norm = ', num2str(Norm6)]);

% 3. Validate results with norm()
disp(['Infinity-norm = ', num2str(norm(G, inf))]);

%% 1.2 Norms of MIMO systems

clc
close all
clear all
A = [-8.0, 2.5, -2.5; 7.8, -3.2, 13.8; 0.1, -0.5, -5.7];
B = [1, 0; 1, -1; 0, -2];
C = [0, 0, 1; 1, -1, 0];
D = zeros(2);
Gss = ss(A, B, C, D);  % state-space
G = tf(Gss);

% 1.2.1 2-NORM
% 1. Frequency response method
omega = linspace(0,5e4,1e5);    % relevant interval
n = length(A);
N = length(omega);
G_jomega = zeros(size(C,1), size(B,2), N);  % allocation of 3 matrix
trace_GstarG = zeros(1, N);

for k = 1:N
    G_jomega(:,:,k) = C * ((1j*omega(k)*eye(n) - A) \ B) + D;
    trace_GstarG(k) = trace(G_jomega(:,:,k)'*G_jomega(:,:,k));
end

Norm7 = sqrt(1/pi*trapz(omega, trace_GstarG));
disp(['2-norm = ', num2str(Norm7)]);

% 2. State-space method
[A,B,C,D] = ssdata(G);
L = are(A', zeros(size(A)), B*B');
Norm8 = sqrt(trace(C*L*C'));
disp(['2-norm = ', num2str(Norm8)]);

% 3. Validate results with norm()
disp(['2-norm = ', num2str(norm(G, 2))]);

% 1.2.2 INFINITY-NORM
% 1. Frequency response method
sigma = zeros(1, N);
for k = 1:N
    sigma(k) = sqrt(max(eig(G_jomega(:,:,k)'*G_jomega(:,:,k))));
end
Norm9 = max(sigma);
disp(['Infinity-norm = ', num2str(Norm9)]);

% 2. Bounded real lemma (iterative bisection algorithm)
gamma_u = 10;
gamma_l = 1e-5;     % boundaries
eps = 1e-8;
flag = true;
niter = 0;

while flag
    if (gamma_u-gamma_l)/gamma_l < eps
        flag = false;
        Norm10 = mean([gamma_u, gamma_l]);
    else
    gamma = mean([gamma_u, gamma_l]);
    H = [A, gamma^(-2)*B*B'; -C'*C, -A'];
    if any(abs(real(eig(H)))<1e-3)
        gamma_l = gamma;
    else
        gamma_u = gamma;
    end
    end
    niter = niter + 1;
end
disp(['Infinity-norm = ', num2str(Norm10)]);

% 3. Validate results with norm()
disp(['Infinity-norm = ', num2str(norm(G, inf))]);

%% 1.3 Uncertainty modeling

clc
close all
clear all

% Parameters definition
a = ureal('a', 2, 'PlusMinus', 1);
b = ureal('b', 3, 'PlusMinus', 2);
c = ureal('c', 10, 'PlusMinus', 3);

% Define uncertain LTI system and plotting of step, bode, nyquist
num = a;
den = [1 b c];
G = tf(num, den);

% Use usample to generate two multimodel uncertainty sets (20, 200 samples)
sample1 = 20;
sample2 = 200;
G_samples1 = usample(G, sample1);
G_samples2 = usample(G, sample2);

% Creation of G_nominal (according to the nominal values of parameters)
G_nominal = tf(a.NominalValue, [1 b.NominalValue c.NominalValue]);

orders = 1:3; % Array of model orders to test

for order = orders

    % Estimate uncertainty models for the given order
    [usys1, info1] = ucover(G_samples1, G_nominal, order);
    [usys2, info2] = ucover(G_samples2, G_nominal, order);

    % Plot Bode response of uncertainty filters
    figure()
    bode(info1.W1, 'k', info2.W1, 'g');
    legend('W2 with 20 samples', 'W2 with 200 samples')
    title(['Order ', num2str(order)])
    set(gcf, 'Position', get(0, 'Screensize')); % Maximize figure
    saveas(gcf, fullfile('Images_CE1', ['Bode_Uncertainty_Order_', num2str(order), '.jpg']))

    % Plot Bode response of the multiplicative model (20 samples)
    figure()
    bode(usys1, G)
    legend('Multiplicative model (20 samples)', 'G')
    title(['Order ', num2str(order)])
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, fullfile('Images_CE1', ['Bode_Multiplicative_20_Order_', num2str(order), '.jpg']))

    % Plot Bode response of the multiplicative model (200 samples)
    figure()
    bode(usys2, G)
    legend('Multiplicative model (200 samples)', 'G')
    title(['Order ', num2str(order)])
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, fullfile('Images_CE1', ['Bode_Multiplicative_200_Order_', num2str(order), '.jpg']))

    % Plot step response for 20 samples
    figure()
    step(usys1, G)
    legend('Multiplicative model (20 samples)', 'G')
    title(['Order ', num2str(order)])
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, fullfile('Images_CE1', ['Step_20_Order_', num2str(order), '.jpg']))

    % Plot step response for 200 samples
    figure()
    step(usys2, G)
    legend('Multiplicative model (200 samples)', 'G')
    title(['Order ', num2str(order)])
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, fullfile('Images_CE1', ['Step_200_Order_', num2str(order), '.jpg']))

    % Plot Nyquist diagram for 20 samples
    figure()
    nyquist(usys1, G)
    legend('Multiplicative model (20 samples)', 'G')
    title(['Order ', num2str(order)])
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, fullfile('Images_CE1', ['Nyquist_20_Order_', num2str(order), '.jpg']))

    % Plot Nyquist diagram for 200 samples
    figure()
    nyquist(usys2, G)
    legend('Multiplicative model (200 samples)', 'G')
    title(['Order ', num2str(order)])
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, fullfile('Images_CE1', ['Nyquist_200_Order_', num2str(order), '.jpg']))
end
