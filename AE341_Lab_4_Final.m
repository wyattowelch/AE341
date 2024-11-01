% ------------------------------------------
% 
% AE341 Lab 4 - Wyatt Welch
%
% ------------------------------------------

%clc, clear all

% Data Gathered

v1 = [.1 .15 .15 .15 .25 .25 .4 .4 .4 .6 .6 .8 .8] / 1e3; % m^3
v2 = [.15 .25 .25 .25 .4 .4 .6 .6 .6 .9 .9 1  1] / 1e3; % m^3
t1 = [33.7 42.4 37.2 32.6 46.2 37.9 49.7 40.5 ...
        35.7 42.1 33.8 35.9 34.4]; % s
t2 = [51.0 70.2 62.9 54.0 74.2 61.3 74.7 61.7 ...
        53.3 63.3 51.3 45.3 42.9]; % s
hl = [78.1 97.7 122.1 152.6 190.7 286.1 429.2   ...
        642.3 963.4 1445.2 2176.7 3251.6 3652.3] / 1e3; %m

mu = 1e-3;
d = .003;
A = 7.06858e-6; 
rho = 998;
Length = .524;

% Calculate Re, f

Q1 = v1 ./ t1; % m^3/s
Q2 = v2 ./ t2; % m^3/s
Qa = .5 * (Q1 + Q2); % m^3/s

v = Qa ./ A; % m/s

Re_d = (rho * v * d) / mu;
f_d = 64 ./ Re_d;

% Head Loss v Velocity Log-Log Graph

xsum = -9.4723;
ysum = 1.9694;
xmean = -.7286;
ymean = .1515;
SS = 21.9992;
SP = 11.3888;
bfit = SP/SS;
afit = ymean - bfit * xmean;
yfit = bfit * log(hl) + afit;

figure()
scatter(log(hl), log(v))
hold on, grid on
plot (log(hl),yfit)
xlabel('Log(hL) [m]')
ylabel('Log(v) [m/s]')
title('Log(hL) v. Log(v)')
legend('Data', 'Best Fit', 'Location', 'northwest')



% PLOTTING MOODY DIAGRAM

figure()
x = [-1.21 -1.23 -1.17 -1.12 -1.05 -1.05 -0.95 -0.87  ...
        -0.77 -0.68 -0.61 -0.47 -0.27];
y = [-0.493 -0.512 -0.451 -0.406 -0.326 -0.326 -0.229 ...
        -0.150 -0.053 0.037 0.109 0.248 0.452];
N = length(x);

% Correlation Coefficient
correlation = (N * sum(x .* y) - sum(x) * sum(y)) / ...
    (sqrt(N * sum(x .^ 2) - sum(x) ^ 2) * ...
    sqrt(N * sum(y .^ 2) - sum(y) ^ 2));
disp(correlation);

% Relative Wall Roughness
e_D = [0 2e-4 5e-4 0.001 0.002 0.005 0.01 0.015 0.02 ...
        0.03 0.04 0.05];
L = length(e_D);

% Re Range
a = 2;
k = 3;
Re = [];
while k < 8.5
    Re = [Re, a * 10^k];
    a = a + 1;
    if a == 10
        a = 1;
        k = k + 1;
    end
end
N = length(Re);
f = zeros(1, N); 

% Experimental
Re_data = Re_d;
f_data = f_d;
H1 = semilogx(Re_data, f_data, 'o', 'MarkerFaceColor', 'r'); % Experimental points
hold on

% Fit experimental Colebrook-White function
x1 = f_data;
y1 = Re_data;
func = @(p, x) 9.35 ./ (sqrt(x) .* (-p + sqrt(10 .^ (1.14 - 1 ./ sqrt(x)))));
P = fminsearch(@(p) norm(y1 - func(p, x1)), 0.01);

f3 = zeros(1, N);
syms x
for i = 1:N
    equation3 = (1 / sqrt(x)) - 1.14 + 2 * log10(P + 9.35 / (Re(i) * sqrt(x))) == 0;
    f3(i) = solve(equation3, x);
end
roughness = string(P);
H2 = semilogx(Re, f3, 'k', 'linewidth', 3); % Fitting curve
hold on

% Calculate friction factors for Moody Diagram
for k = 1:L
    for i = 1:N
        equation = (1 / sqrt(x)) - 1.14 + 2 * log10(e_D(k) + 9.35 / (Re(i) * sqrt(x))) == 0;
        f(i) = solve(equation, x);
    end
    semilogx(Re, f, 'linewidth', 1.5);
    hold on
    grid on
    if e_D(k) == 0
        txt = strcat('Smooth pipe: e/D =', string(e_D(k)));
        text(2 * 1e9, f(end), txt, 'FontSize', 10);
    else
        txt = strcat('e/D =', string(e_D(k)));
        text(2 * 1e9, f(end), txt, 'FontSize', 10);
    end
end

% Plot friction factor for laminar flow
Re_lam = 600:1400:2000;
f_lam = 64 ./ Re_lam;
semilogx(Re_lam, f_lam, '-ko', 'linewidth', 2);
text(1000, 0.095, 'f=64/Re Laminar flow', 'FontSize', 10);

% Transition region
yl = [0.03 0.09];
xl = [2000 3000];
xBox = [xl(1), xl(1), xl(2), xl(2), xl(1)];
yBox = [yl(1), yl(2), yl(2), yl(1), yl(1)];
patch(xBox, yBox, 'black', 'FaceColor', '#525260', 'FaceAlpha', 0.1);
text(3e3, 0.08, '\leftarrow Transition Region');

% Adjustments
xlim([100 7e10]);
ylim([0 0.1]);
xlabel('\textbf{\emph Re}', 'fontsize', 16, 'Interpreter', 'latex');
ylabel('\textbf{\emph f}', 'fontsize', 16, 'Interpreter', 'latex');
title('\textbf{Moody Diagram}', 'fontsize', 16, 'Interpreter', 'latex');
legend([H1, H2], ['Experimental data', strcat('Fitted data & roughness =', roughness)])

