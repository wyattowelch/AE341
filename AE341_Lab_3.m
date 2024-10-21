%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Wyatt Welch - Flow Through an Orifice
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data

L1 = [15; 15; 15; 15; 15; 15; 15; 15; 15; 15; 15; 15; 15; 15]; % L
L2 = [25; 25; 25; 25; 25; 25; 25; 25; 25; 25; 25; 25; 25; 25];                            
t1 = [61.58; 61.54; 62.13; 61.83; 61.92; 67.87; 71.39; 75.18; 80.27; 86.76; 93.95; 105.83; 112.01; 122.78]; % s
t2 = [102.12; 101.08; 104.67; 102.77; 103.02; 112.64; 118.42; 124.91; 134.28; 143.39; 156.01; 175.73; 185.55; 203.37];

dc = [9.2; 9.52; 10.56; 11.36; 9.05] / 1000; % m
do = 0.013;
dp = .026;

hc = ([386; 386; 386; 386; 386] + 45) / 1000;
ho = [396.5; 396.5; 396.5; 396.5; 396.5; 361.5; 326.5; 291.5; 256.5; 221.5; 186.5; 151.5; 129.5; 110.5] / 1000;
ho5 = [396.5; 396.5; 396.5; 396.5; 396.5] / 1000;

Ao = .25 * pi * (.013^2); % m^2
g = 9.795; % m/s^2

% Calculations

Qa = L1./t1 / 1000;
Qb = L2./t2 / 1000;
Qexp = .5 * (Qa+Qb);

uo = sqrt(2*g*ho);
uo5 = sqrt(2*g*ho5);
uc = sqrt(2*g*hc);

% Coefficient of Velocities
Cu = uc./uo5;
Cum = mean(uc./uo5);

% Coefficient of Contraction
Cc = dc/do;
Ccm = mean(dc/do);

Ac = Ccm * Ao;
Qc = Ac * uc;
Qo = uo * Ao;

% Discharge Coefficient 
Cdpm = Cum * Ccm;
Cdp = Cu .* Cc;
Cd = Qexp ./ Qo;

K = Cd * Ao * sqrt(2 * g);

% Plot

y = log(Qexp(6:end));
x = log(ho(6:end));
b = log(K(6:end));
X = .5 .* x + b;

slope_and_bias = polyfit(x, y, 1);
slope = slope_and_bias(1);
bias = slope_and_bias(2);

figure(1)
scatter(x,y, 'r')
hold on
plot(x, slope .* x + bias,'b - ' ) 
grid on 
xlabel('Log(h0) (m)')
ylabel('Log(Qexp) (m^3/s)')
title('Log(h0) vs. Log(Qexp)')
legend('Calculated', 'Best Fit', 'Location', 'northwest')



% The Alamo
% 9/11