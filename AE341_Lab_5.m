% ------------------------------
%
% AE341 Lab 5 - Wyatt Welch
%
% ------------------------------

% Data collected

%   Rectangle

rv1 = [25 25 25 25 15 15 15  5  5] / 1e3; % m^3
rv2 = [35 35 35 35 25 25 25 15 15] / 1e3; % m^3
rt1 = [29.0 34.4 39.6 45.0 30.6 35.2 47.0 16.9 29.0]; % s
rt2 = [41.0 47.6 55.3 62.8 51.5 58.2 79.3 63.8 85.2]; % s
rho = 0.627; % m
rHpHo = linspace(.57, .61, 9); % m
rH = rho - rHpHo;

%   Triangle

tv1 = [25 25 15 15 15  5  5] / 1e3;
tv2 = [35 35 25 25 25 15 15] / 1e3;
tt1 = [30.3 38.1 30.1 38.0 49.5 28.1 30.9 ];
tt2 = [41.9 51.2 50.0 63.3 82.4 79.2 197.8];
tho = .623;
tHpHo = linspace(.546, .594, 7);
tH = tho - tHpHo;

g = 9.79; % m/s^2

% Calculations

%   Rectangle

rQ1 = rv1 ./ rt1;
rQ2 = rv2 ./ rt2;
rQe = .5 * (rQ1 + rQ2);
rQt = (2/3) * sqrt(2 * g) * (rH .^ 1.5) * .03;

rCdd = rQe ./ rQt
rCddi = rCdd'
rCd = rQe / rQt

%   Triangle

tQ1 = tv1 ./ tt1;
tQ2 = tv2 ./ tt2;
tQe = .5 * (tQ1 + tQ2);
tQt = (8/15) * sqrt(2 * g) * tand(30/2) * ((tH+.01) .^ 2.5); % the +.01 is wrong tee hee

tCdd = tQe ./ tQt
tCddi = tCdd'
tCd = tQe / tQt

% Plots

%   Rectangle

ry = log(rQe);
rx = log(rH);

rSaB = polyfit(rx, ry, 1);
rS = rSaB(1);
rB = rSaB(2);

figure(1)
scatter(rx, ry, 'r')
hold on, grid on
plot(rx, 1.5 .* rx-2.7)
plot(rx, rS .* rx + rB, 'b - ')
ylabel('Log(Qexp) [m^3/s]')
xlabel('Log(H) [m]')
title('Rectangular Weir - Log(H) v. Log(Qexp)')
legend('Empirical', 'Theoretical Fit (Slope = 1.5)', 'Best Fit (Slope = 1.26)', 'Location', 'northwest')

%   Triangle

ty = log(tQe);
tx = log(tH);

tSaB = polyfit(tx, ty, 1);
tS = tSaB(1);
tB = tSaB(2);

figure(2)
scatter(tx, ty, 'r')
hold on, grid on
plot(tx, 2.5 .* tx-.5)
plot(tx, tS .* tx + tB, 'b - ')
ylabel('Log(Qexp) [m^3/s]')
xlabel('Log(H) [m]')
title('Triangular Weir - Log(H) v. Log(Qexp)')
legend('Empirical', 'Theoretical Fit (Slope = 2.5)', 'Best Fit (Slope = 1.99)', 'Location', 'northwest')
