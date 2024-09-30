clc, clear

% ----------------------------
% Data
% ----------------------------

nrun = [1;2;3;4;5;6;7;8;9;10];

Ajockeyd = [10;20;30;40;50;60;70;80;90;100];
Ajockeyd = Ajockeyd / 1000;
Alitre1 = [15;15;15;15;15;15;15;15;25;25];
Alitre2 = [25;25;25;25;25;25;25;25;35;35];

Bjockeyd = [20;40;60;80;100;120;140;160;180;200];
Bjockeyd = Bjockeyd / 1000;
Blitre1 = [15;15;15;15;15;15;15;15;15;15];
Blitre2 = [25;25;25;25;25;25;25;25;25;25];

Asec1 = [67.27;54.14;50.98;42.26;41.78;39.02;32.12;34.93;53.35;50.77];
Asec2 = [115.67;90.36;82.89;70.08;66.29;61.83;52.55;58.31;72.36;69.23];
AsecSplit = [48.40;36.22;31.91;27.82;24.51;22.81;20.43;23.38;19.01;18.46];

Bsec1 = [76.31;57.68;50.29;45.91;41.04;37.66;36.07;35.48;34.12;34.34];
Bsec2 = [128.69;96.24;80.64;72.4;65.2;59.83;56.46;54.99;51.28;51.65];

Jock = .6;
NozA = pi*((.01/2)^2);
Den = 998;
SW = 9975.41;
G = 9.795;
h = .035;
ABeta = 90;
BBeta = 180;
n = 10;

% ----------------------------
% Calculations
% ----------------------------

%A
Au1 = (Alitre1/NozA) ./ Asec1;
Au2 = (Alitre2/NozA) ./ Asec2;
Au = (Au1 + Au2) / 2;
Au0 = sqrt((Au .^2) - 2 * G * h) * 1e-6;
AmDot = Den * NozA * Au;
AFt = (AmDot .* Au0) .* (1 - cosd(ABeta));
AFe = 4 * G * Ajockeyd;
AmDu0 = AmDot .* Au0;
cc_Ae = (n * sum(Au0 .* AFe) - sum(Au0) * sum(AFe)) / (sqrt(n * sum(Au0 .^2) ...
    - sum(Au0) ^2) * sqrt(n * sum(AFe .^2) - sum(AFe) ^2))

%B
Bu1 = (Blitre1/NozA) ./ Bsec1;
Bu2 = (Blitre2/NozA) ./ Bsec2;
Bu = (Bu1 + Bu2) / 2;
Bu0 = sqrt((Bu .^2) - 2 * G * h) * 1e-6;
BmDot = Den * NozA * Bu;
BFt = (BmDot .* Bu0) .* (1 - cosd(BBeta));
BFe = 4 * G * Bjockeyd;
BmDu0 = BmDot .* Bu0;
cc_Ae = (n * sum(Bu0 .* BFe) - sum(Bu0) * sum(BFe)) / (sqrt(n * sum(Bu0 .^2) ...
    - sum(Bu0) ^2) * sqrt(n * sum(BFe .^2) - sum(BFe) ^2))

% ----------------------------
% Plot
% ----------------------------

%A
figure(1)
scatter(AmDu0,AFt,"red","o","filled")
hold on
scatter(AmDu0,AFe,"blue","o","filled")
Ap = polyfit(AmDu0,AFe,1);
Af1 = polyval(Ap,AmDu0);
plot(AmDu0,Af1,'-','LineWidth',2)
grid on
xlabel('Mass Flow Rate x Velocity at Deflector (kgm/s^2)')
ylabel('Theoretical and Experimental Force (Ns^2/m^2)')
title("90 deg. Deflector: Mass Flow Rate and Velocity vs. Forces")
legend('Theoretical Force', 'Experimental Force', 'Least Square Fit', 'Location', 'NorthWest')
hold off

%B
figure(2)
scatter(BmDu0,BFt,"red","o","filled")
hold on
scatter(BmDu0,BFe,"blue","o","filled")
Bp = polyfit(BmDu0,BFe,1);
Bf1 = polyval(Bp,BmDu0);
plot(BmDu0,Bf1,'-','LineWidth',2)
grid on
xlabel('Mass Flow Rate x Velocity at Deflector (kgm/s^2)')
ylabel('Theoretical and Experimental Force (Ns^2/m^2)')
title("180 deg. Deflector: Mass Flow Rate and Velocity vs. Forces")
legend('Theoretical Force', 'Experimental Force', 'Least Square Fit', 'Location', 'NorthWest')
hold off