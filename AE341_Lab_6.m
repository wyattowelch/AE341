% ------------------------------------------
% 
% AE341 Lab 6 - Wyatt Welch
%
% ------------------------------------------

% Data Gathered 

h12 = [140 125 110 95 80 65 50 35 20] / 1000;
h34 = [90  75  65  60 50 40 35 25 10] / 1000;
h56 = [35  30  26  25 20 15 15 10 5 ] / -1000;
h78 = [90  80  70  65 55 40 35 25 15] / 1000;
h90 = [115 95  85  80 70 50 40 30 10] / 1000;

V1 = [25 25 15 15 15 15 15 15 15] / 1000;
V2 = [35 35 25 25 25 25 25 25 25] / 1000;
t1 = [45.94 49.15 31.85 33.08 37.08 41.14 45.96 54.64 81.48];
t2 = [63.20 68.80 52.71 54.84 61.20 68.12 77.90 91.37 135.3];

D1 = .0225;
A1 = (D1/2)^2 * pi;
D2 = .0296;
A2 = (D2/2)^2 * pi;
Rs = .0125;
Rl = .0534;
g = 9.806;

% Calculations

Q1 = V1 ./ t1;
Q2 = V2 ./ t2;
Q = .5 .* (Q1 + Q2);

v1 = Q ./ A1;
v2 = Q ./ A2;

v12g = (v1 .^2) ./ (2 * g);
v22g = (v2 .^2) ./ (2 * g);

hL56 = h12 + v12g - v22g;
hL78 = h78 + v22g - v12g;

Re = (v1 .* D1) ./ 1e-6;
K12 = h12 ./ v12g;
K34 = h34 ./ v12g;
K56 = hL56 ./ v12g;
K78 = hL78 ./ v12g;
K90 = h90 ./ v12g;

% Best Fits

c1 = polyfit(h12,v12g,1);
xf1 = linspace(min(h12),max(h12), 9);
yf1 = polyval(c1, xf1);

c2 = polyfit(h34,v12g,1);
xf2 = linspace(min(h34),max(h34), 9);
yf2 = polyval(c2, xf2);

c3 = polyfit(h90,v12g,1);
xf3 = linspace(min(h90),max(h90), 9);
yf3 = polyval(c3, xf3);

c4 = polyfit(hL56,v12g,1);
xf4 = linspace(min(hL56),max(hL56), 9);
yf4 = polyval(c4, xf4);

c5 = polyfit(hL78,v12g,1);
xf5 = linspace(min(hL78),max(hL78), 9);
yf5 = polyval(c5, xf5);

c6 = polyfit(K12,Re,1);
xf6 = linspace(min(K12),max(K12), 9);
yf6 = polyval(c6, xf6);

c7 = polyfit(K34,Re,1);
xf7 = linspace(min(K34),max(K34), 9);
yf7 = polyval(c7, xf7);

c8 = polyfit(K56,Re,1);
xf8 = linspace(min(K56),max(K56), 9);
yf8 = polyval(c8, xf8);

c9 = polyfit(K78,Re,1);
xf9 = linspace(min(K78),max(K78), 9);
yf9 = polyval(c9, xf9);

c0 = polyfit(K90,Re,1);
xf0 = linspace(min(K90),max(K90), 9);
yf0 = polyval(c0, xf0);

% Plots

figure(1)
scatter(h12, v12g, 25, 'r*')
hold on, grid on
scatter(h34, v12g, 25, 'bo')
scatter(h90, v12g,25,'g+')
plot(xf1,yf1,"color", 'r')
plot(xf2,yf2,"color", 'b')
plot(xf3,yf3,"color", 'g')
legend('Mitre', 'Elbow', 'Bend', 'Location', "northwest")
xlabel('Total Head Loss [m]')
ylabel('v^2 / 2g [m]')
title('Total Head Loss v. v^2 / 2g')

figure(2)
scatter(hL56, v12g, 25, 'r*')
hold on, grid on
plot(xf4,yf4,"color", 'r')
legend('Enlargement', 'Location', "northwest")
xlabel('Total Head Loss [m]')
ylabel('v^2 / 2g [m]')
title('Total Head Loss v. v^2 / 2g')

figure(3)
scatter(hL78, v12g, 25, 'r*')
hold on, grid on
plot(xf5, yf5, 'color', 'r')
legend('Contraction', 'Location', "northwest")
xlabel('Total Head Loss [m]')
ylabel('v^2 / 2g [m]')
title('Total Head Loss v. v^2 / 2g')

figure(4)
scatter(K12,Re,'r*')
hold on, grid on
scatter(K34,Re,'bo')
scatter(K56,Re,"g+")
scatter(K78,Re,'black','x')
scatter(K90,Re,'magenta','d')
plot(xf6, yf6, 'color', 'r')
plot(xf7, yf7, 'color', 'b')
plot(xf8, yf8, 'color', 'g')
plot(xf9, yf9, 'color', 'black')
plot(xf0, yf0, 'color', 'magenta')
legend('Mitre', 'Elbow', 'Enlargement','Contraction', 'Bend', ...
    'Location', "northwestoutside")
xlabel('Loss Coefficient (K)')
ylabel('Reynold''s'' Number (Re)')
title("K v. Re")