%% Exercise 2

close all;
% Item 1

S = [-1 1 -2 1; 1 1 1 -1; 1 3 2 -1]
% cancel out x1 in the 2nd and 3rd rows
S(2,:) = S(2,:) + S(1,:)
S(3,:) = S(3,:) + S(1,:)

%% S =
%% 
%%   -1   1  -2   1
%%    0   2  -1   0
%%    0   4   0   0


% cancel out x2 in the last row
S(3,:) = S(3,:) - 2* S(2,:)

%% S = 
%% 
%%   -1   1  -2   1
%%    0   2  -1   0
%%    0   0   2   0
   
%% S is now in echelon form
% 2x3 = 0
% 2x2 - x3 = 0 -> x2 = 0
% -x1 + x2 - 2x3 = 1 ->  -x1 = 1
% x1 = -1

%% Item 2 LU Decomposition
A = [-1 1 -2; 1 1 1; 1 3 2]
b = [1; -1; -1]
L = eye(3);
% set the multipliers
L(2, 1) = -1;
L(3, 1) = -1;
L(3, 2) = 2;
U = [-1 1 -2; 0 2 -1; 0 0 2];

% Solve Ly = b
%% L =
%% 
%%    1   0   0
%%   -1   1   0
%%   -1   2   1
% y1 = b(1) = 1
% -y1+ y2 = b(2) = -1 -> y2=0
% -y1+ 2y2 + y3 = b(3) = -1 -> y3 = 0 

% Solve Ux = y
% [-1 1 -2; 0 2 -1; 0 0 2]x = [1; 0; 0]
% 2x3 = y3 = 0 -> x3 = 0
% 2x2 - x3 = y2 = 0 -> x2 = 0
% -x1 + x2 - 2x3 = y1 = 1 -> x1 = -1 

%% Item 3
A = [1 -3 5; 2 -4 3; 0 1 -1]
b = [1; -1; 3]
x_left_div = A \ b
x_row_reduce = rref([A b])
% The results of row reduce and left division are the same
% As mentioned in 2.1.3 LU decomposition is the basis of left division

%% Item 4
[L U P] = lu(A)
P*b
%% Solve Ly = Pb
%% P*b =
%% 
%%   -1
%%    1
%%    3
%% L =
%% 
%%    1.00000   0.00000   0.00000
%%    0.50000   1.00000   0.00000
%%    0.00000  -1.00000   1.00000
% y1 = -1
% 0.5y1 + y2 = 1 -> y2 = 1.5
% -y2 + y3 = 3 -> y3 = 4.5
%% Solve Ux = y
%% U =
%% 
%%    2.00000  -4.00000   3.00000
%%    0.00000  -1.00000   3.50000
%%    0.00000   0.00000   2.50000
% 2.5x3 = y3 = 4.5 -> x3 = 1.8
% -x2 + 3.5x3 = y2 = 1.5 -> x2 = 4.8
% 2x1 - 4x2 + 3x3 = y1 = -1 -> x1 = 6.4

%% Item 5
A = [1 -3 5; 2 -4 3; 0 1 -1]
b = [1; -1; 3]
A(1,:) = 3*A(2,:) - 4*A(3,:)
x_left_div = A \ b
x_row_reduce = rref([A b])
%% warning: matrix singular to machine precision
%% warning: called from
%%     lab2 at line 97 column 12
%% x_left_div =
%% 
%%    1.213675
%%   -0.034188
%%   -0.572650
%% 
%% x_row_reduce =
%% 
%%    1.00000   0.00000  -0.50000   0.00000
%%    0.00000   1.00000  -1.00000   0.00000
%%    0.00000   0.00000   0.00000   1.00000
% Left division result is no longer the same as the row reduce result

%% Item 6 Large system
randM = floor(rand(25,26)*10);
x_left_div = randM(:, 1:25) \ randM(:,26)
x_row_reduce = rref(randM)(:,26)

%% Item 7
time = linspace(0, 120, 13);
altitude = [
  7.00
  938.00
  4160.00
  9872.00
  17635.00
  26969.00
  37746.00
  50548.00
  66033.00
  83966.00
  103911.00
  125512.00
  147411.00
];
% Solve the quadratic polynomial
% A(t) = a1*t^2 + a2*t + a3
T = [time'.^2 time' ones(numel(time),1)];
a = T \ altitude
% Plot the parabola
altitude_fit = a(1).*time.^2 + a(2).*time + a(3);
figure;plot(time, altitude, 'o', time, altitude_fit , 'r-', 'linewidth', 2);
grid on;
legend('Data Values', 'Left Division Fit');
% First derivative: velocity v(t)  = 2*a1*t + a2
% velocity at 2 mins into flight is v(t=2mins=120sec) 
v2mins = 2*a(1)*120 
% v2mins =  2407.5 ft/s 
% Second derivative: acceleration a(t) = 2*a1

%% Item 8
x_year = [
  1900
  1910
  1920
  1930
  1940
  1950
  1960
  1970
  1980
  1990
  2000
  2010
];
y_population = [
  1650
  1750
  1860
  2070
  2300
  2525
  3018
  3682
  4440
  5310
  6127
  6930
];
% 8a Plot year vs Y=ln(y_population)
Y = log(y_population);
figure; plot(x_year, Y); 
% 8b Use polyfit
P=polyfit(x_year, Y,1);
Y_polyfit = polyval(P, x_year);
figure; plot(x_year, Y, x_year, Y_polyfit);
grid on;
title('Semilog plot');
xlabel('Year'); ylabel('ln(population)');
legend('original', 'polyfit');
%% Y = P(1)x + P(2) -> A = P(2)
% 8c 
C = exp(P(2));
y_m = C*exp(P(1)*x_year);
figure; plot(x_year, y_population, x_year, y_m);
grid on;
title('Exponential plot');
legend('original', 'polyfit');
xlabel('Year'); ylabel('population in millions');

%8d

% Y = P(1)x_year + P(2) = log(7e9)
year_7b = ( log(7000) - P(2) ) / P(1)
% year_7b =  2014.1

%% Item 9
D = [
-4 -3 -3 -2 -2 2 3 3 -4;
-1 1 2 2 3 3 2 -1 -1
];
x=D(1,:);
y=D(2,:);
% Rotate with R1 by 30 deg
theta1 = 30 * pi / 180;
R1 = [ cos(theta1) -sin(theta1); sin(theta1) cos(theta1) ] ; 
rotated30D = R1*D;

% Rotate with R2 by 135 deg
theta2 = 135 * pi / 180;
R2 = [ cos(theta2) -sin(theta2); sin(theta2) cos(theta2) ] ; 
rotated135D = R2*D;
figure;
plot(D(1,:), D(2,:), 'bo-', rotated30D(1,:), rotated30D(2,:), 'ro-',
  rotated135D(1,:), rotated135D(2,:), 'go-');
title('Matrix Rotations');
legend('original', 'rotated 90 degrees', 'rotated 135 degrees');
grid on;

%% Item 10
% dilate D by 2x
T = [2 0; 0 2];
dilatedD = T*D;
% reflect D along x, where angle l=0deg wrt to x axis
Ref = [cos(0) sin(0) ; sin(0) -cos(0)];
reflectedD = Ref*D;
figure;
plot(D(1,:), D(2,:), 'bo-', dilatedD(1,:), dilatedD(2,:), 'ro-',
  reflectedD(1,:), reflectedD(2,:), 'go-');
title('Matrix Transformations');
legend('original', 'expand 2x', 'reflect in x axis');
grid on;

%% Item 11
T=eye(3);
T(1,3) = -3; % 3 units left
T(2,3) = 2; % 2 units up
%% T =
%% 
%%    1   0  -3
%%    0   1   2
%%    0   0   1
% D1 -> append ones below y row of D
D1=[D(1,:); D(2,:);ones(1,numel(D(1,:)))]; 
translatedD = T*D1;
figure;
plot(D(1,:), D(2,:), 'bo-', translatedD(1,:), translatedD(2,:), 'ro-');
title('Matrix translation');
legend('original', 'translated');
grid on;

%% Item 12
D = [ 1 1 3 3 2 1 3 ; 2 0 0 2 3 2 2 ];
D1 = [D;ones(1,numel(D(1,:)))]; 
theta = 90 * pi / 180;
R = [ cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1] ; 
T=eye(3);
T(1,3) = 2;
T(2,3) = 1;
combinedT=inv(T)*R*T;
rotatedD = combinedT*D1;

figure;
plot(2,1, 'g*', D(1,:), D(2,:), 'bo-', rotatedD(1,:), rotatedD(2,:), 'ro-');
title('Rotate House Graph about (2,1)');
legend('(2,1)', 'original', 'rotated about (2,1)');
grid on;
axis([-10 5 -1 5], 'equal')
