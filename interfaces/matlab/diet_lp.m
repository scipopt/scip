% Stigler's famous diet problem, as given in, e.g., George Dantzig, Chapter 27.1. In Linear Programming and
% Extensions. Princeton University Press, Princeton, New Jersey, 1963.
% In this version, the nutrition values are scaled per dollar spent to by a certain type of food.
% See also the GAMS model library: http://www.gams.com/modlib/libhtml/diet.htm

clear;

% the (dense) matrix with the nine nutrition values for each of the 20 foods
matrix = [
44.7,   36,  8.4, 20.6,  7.4, 15.7, 41.7,   2.2,  4.4, 5.8, 2.4,  2.6,  5.8, 14.3,  1.1,    9.6,  8.5, 12.8, 17.4, 26.9;...
1411,  897,  422,   17,  448,  661,  0.0,   333,  249, 705, 138,  125,  166,  336,   106,   138,   87,   99, 1055, 1691;...
 2.0,  1.7, 15.1,  0.6, 16.4,    1,  0.0,   0.2,  0.3, 6.8, 3.7,    4,  3.8,  1.8,   0.0,   2.7,  1.7,  2.5,  3.7, 11.4;...
 365,   99,    9,    6,   19,   48,  0.0,   139,   37,  45,  80,   36,   59,  118,   138,    54,  173,  154,  459,  792;...
 0.0, 30.9,   26, 55.8, 28.1,  0.0,  0.2, 169.2,  0.0, 3.5,  69,  7.2, 16.6,  6.7, 918.4, 290.7, 86.8, 85.7,  5.1,  0.0;...
55.4, 17.4,    3,  0.2,  0.8,  9.6,  0.0,   6.4, 18.2,   1, 4.3,    9,  4.7, 29.4,   5.7,   8.4,  1.2,  3.9, 26.9, 38.4;...
33.3,  7.9, 23.5,  0.0, 10.3,  8.1,  0.5,  50.8,  3.6, 4.9, 5.8,  4.5,  5.9,  7.1,  13.8,   5.4,  4.3,  4.3, 38.2, 24.6;...
 441,  106,   11,  0.0,    4,  471,    5,   316,   79, 209,  37,   26,   21,  198,    33,    83,   55,   65,   93,  217;...
 0.0,  0.0,   60,  0.0,  0.0,  0.0,  0.0,   525,  0.0, 0.0, 862, 5369, 1184, 2522,  2755,  1912,   57,  257,  0.0,  0.0
];

% nine linear constraints ensure that the daily requirements of nutrients are fulfilled
% all constraints are of the type a_1*x_1+...+a20*x_20 >= b
% In SCIP, constraints are represented as left hand side <= linear sum <= right hand side
% Hence, the left hand sides are all finite, the right hand sides are plus infinity
lhs = [3,70,0.8,12,5,1.8,2.7,18,75]';
rhs = [1e+20,1e+20,1e+20,1e+20,1e+20,1e+20,1e+20,1e+20,1e+20]';

% the lower bounds of all variables are zero (you cannot buy negative food),
% the upper bounds are plus infinity (you can buy arbitrarily many food)
lb = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]';
ub = [1e+20,1e+20,1e+20,1e+20,1e+20,1e+20,1e+20,1e+20,1e+20,1e+20,1e+20,1e+20,1e+20,1e+20,1e+20,1e+20,1e+20,1e+20,1e+20,1e+20]';

% the nutrition values in the matrix are scaled per dollar spent, hence, the objective is all one
% all variables are continuous (you can buy arbitrary fractions of food)
% the goal is to minimize the money spent, while all nutrition requirements are fulfilled
obj = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]';
vartype = ['c','c','c','c','c','c','c','c','c','c','c','c','c','c','c','c','c','c','c','c'];
objsense = 'min';

% call SCIP
[bestsol, objval] = scip(matrix, lhs, rhs, obj, lb, ub, vartype, objsense);
