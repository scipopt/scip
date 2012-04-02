%* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%*                                                                           *
%*                  This file is part of the program and library             *
%*         SCIP --- Solving Constraint Integer Programs                      *
%*                                                                           *
%*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            *
%*                            fuer Informationstechnik Berlin                *
%*                                                                           *
%*  SCIP is distributed under the terms of the ZIB Academic Licence.         *
%*                                                                           *
%*  You should have received a copy of the ZIB Academic License              *
%*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
%*                                                                           *
%* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

%@file    scip.m
%@brief   SCIP installation file for MATLAB interface
%@author  Ambros Gleixner
%@author  Timo Berthold

function[bestsol, objval] = scip(matrix, lhs, rhs, obj, lb, ub, vartype, objsense)

% copyright %
fprintf('SCIP-MATLAB interface\n');
fprintf('Copyright (c) 2002-2012 Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)\n\n');

% syntax description %
usage = '\n';
usage = [usage 'syntax: [bestsol, objval] = scip(matrix, lhs, rhs, obj, lb, ub, vartype, objsense)\n'];
usage = [usage 'input parameters are\n'];
usage = [usage '   matrix   : real matrix with constraint nonzeros\n'];
usage = [usage '   lhs      : real vector with left-hand sides of constraints\n'];
usage = [usage '   rhs      : real vector with right-hand sides of constraints\n'];
usage = [usage '   obj      : real vector with objective coefficient of variables\n'];
usage = [usage '   lb       : real vector with lower bounds of variables\n'];
usage = [usage '   ub       : real vector with upper bounds of variables\n'];
usage = [usage '   vartype  : character vector with types of variables (b: binary, i: general integer, c: continuous)\n'];
usage = [usage '   objsense : value indicating optimization sense (+1: minimize, -1: maximize)\n'];
usage = [usage 'output parameters are\n'];
usage = [usage '   bestsol  : real vector containing best solution found\n'];
usage = [usage '   objval   : objective function value of best solution found\n'];
usage = [usage '\n'];

% output is undefined for an invalid call %
bestsol = NaN;
objval = NaN;

% we do not define default values for input parameters %
if( nargin ~= 8 )
   fprintf('\ninvalid call: provide all input parameters\n');
   fprintf(usage);
   return;
end

% display input parameters %
matrix
lhs
rhs
obj
lb
ub
vartype
objsense

% matrix must be two-dimensional %
if( ndims(matrix) ~= 2 )
   fprintf('\ninvalid call: ndims(matrix) ~= 2\n');
   fprintf(usage);
   return;
end

% get number of variables and constraints %
[nconss,nvars] = size(matrix);

% check that vectors are one-dimensional %
if( ~isvector(lhs) )
   fprintf('\ninvalid call: ~isvector(lhs)\n');
   fprintf(usage);
   return;
elseif( ~isvector(rhs) )
   fprintf('\ninvalid call: ~isvector(rhs)\n');
   fprintf(usage);
   return;
elseif( ~isvector(obj) )
   fprintf('\ninvalid call: ~isvector(obj)\n');
   fprintf(usage);
   return;
elseif( ~isvector(lb) )
   fprintf('\ninvalid call: ~isvector(lb)\n');
   fprintf(usage);
   return;
elseif( ~isvector(ub) )
   fprintf('\ninvalid call: ~isvector(ub)\n');
   fprintf(usage);
   return;
elseif( ~isvector(vartype) )
   fprintf('\ninvalid call: ~isvector(vartype)\n');
   fprintf(usage);
   return;
end

% check that vectors have consistent length %
if( length(lhs) ~= nconss )
   fprintf('\ninvalid call: length(lhs) ~= size(matrix,1)\n');
   fprintf(usage);
   return;
elseif( length(rhs) ~= nconss )
   fprintf('\ninvalid call: length(rhs) ~= size(matrix,1)\n');
   fprintf(usage);
   return;
elseif( length(obj) ~= nvars )
   fprintf('\ninvalid call: length(obj) ~= size(matrix,2)\n');
   fprintf(usage);
   return;
elseif( length(lb) ~= nvars )
   fprintf('\ninvalid call: length(lb) ~= size(matrix,2)\n');
   fprintf(usage);
   return;
elseif( length(ub) ~= nvars )
   fprintf('\ninvalid call: length(ub) ~= size(matrix,2)\n');
   fprintf(usage);
   return;
elseif( length(vartype) ~= nvars )
   fprintf('\ninvalid call: length(vartype) ~= size(matrix,2)\n');
   fprintf(usage);
   return;
end

% check that input values have correct type %
if( ~isreal(matrix) )
   fprintf('\ninvalid call: ~isreal(matrix)\n');
   fprintf(usage);
   return;
elseif( ~isreal(lhs) )
   fprintf('\ninvalid call: ~isreal(lhs)\n');
   fprintf(usage);
   return;
elseif( ~isreal(rhs) )
   fprintf('\ninvalid call: ~isreal(rhs)\n');
   fprintf(usage);
   return;
elseif( ~isreal(obj) )
   fprintf('\ninvalid call: ~isreal(obj)\n');
   fprintf(usage);
   return;
elseif( ~isreal(lb) )
   fprintf('\ninvalid call: ~isreal(lb)\n');
   fprintf(usage);
   return;
elseif( ~isreal(ub) )
   fprintf('\ninvalid call: ~isreal(ub)\n');
   fprintf(usage);
   return;
elseif( ~ischar(vartype) )
   fprintf('\ninvalid call: ~ischar(vartype)\n');
   fprintf(usage);
   return;
elseif( ~strcmp(objsense,'min') && ~strcmp(objsense,'max') )
   fprintf('\ninvalid call: objsense must be min or max\n');
   fprintf(usage);
   return;
else
   for v = 1:length(vartype)
      if( vartype(v) ~= 'b' && vartype(v) ~= 'i' && vartype(v) ~= 'c' )
         fprintf('\ninvalid call: vartype(%d) must be one of {b,i,c}\n', v);
         fprintf(usage);
         return;
      end
   end
end

% call SCIP %
[bestsol, objval] = matscip(matrix, lhs, rhs, obj, lb, ub, vartype, objsense);

bestsol
objval
