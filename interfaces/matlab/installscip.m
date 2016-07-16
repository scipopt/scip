%* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%*                                                                           *
%*                  This file is part of the program and library             *
%*         SCIP --- Solving Constraint Integer Programs                      *
%*                                                                           *
%*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            *
%*                            fuer Informationstechnik Berlin                *
%*                                                                           *
%*  SCIP is distributed under the terms of the ZIB Academic Licence.         *
%*                                                                           *
%*  You should have received a copy of the ZIB Academic License              *
%*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
%*                                                                           *
%* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

%@file    installscip.m
%@brief   SCIP installation file for MATLAB interface
%@author  Ambros Gleixner
%@author  Timo Berthold
%@author  Marc Pfetsch

fprintf('Starting installer for SCIP-Matlab interface.\n');
fprintf('Copyright (c) 2002-2016 Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)\n\n');

% determine directory char depending on environment
if ispc
   dirchar = '\';
else
   dirchar = '/';
end

% extract path from current position
pathtoscip = mfilename('fullpath');
slashpos = strfind(pathtoscip, dirchar);
pathtoscip = pathtoscip(1:slashpos(end-2));

% remove trailing slash
if pathtoscip(end) == dirchar
  pathtoscip(end) = [];
end

fprintf('Checking SCIP base path <%s> ...\n', pathtoscip);

% check whether we are in the correct position
if ( ~isempty(dir([pathtoscip '/src/scip/scip.h'])) )
  fprintf('Base path seems to be correct.\n');
else
  fprintf('Could not find file <%s>.\n', [pathtoscip '/src/scip/scip.h']);

  % if path did not work, ask for SCIP directory
  pathtoscip = '';
  while( ~ischar(pathtoscip) || isequal(pathtoscip, '') )
    pathtoscip = input('Specify path to SCIP installation: ', 's');

    if( ~ischar(pathtoscip) )
      fprintf('Incorrect path specified. Try again or press CTRL-C to abort.\n');
    elseif( isequal(pathtoscip, '') )
      fprintf('Empty path specified. Try again or press CTRL-C to abort.\n');
    end
  end
end

% detect SCIP source directory
fprintf('\nDetecting include directory ...');

pathtosrc = [pathtoscip '/src'];
pathtosciph = [pathtosrc '/scip/scip.h'];

if( ~isempty(dir(pathtosciph)) )
  fprintf('\t using <%s>.\n', pathtosrc);
else
   while( isempty(dir(pathtosciph)) )
      fprintf('Could not find header file <scip/scip.h> in <%s>.\n', pathtosrc);
      pathtosrc = input('Specify SCIP source directory manually or press CTRL-C to abort: ', 's');
      pathtosciph = [pathtosrc '/scip/scip.h'];
   end
end

% detect SCIP library
fprintf('Detecting SCIP library ...');

pathtolibscip = [pathtoscip '/lib/libscip.a'];

if( ~isempty(dir(pathtolibscip)) )
   fprintf('\t using SCIP library <%s>.\n', pathtolibscip);
else
   while( isempty(dir(pathtolibscip)) )
      fprintf('Could not find library <%s>.\n', pathtolibscip);
      pathtolibscip = input('Specify libscip.a manually or press CTRL-C to abort: ', 's');
   end
end

% detect LPI library
fprintf('Detecting LPI library ...');

pathtoliblpi = [pathtoscip '/lib/liblpispx.a'];

if( ~isempty(dir(pathtoliblpi)) )
   fprintf('\t using LPI library <%s>.\n', pathtoliblpi);
else
   while( isempty(dir(pathtoliblpi)) )
      fprintf('Could not find library <%s>.\n', pathtoliblpi);
      pathtoliblpi = input('Specify liblpispx.a manually or press CTRL-C to abort:', 's');
   end
end

% detect NLPI library
fprintf('Detecting NLPI library ...');

pathtolibnlpi = [pathtoscip '/lib/libnlpi.cppad.a'];

if( ~isempty(dir(pathtolibnlpi)) )
   fprintf('\t using NLPI library <%s>.\n', pathtolibnlpi);
else
   while( isempty(dir(pathtolibnlpi)) )
      fprintf('Could not find library <%s>.\n', pathtolibnlpi);
      pathtolibnlpi = input('Specify libnlpi.cppad.a manually or press CTRL-C to abort:', 's');
   end
end

% detect SoPlex library
fprintf('Detecting SoPlex library ...');

pathtolibsoplex = [pathtoscip '/lib/libsoplex.a'];

if ( ~isempty(dir(pathtolibsoplex)) )
   fprintf('\t using SoPlex library <%s>.\n', pathtolibsoplex);
else
  % try further possiblilties - start with linux x64 systems
  pathtolibsoplex = [pathtoscip '/lib/libsoplex.linux.x86_64.gnu.opt.a'];

  if ( ~isempty(dir(pathtolibsoplex)) )
    fprintf('\t using SoPlex library <%s>.\n', pathtolibsoplex);
  else
    while ( isempty(dir(pathtolibsoplex)) )
      fprintf('Could not find library <%s>.\n', pathtolibsoplex);
      pathtolibsoplex = input('Specify libsoplex.a manually or press CTRL-C to abort: ', 's');
    end
  end
end

% compile
flags = ['-I"' pathtosrc '"' ' matscip.c ' pathtolibscip ' ' pathtoliblpi ' ' pathtolibnlpi ' ' pathtolibsoplex];
fprintf('\nCompiling with flags <%s>\n', flags);
eval(['mex ' flags]);
