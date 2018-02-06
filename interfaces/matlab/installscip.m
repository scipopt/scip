%* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%*                                                                           *
%*                  This file is part of the program and library             *
%*         SCIP --- Solving Constraint Integer Programs                      *
%*                                                                           *
%*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            *
%*                            fuer Informationstechnik Berlin                *
%*                                                                           *
%*  SCIP is distributed under the terms of the ZIB Academic License.         *
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

% setting whether we want to use ZIMPL - modify, if ZIMPL should be used
withzimpl = 1;

fprintf('Starting installer for SCIP-Matlab interface.\n');
fprintf('Copyright (C) 2002-2018 Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)\n\n');

% determine whether we are using OCTAVE
isoctave = (exist ('OCTAVE_VERSION', 'builtin') > 0);

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

fprintf('Detecting SCIP base path <%s> ...\n', pathtoscip);

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
   error('\nCould not find header file <scip/scip.h> in <%s>.', pathtosrc);
end

% detect SCIP library
fprintf('Detecting SCIP library ...');

pathtolibscip = [pathtoscip '/lib/libscip.so'];

if( ~isempty(dir(pathtolibscip)) )
   fprintf('\t using SCIP library <%s>.\n', pathtolibscip);
else
   error('\nCould not find library <%s>.', pathtolibscip);
end

% detect LPI library
fprintf('Detecting LPI library ...');

pathtoliblpi = [pathtoscip '/lib/liblpispx1.so'];
liblpi = '-llpispx1';

if( ~isempty(dir(pathtoliblpi)) )
   fprintf('\t using LPI library <%s>.\n', pathtoliblpi);
else
   % try old SoPlex interface
   pathtoliblpi = [pathtoscip '/lib/liblpispx.so'];

   if ( ~isempty(dir(pathtoliblpi)) )
      fprintf('\t using LPI library <%s>.\n', pathtoliblpi);
      liblpi = '-llpispx';
   else
      % try new SoPlex interface
      pathtoliblpi = [pathtoscip '/lib/liblpispx2.so'];

      if ( ~isempty(dir(pathtoliblpi)) )
	 fprintf('\t using LPI library <%s>.\n', pathtoliblpi);
	 liblpi = '-llpispx2';
      else
	 error('\nCould not find library <%s>.', pathtoliblpi);
      end
   end
end

% detect NLPI library
fprintf('Detecting NLPI library ...');

pathtolibnlpi = [pathtoscip '/lib/libnlpi.cppad.so'];

if( ~isempty(dir(pathtolibnlpi)) )
   fprintf('\t using NLPI library <%s>.\n', pathtolibnlpi);
else
   error('\nCould not find library <%s>.', pathtolibnlpi);
end

% detect SoPlex library
fprintf('Detecting SoPlex library ...');

pathtolibsoplex = [pathtoscip '/lib/libsoplex.so'];

% determine type of computer
arch = computer;

if ( ~isempty(dir(pathtolibsoplex)) )
   fprintf('\t using SoPlex library <%s>.\n', pathtolibsoplex);
else
  % if path has not been set, we try to determine SoPlex library from architecture
  % we guess that we are using a gnu compiler
  switch arch
  case 'PCWIN'
     error('\nDo not know how to determine soplex library on Windows.');
  case 'PCWIN64'
     error('\nDo not know how to determine soplex library on Windows.');
  case 'GLNXA64'
     pathtolibsoplex = [pathtoscip '/lib/libsoplex.linux.x86_64.gnu.opt.so'];
     libsoplex = ' -lsoplex.linux.x86_64.gnu.opt.so';
  case 'x86_64-pc-linux-gnu'
     pathtolibsoplex = [pathtoscip '/lib/libsoplex.linux.x86_64.gnu.opt.so'];
     libsoplex = ' -lsoplex.linux.x86_64.gnu.opt';
  case 'MACI64'
     error('\nDo not know how to determine soplex library on Macs.');
  otherwise
     error('\nUnsupported MATLAB architecture %s.', arch);
  end

  if ( ~isempty(dir(pathtolibsoplex)) )
     fprintf('\t using SoPlex library <%s>.\n', pathtolibsoplex);
  else
     error('\nCould not find library <%s>.', pathtolibsoplex);
  end
end

% possibly detect Zimpl library
if ( withzimpl == 1 )
   fprintf('Detecting Zimpl library ...');

   pathtolibzimpl = [pathtoscip '/lib/libzimpl.so'];
   libszimpl = ' -lzimpl';

   if ( ~isempty(dir(pathtolibzimpl)) )
      fprintf('\t using Zimpl library <%s>.\n', pathtolibzimpl);
   else
      % if path has not been set, we try to determine Zimpl library from architecture
      % we guess that we are using a gnu compiler
      switch arch
	 case 'PCWIN'
	    error('\nDo not know how to determine Zimpl library on Windows.');
	 case 'PCWIN64'
	    error('\nDo not know how to determine Zimpl library on Windows.');
	 case 'GLNXA64'
	    pathtolibzimpl = [pathtoscip '/lib/libzimpl.linux.x86_64.gnu.opt.so'];
	    libzimpl = ' -lzimpl.linux.x86_64.gnu.opt.so';
	 case 'x86_64-pc-linux-gnu'
	    pathtolibzimpl = [pathtoscip '/lib/libzimpl.linux.x86_64.gnu.opt.so'];
	    libzimpl = ' -lzimpl.linux.x86_64.gnu.opt';
	 case 'MACI64'
	    error('\nDo not know how to determine soplex library on Macs.');
	 otherwise
	    error('\nUnsupported MATLAB architecture %s.', arch);
      end

      if ( ~isempty(dir(pathtolibzimpl)) )
	 fprintf('\t using Zimpl library <%s>.\n', pathtolibzimpl);
      else
	 error('\nCould not find library <%s>.', pathtolibzimpl);
      end
   end
else
   libzimpl = '';
end

% compile
if ( ~isoctave )
   pathtosciplib = [pathtoscip '/lib'];
   flags = ['-I' pathtosrc ' -L' pathtosciplib ' -lscip -lnlpi.cppad.so ' liblpi libzimpl libsoplex ' -lgmp -lz LDFLAGS=''$LDFLAGS -Wl,-rpath,' pathtosciplib ''''];
   fprintf('\nRunning:\nmex %s matscip.c\n', flags);
   eval(['mex ' flags ' matscip.c']);
else
   pathtosciplib = [pathtoscip '/lib'];
   flags = ['"-Wl,-rpath=' pathtosciplib '" -I' pathtosrc ' -L' pathtosciplib ' -lscip -lnlpi.cppad ' liblpi libzimpl libsoplex ' -lgmp -lz'];
   fprintf('\nRunning:\nmkoctfile --mex %s matscip.c\n', flags);
   eval(['mkoctfile --mex ' flags ' matscip.c'], 'printf ("The following error occurred:\n%s\n", lasterr());');
end

if ( isempty(dir(['matscip.' mexext])) )
   fprintf('\nBuild failed.\n');
   fprintf('\nPossibly you have to compile SCIP with make USRCFLAGS="-fPIC" USRCXXFLAGS="-fPIC"\n');
   if ( isoctave )
      fprintf('\nCurrently the error reporting of mkoctfile does not seem to work in Octave.\n');
      fprintf('If there are any problems, try to manually enter the above command.\n');
   end
else
   fprintf('\nBuild complete.\n');
end
