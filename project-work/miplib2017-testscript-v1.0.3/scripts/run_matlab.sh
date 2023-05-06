#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2017            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

SOLVER=$1
BINNAME=$2
NAME=$3
TIMELIMIT=$4
MEMLIMIT=$5   # The memory limit (in MB) # currently unused
SOLFILE=$6
THREADS=$7
MIPGAP=$8

TMPBASE=tmp
TMPFILE=tmp.m

FILENAME=$(basename "$NAME")
EXTENSION="${FILENAME##*.}"
PROBNAME="${FILENAME%.*}"
if [ "$EXTENSION" = "gz" ] ; then
   gunzip -c $NAME > $PROBNAME
else
   $PROBNAME=$NAME
fi

echo > $TMPFILE
echo > $SOLFILE

echo try                                                                          >>$TMPFILE
echo "fprintf('Release R%s\n',version('-release'))"                               >>$TMPFILE
echo "p = mpsread('$PROBNAME')"                                                   >>$TMPFILE

echo "p.options.RelativeGapTolerance = $MIPGAP;"                                  >>$TMPFILE
echo "p.options.MaxTime = $TIMELIMIT;"                                            >>$TMPFILE
echo "p.options.MaxNodes = Inf;"                                                  >>$TMPFILE
echo "p.options.ObjectiveImprovementThreshold = 0.0;"                             >>$TMPFILE
echo "p.options"                                                                  >>$TMPFILE
echo "[x,fval,exitflag,output] = intlinprog(p);"                                  >>$TMPFILE

echo "output"                                                                     >>$TMPFILE
echo "fprintf('BBnodes %d\n',output.numnodes);"                                   >>$TMPFILE
echo "if exitflag > 0"                                                            >>$TMPFILE
echo "fprintf('PrimalBound %.16g\n',fval);"                                       >>$TMPFILE
echo "fprintf('DualBound %.16g\n',fval-output.absolutegap);"                      >>$TMPFILE
echo "end"                                                                        >>$TMPFILE

echo "ncols = size(p.lb,1);"                                                      >>$TMPFILE
echo "fileid = fopen('$PROBNAME');"                                               >>$TMPFILE
echo "tline = fgetl(fileid);"                                                     >>$TMPFILE
echo "while length(tline) == 0 || (ischar(tline) && ~strcmp(tline(1),'C'))"       >>$TMPFILE
echo "    tline = fgetl(fileid);"                                                 >>$TMPFILE
echo "end"                                                                        >>$TMPFILE
echo "tline = fgetl(fileid);"                                                     >>$TMPFILE
echo "columnname = strings(ncols,1);"                                             >>$TMPFILE
echo "prevname = \"\";"                                                           >>$TMPFILE
echo "i = 0;"                                                                     >>$TMPFILE
echo "while length(tline) == 0 || (ischar(tline) && strcmp(tline(1),' ')) "       >>$TMPFILE
echo "    if length(tline) > 0 && ~strcmp(tline(1),'*')"                          >>$TMPFILE
echo "        fields = textscan(tline, ' %s');"                                   >>$TMPFILE
echo "        thisname = string(fields{1}{1});"                                   >>$TMPFILE
echo "        thistype = string(fields{1}{2});"                                   >>$TMPFILE
echo "        if ~strcmp(thisname,prevname) && ~strcmp(thistype,\"'MARKER'\")"    >>$TMPFILE
echo "            i = i + 1;"                                                     >>$TMPFILE
echo "            columnname(i) = thisname;"                                      >>$TMPFILE
echo "            prevname = thisname;"                                           >>$TMPFILE
echo "        end"                                                                >>$TMPFILE
echo "    end"                                                                    >>$TMPFILE
echo "    tline = fgetl(fileid);"                                                 >>$TMPFILE
echo "end"                                                                        >>$TMPFILE
echo "assert(ncols == i)"                                                         >>$TMPFILE
echo "fclose(fileid);"                                                            >>$TMPFILE
echo "solfileid = fopen('$SOLFILE','w');"                                         >>$TMPFILE
echo "if exitflag > 0"                                                            >>$TMPFILE
echo "    fprintf(solfileid,'%-20s %.16g\n','=obj=',fval);"                       >>$TMPFILE
echo "    for j=1:size(x)"                                                        >>$TMPFILE
echo "        fprintf(solfileid,'%-20s %.16g\n',columnname(j),x(j));"             >>$TMPFILE
echo "    end"                                                                    >>$TMPFILE
echo "else"                                                                       >>$TMPFILE
echo "    fprintf(solfileid,'=infeas=\n');"                                       >>$TMPFILE
echo "end"                                                                        >>$TMPFILE
echo "fclose(solfileid);"                                                         >>$TMPFILE
echo "catch ME"                                                                   >>$TMPFILE
echo "disp(ME.message)"                                                           >>$TMPFILE
echo "end"                                                                        >>$TMPFILE
echo "exit"                                                                       >>$TMPFILE

MATLABOPTIONS="-nojvm -nodisplay -nosplash"

$BINNAME  $MATLABOPTIONS -r $TMPBASE

# remove tmp file
rm $TMPFILE
if [ "$NAME" != "$PROBNAME" ]; then
   rm $PROBNAME
fi
