# $Id: check_cplex.sh,v 1.1 2002/10/23 14:31:35 bzfpfend Exp $
CPLEXBIN=cplex75
BINNAME=$CPLEXBIN.$2
TSTNAME=`basename $1 .test`

OUTFILE=results/check.$TSTNAME.$BINNAME.out
ERRFILE=results/check.$TSTNAME.$BINNAME.err
TEXFILE=results/check.$TSTNAME.$BINNAME.tex
TMPFILE=results/check.$TSTNAME.$BINNAME.tmp

date >$OUTFILE
date >$ERRFILE

for i in `cat $1`
do
    echo @01 $i ===========
    echo @01 $i =========== >>$ERRFILE
    if [ -f $i ]
    then
	echo read $i                          > $TMPFILE
	echo set mip limits nodes 10000000      >> $TMPFILE
	echo set timelimit 3600                 >> $TMPFILE
	echo set clocktype 1                    >> $TMPFILE
	echo set mip interval 10000             >> $TMPFILE
	echo set mip tolerances absmipgap 1e-10 >> $TMPFILE
	echo set mip tolerances mipgap 0.0      >> $TMPFILE
	echo optimize                           >> $TMPFILE
	echo quit                               >> $TMPFILE
	$CPLEXBIN < $TMPFILE 2>>$ERRFILE
    else
	echo FILE NOT FOUND
	echo FILE NOT FOUND >>$ERRFILE
    fi
    echo =ready=
done | tee -a $OUTFILE

rm $TMPFILE

date >>$OUTFILE
date >>$ERRFILE

awk -f check_cplex.awk $OUTFILE > $TEXFILE
 
