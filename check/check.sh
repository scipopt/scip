# $Id: check.sh,v 1.3 2003/11/21 09:18:58 bzfpfend Exp $
TSTNAME=$1
BINNAME=$2
SETNAME=$3
BINID=$4
TIMELIMIT=$5

OUTFILE=results/check.$TSTNAME.$BINID.$SETNAME.out
ERRFILE=results/check.$TSTNAME.$BINID.$SETNAME.err
RESFILE=results/check.$TSTNAME.$BINID.$SETNAME.res
TEXFILE=results/check.$TSTNAME.$BINID.$SETNAME.tex
TMPFILE=results/check.$TSTNAME.$BINID.$SETNAME.tmp
SETFILE=settings/$SETNAME.set

uname -a >$OUTFILE
uname -a >$ERRFILE
date >>$OUTFILE
date >>$ERRFILE

for i in `cat $TSTNAME.test`
do
    if [ -f $i ]
    then
        echo @01 $i ===========
	echo @01 $i ===========                >> $ERRFILE
	echo read $i                           >  $TMPFILE
	echo set load $SETFILE                 >> $TMPFILE
	echo set limits nodelimit 10000000     >> $TMPFILE
	echo set limits timelimit $TIMELIMIT   >> $TMPFILE
	echo set timing clocktype 1            >> $TMPFILE
	echo set display verblevel 4           >> $TMPFILE
	echo set display dispfreq 10000        >> $TMPFILE
	echo optimize                          >> $TMPFILE
	echo display statistics                >> $TMPFILE
	echo quit                              >> $TMPFILE
	../$2 < $TMPFILE 2>>$ERRFILE
	echo =ready=
    else
	echo @02 FILE NOT FOUND: $i ===========
	echo @02 FILE NOT FOUND: $i =========== >>$ERRFILE
    fi
done | tee -a $OUTFILE

rm $TMPFILE

date >>$OUTFILE
date >>$ERRFILE

gawk -f check.awk -vTEXFILE=$TEXFILE $TSTNAME.solu $OUTFILE | tee $RESFILE
