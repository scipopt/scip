# $Id: check.sh,v 1.6 2004/01/22 14:42:25 bzfpfend Exp $
TSTNAME=$1
BINNAME=$2
SETNAME=$3
BINID=$4
TIMELIMIT=$5
NODELIMIT=$6
MEMLIMIT=$7

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
	echo set load $SETFILE                 >  $TMPFILE
	echo set limits timelimit $TIMELIMIT   >> $TMPFILE
	echo set limits nodelimit $NODELIMIT   >> $TMPFILE
	echo set limits memlimit $MEMLIMIT     >> $TMPFILE
	echo set timing clocktype 1            >> $TMPFILE
	echo set display verblevel 3           >> $TMPFILE
	echo set display dispfreq 10000        >> $TMPFILE
	echo read $i                           >> $TMPFILE
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
