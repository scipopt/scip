# $Id: check.sh,v 1.2 2003/11/20 17:42:00 bzfpfend Exp $
BINNAME=`basename $2`
TSTNAME=`basename $1 .test`
SETNAME=$3

OUTFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.out
ERRFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.err
RESFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.res
TEXFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.tex
TMPFILE=results/check.$TSTNAME.$BINNAME.$SETNAME.tmp
SETFILE=settings/$SETNAME.set

date >$OUTFILE
date >$ERRFILE

for i in `cat $1`
do
    if [ -f $i ]
    then
        echo @01 $i ===========
	echo @01 $i ===========                >> $ERRFILE
	echo read $i                           >  $TMPFILE
	echo set load $SETFILE                 >> $TMPFILE
	echo set limits nodelimit 10000000     >> $TMPFILE
	echo set limits timelimit 3600         >> $TMPFILE
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
 
