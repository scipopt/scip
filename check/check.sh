# $Id: check.sh,v 1.1 2002/10/23 14:31:35 bzfpfend Exp $
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
    echo @01 $i ===========
    echo @01 $i =========== >>$ERRFILE
    if [ -f $i ]
    then
	cp $SETFILE $TMPFILE
	echo read $i                    >> $TMPFILE
	echo set limit bab 10000000     >> $TMPFILE
	echo set limit time 3600        >> $TMPFILE
	echo set time cpu               >> $TMPFILE
	echo set info level 2           >> $TMPFILE
	echo set info frequ 10000       >> $TMPFILE
	echo optimize                   >> $TMPFILE
	echo display stat               >> $TMPFILE
	echo quit                       >> $TMPFILE
	../$2 < $TMPFILE 2>>$ERRFILE
    else
	echo FILE NOT FOUND
	echo FILE NOT FOUND >>$ERRFILE
    fi
    echo =ready=
done | tee -a $OUTFILE

rm $TMPFILE

date >>$OUTFILE
date >>$ERRFILE

gawk -f check.awk -vTEXFILE=$TEXFILE $TSTNAME.solu $OUTFILE | tee $RESFILE
 
