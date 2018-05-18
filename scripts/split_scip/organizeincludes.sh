#! /bin/bash
echo "organizing includes..."
for i in newfiles/scip/*.c
do
   echo $i
   ./list_includes.py $i > tmp.list
   echo "#include \"scip/$(basename $i .c).h\"" >> tmp.list
   cat tmp.list | grep scip/scip_ | sort -u > tmp2.list
   mv tmp2.list tmp.list
   sed -e '/scip_bandit.h/{:a; N; /scip_var.h/!ba; r tmp.list' -e 'd;}' $i  > tmpnewfile.c

   mv tmpnewfile.c $i
done
