#! /bin/bash
echo "organizing includes..."
for i in newfiles/scip/*.c
do
   echo $i
   ./list_includes.py $i > tmp.list
   echo "#include \"scip/$(basename $i .c).h\"" >> tmp.list
   sort tmp.list -o tmp.list

   sed -e '/ctype/{:a; N; /pub_expr/!ba; r tmp.list' -e 'd;}' $i  > tmpnewfile.c

   mv tmpnewfile.c $i
done
