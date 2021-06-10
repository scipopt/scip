#!/usr/bin/awk
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*                                                                           *
#*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
#
BEGIN {
    printf("FILE                    EMPTY     FREE     SING     AGGR    VARBD  SETPART  SETPACK   SETCOV     CARD  INVKNAP   EQKNAP  BINPACK     KNAP  INTKNAP   MIXBIN      GEN\n");
    out = 0;
}

// {
    if( out == 1 )
    {
	printf("%20s", shortprob);
	print;
	shortprob = "";
	out = 0;
    }
}

/read / {
    if( shortprob != "" )
    {
	printf("%20s  ABORTED\n", shortprob);
    }

   filename = $3;

   n  = split (filename, a, "/");
   m = split(a[n], b, ".");
   prob = b[1];
   if( b[m] == "gz" || b[m] == "z" || b[m] == "GZ" || b[m] == "Z" )
      m--;
   for( i = 2; i < m; ++i )
      prob = prob "." b[i];

   if( length(prob) > 18 )
      shortprob = substr(prob, length(prob)-17, 18);
   else
      shortprob = prob;
}

/MIXBIN/ {
    out = 1;
}


