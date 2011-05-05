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


