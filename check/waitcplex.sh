#!/bin/sh
while true
do
    AVAILABLE=`/home/optimi/adm_cple/ilm-2.6/ilmlist | gawk --source '
BEGIN { incplex = 0; available = -1; }
/^\* PRODUCT/ {
    if( $3 == "CPLEX:" )
	incplex = 1;
    else
	incplex = 0;
    }
/^available tokens:/ {
    if( incplex )
	available = $3;
}
END {
    printf("%d", available);
}'`

    echo available CPLEX tokens: $AVAILABLE

    if [ $AVAILABLE -ge 2 ]
    then
	break
    fi

    sleep 30
done
