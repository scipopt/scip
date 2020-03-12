BEGIN {
    inreport = 0;
    intreeestimation = 0;
}
/^Report/{
    inreport = 1;
    reportnum = $2;
}

# Time Elapsed: 17.68
/^Time Elapsed/ {
    timeelapsed = $3;
}

#   Tree Data        : 17 nodes (10 visited, 8 inner, 2 leaves, 7 open), weight: 0.0117

/^  Tree Data|^Estimation Tree/ {
    report["visited"] = substr($6, 2);
    report["inner"] = $8
    report["leaves"] = $10;
    report["open"] = $12;
    weight = $15
}
#Estimation Tree    : 4859 nodes (4852 visited, 2429 inner, 2423 leaves, 7 open), weight: 0.9909 completed 0.9461


/^Tree Estimation|^Estimations/ {
    intreeestimation = 1;
}

/^  / {
    if( intreeestimation )
    {
        name =  substr($0, 3, index($0, ":") - 3)
        do{
            gsub(" +", "", name);
        } while (index(name, " ") >= 1)

        if( $2 == ":" )
        {
            estim[name] = $3
            value[name] = $4
            trend[name] = $5
            resolution[name] = $6
            smooth[name] = $7
        }
        else
        {
            estim[name] = $4
            value[name] = $5
            trend[name] = $6
            resolution[name] = $7
            smooth[name] = $8
        }

    }
}

# Tree Estimation    :       estim       value       trend  resolution
#   wbe              :         255           -           -           -
#   tree profile     :          -1           -           -           -
#   gap              :          20     0.00000     0.00000           1
#   tree-weight      :         270     0.01172     0.00740           1
#   leaf-frequency   :           4     0.15000     0.51706           1
#   ssg              :         403     1.00000    -0.00500           1
# End of Report 1

/End of Report/ {
    for( key in estim )
    {
        printf "%s,", FILENAME;
        printf "%s,", key;
        printf "%d,%.2f,", reportnum, timeelapsed;
        printf "%d,%d,%d,", report["visited"], report["leaves"], report["open"];
        printf "%f,", weight;
        printf "%d,", estim[key];
        if( value[key] != "-" )
            printf "%f,", value[key];
        else
            printf "NA,"
        if( trend[key] != "-" )
            printf "%f,", trend[key];
        else
            printf "NA,"
        if( resolution[key] != "-" )
            printf "%d,", resolution[key];
        else
            printf "NA,"
        if( smooth[key] != "-" )
            printf "%d", smooth[key];
        else
            printf "NA"
        printf "\n";

    }
    intreeestimation = 0
}