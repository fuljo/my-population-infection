#!/bin/bash

#------------------------------------------------------------------------------
# Measure execution time with increasing number of countries
#------------------------------------------------------------------------------

program=$0;
exec="../src/my-population-infection";
worldsize=4096;
individuals=30000;

function usage {
    echo "Usage: $program min increment max";
}

function fmtduration() {
    echo -n "$(($1/3600)):$(($1%3600/60)):$(($1%60))";
}

if [ -z $1 ] || [ -z $2 ] || [ -z $3 ] ; then
    usage;
    exit 1;
fi

if [ $1 -gt $3 ]; then
    echo "ncountries_min must be smaller or equal to ncountries_max";
    exit 1;
fi

outfile="profile_countries_$1_$3.csv";

# Create directory for results
mkdir -p ./results;

# Write header of csv file
echo "rows,cols,countries,individuals,user,system,elapsed,memory" > $outfile;

# Iterate over countries
for countries in $(seq $1 $2 $3)
do
    echo "Simulating with $countries processes... ";
    # Save start time
    start=$(date +%s);
    # Perform the simulation
    /usr/bin/time -o $outfile -a -f "1,$countries,$countries,$individuals,%U,%S,%e,%M" \
        mpirun -np $countries --oversubscribe $exec \
            -N $(($individuals)) \
            -I $(($individuals / 100)) \
            -W $(($worldsize - $worldsize % $countries)) \
            -L $worldsize \
            -w $(($worldsize / $countries)) \
            -l $worldsize \
            -v 1.4 \
            -d 2 \
            --t-infection=$((1 * 60)) \
            --t-recovery=$((1 * 3600)) \
            --t-immunity=$((4 * 3600)) \
            --sim-step=60 \
            --sim-length=1 \
            --log-level WARN \
            ;
    # Save end time and calculate duration
    end=$(date +%s);
    duration=$(($end - $start));

    # Write to console and csv
    echo "Duration $duration seconds ( $(fmtduration $duration) )";
done

# Clean up results of simulations
rm -rf ./results;