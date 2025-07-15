#!/bin/bash

#Get the current working directory
WORK_DIR=$(pwd)
echo "Info: Working directory: $WORK_DIR"

#Change to the perf_runs directory
cd perf_runs 
if [ $? -ne 0 ]; then
    echo "Error: Failed to change directory to perf_runs."
    exit 1
fi

NUM_THREADS="-num_threads="
SUFF=""
EXE="./blastp"

#Check if the script is called with multithreaded option
if [ "$1" == "--mt" ]; then
    echo "Info: Running performance tests in multithreaded mode."
    NUM_THREADS+="16"
    SUFF="_mt"
elif [ "$1" == "--st" ]; then
    echo "Info: Running performance tests in single-threaded mode."
    NUM_THREADS+="1"
    SUFF="_st"
elif [ "$1" == "--blaze" ]; then
    echo "Info: Running BLAZE performance tests."
    NUM_THREADS+="16"
    SUFF="_blaze"
    EXE="./blaze_blastp"
else
    echo "Error: Invalid option. Use --mt for multithreaded or --st for single-threaded."
    exit 1
fi

#Check if the script is called with database option with the path
if [ -z "$2" ]; then
    echo "Error: No database path provided. Please provide the path to the database."
    exit 1
else
    DB_PATH="$2"
    # Ensure that the database file exists with a wildcard check as it might have an extension
    shopt -s nullglob
    db_files=(${DB_PATH}.*)
    if [ ${#db_files[@]} -eq 0 ]; then
        echo "Error: Database file $DB_PATH does not exist."
        exit 1
    fi
    echo "Info: Using database file: $DB_PATH"
    shopt -u nullglob
fi

time (  $EXE $NUM_THREADS -query NP_670887.fa -db $DB_PATH 1>NP_670887_op$SUFF ) 2> 1_time$SUFF
time (  $EXE $NUM_THREADS -query YP_244862.fa -db $DB_PATH 1>YP_244862_op$SUFF ) 2> 2_time$SUFF
time (  $EXE $NUM_THREADS -query NP_308413.fa -db $DB_PATH 1>NP_308413_op$SUFF ) 2> 3_time$SUFF
time (  $EXE $NUM_THREADS -query NP_931067.fa -db $DB_PATH 1>NP_931067_op$SUFF ) 2> 4_time$SUFF
time (  $EXE $NUM_THREADS -query YP_182226.fa -db $DB_PATH 1>YP_182226_op$SUFF ) 2> 5_time$SUFF
time (  $EXE $NUM_THREADS -query YP_054627.fa -db $DB_PATH 1>YP_054627_op$SUFF ) 2> 6_time$SUFF
time (  $EXE $NUM_THREADS -query YP_209371.fa -db $DB_PATH 1>YP_209371_op$SUFF ) 2> 7_time$SUFF
time (  $EXE $NUM_THREADS -query YP_141079.fa -db $DB_PATH 1>YP_141079_op$SUFF ) 2> 8_time$SUFF
time (  $EXE $NUM_THREADS -query NP_436070.fa -db $DB_PATH 1>NP_436070_op$SUFF ) 2> 9_time$SUFF
time (  $EXE $NUM_THREADS -query YP_000338.fa -db $DB_PATH 1>YP_000338_op$SUFF ) 2> 10_time$SUFF
time (  $EXE $NUM_THREADS -query NP_349135.fa -db $DB_PATH 1>NP_349135_op$SUFF ) 2> 11_time$SUFF
time (  $EXE $NUM_THREADS -query NP_566269.fa -db $DB_PATH 1>NP_566269_op$SUFF ) 2> 12_time$SUFF
time (  $EXE $NUM_THREADS -query NP_894055.fa -db $DB_PATH 1>NP_894055_op$SUFF ) 2> 13_time$SUFF
time (  $EXE $NUM_THREADS -query YP_137100.fa -db $DB_PATH 1>YP_137100_op$SUFF ) 2> 14_time$SUFF
time (  $EXE $NUM_THREADS -query YP_062020.fa -db $DB_PATH 1>YP_062020_op$SUFF ) 2> 15_time$SUFF
time (  $EXE $NUM_THREADS -query NP_864591.fa -db $DB_PATH 1>NP_864591_op$SUFF ) 2> 16_time$SUFF
time (  $EXE $NUM_THREADS -query NP_214548.fa -db $DB_PATH 1>NP_214548_op$SUFF ) 2> 17_time$SUFF
time (  $EXE $NUM_THREADS -query YP_069013.fa -db $DB_PATH 1>YP_069013_op$SUFF ) 2> 18_time$SUFF
time (  $EXE $NUM_THREADS -query NP_872878.fa -db $DB_PATH 1>NP_872878_op$SUFF ) 2> 19_time$SUFF
time (  $EXE $NUM_THREADS -query NP_345015.fa -db $DB_PATH 1>NP_345015_op$SUFF ) 2> 20_time$SUFF
time (  $EXE $NUM_THREADS -query NP_085504.fa -db $DB_PATH 1>NP_085504_op$SUFF ) 2> 21_time$SUFF
time (  $EXE $NUM_THREADS -query NP_174344.fa -db $DB_PATH 1>NP_174344_op$SUFF ) 2> 22_time$SUFF
time (  $EXE $NUM_THREADS -query NP_611552.fa -db $DB_PATH 1>NP_611552_op$SUFF ) 2> 23_time$SUFF
time (  $EXE $NUM_THREADS -query XP_610726.fa -db $DB_PATH 1>XP_610726_op$SUFF ) 2> 24_time$SUFF
time (  $EXE $NUM_THREADS -query YP_189837.fa -db $DB_PATH 1>YP_189837_op$SUFF ) 2> 25_time$SUFF
time (  $EXE $NUM_THREADS -query NP_245243.fa -db $DB_PATH 1>NP_245243_op$SUFF ) 2> 26_time$SUFF
time (  $EXE $NUM_THREADS -query NP_818972.fa -db $DB_PATH 1>NP_818972_op$SUFF ) 2> 27_time$SUFF
time (  $EXE $NUM_THREADS -query NP_633288.fa -db $DB_PATH 1>NP_633288_op$SUFF ) 2> 28_time$SUFF
time (  $EXE $NUM_THREADS -query NP_781168.fa -db $DB_PATH 1>NP_781168_op$SUFF ) 2> 29_time$SUFF
time (  $EXE $NUM_THREADS -query NP_784621.fa -db $DB_PATH 1>NP_784621_op$SUFF ) 2> 30_time$SUFF
time (  $EXE $NUM_THREADS -query NP_952492.fa -db $DB_PATH 1>NP_952492_op$SUFF ) 2> 31_time$SUFF
time (  $EXE $NUM_THREADS -query YP_121326.fa -db $DB_PATH 1>YP_121326_op$SUFF ) 2> 32_time$SUFF
time (  $EXE $NUM_THREADS -query NP_965272.fa -db $DB_PATH 1>NP_965272_op$SUFF ) 2> 33_time$SUFF
time (  $EXE $NUM_THREADS -query NP_706440.fa -db $DB_PATH 1>NP_706440_op$SUFF ) 2> 34_time$SUFF
time (  $EXE $NUM_THREADS -query XP_499151.fa -db $DB_PATH 1>XP_499151_op$SUFF ) 2> 35_time$SUFF
time (  $EXE $NUM_THREADS -query XP_237451.fa -db $DB_PATH 1>XP_237451_op$SUFF ) 2> 36_time$SUFF
time (  $EXE $NUM_THREADS -query NP_733286.fa -db $DB_PATH 1>NP_733286_op$SUFF ) 2> 37_time$SUFF
time (  $EXE $NUM_THREADS -query YP_129633.fa -db $DB_PATH 1>YP_129633_op$SUFF ) 2> 38_time$SUFF
time (  $EXE $NUM_THREADS -query YP_081219.fa -db $DB_PATH 1>YP_081219_op$SUFF ) 2> 39_time$SUFF