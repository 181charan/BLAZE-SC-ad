#!/bin/bash

# Get the current working directory
WORK_DIR=$(pwd)
echo "Info: Working directory: $WORK_DIR"

NCBI_ARCHIVE="ncbi-blast-2.13.0+-x64-linux.tar.gz"

# Check if the NCBI archive exists in the current directory, if so remove it
if [ -f "$NCBI_ARCHIVE" ]; then
    echo "Info: $NCBI_ARCHIVE found in the current directory."
    echo "Info: Removing $NCBI_ARCHIVE..."
    rm "$NCBI_ARCHIVE"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to remove $NCBI_ARCHIVE."
        exit 1
    fi
    echo "Info: $NCBI_ARCHIVE removed successfully."
fi

# Download the NCBI BLAST archive
echo "Info: Downloading NCBI BLAST archive..."
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/ncbi-blast-2.13.0+-x64-linux.tar.gz
if [ $? -ne 0 ]; then
    echo "Error: Failed to download $NCBI_ARCHIVE."
    exit 1
fi
echo "Info: Download complete."

# Extract the NCBI BLAST archive
echo "Info: Extracting $NCBI_ARCHIVE..."
tar -xzf "$NCBI_ARCHIVE"
if [ $? -ne 0 ]; then
    echo "Error: Failed to extract $NCBI_ARCHIVE."
    exit 1
fi
echo "Info: Extraction complete."

# Remove if perf_runs directory exists
if [ -d "perf_runs" ]; then
    echo "Info: Removing existing perf_runs directory..."
    rm -rf perf_runs
    if [ $? -ne 0 ]; then
        echo "Error: Failed to remove perf_runs directory."
        exit 1
    fi
    echo "Info: perf_runs directory removed successfully."
fi

# Create the perf_runs directory
echo "Info: Creating perf_runs directory..."
mkdir perf_runs
if [ $? -ne 0 ]; then
    echo "Error: Failed to create perf_runs directory."
    exit 1
fi

# Move the NCBI BLAST executables to the perf_runs directory
echo "Info: Moving NCBI BLAST executables to perf_runs directory..."
mv ncbi-blast-2.13.0+/bin/blastp perf_runs/
mv ncbi-blast-2.13.0+/bin/makeblastdb perf_runs/

if [ $? -ne 0 ]; then
    echo "Error: Failed to move NCBI BLAST executables."
    exit 1
fi
echo "Info: NCBI BLAST executables moved successfully."

echo "----------------------------------------"
echo "Info: NCBI BLAST setup complete."
echo "----------------------------------------"

# This script assumes that the compiled BLAZE executable has already been compiled successfully

BLAZE_BIN="$WORK_DIR/ncbi-blast-2.13.0+-src/c++/ReleaseMT/bin/blastp"
if [ ! -f "$BLAZE_BIN" ]; then
    echo "Error: BLAZE executable not found at $BLAZE_BIN."
    echo "Please ensure that the BLAZE executable is compiled and available at the specified path."
    exit 1
fi

echo "Info: BLAZE executable found at $BLAZE_BIN."

#Copy the BLAZE executable to the perf_runs directory
cp "$BLAZE_BIN" perf_runs/blaze_blastp
if [ $? -ne 0 ]; then
    echo "Error: Failed to copy BLAZE executable to perf_runs directory."
    exit 1
fi
echo "Info: BLAZE executable copied to perf_runs directory."

# Copy queries into the perf_runs directory
echo "Info: Copying queries to perf_runs directory..."
cp "$WORK_DIR/BLAZE_queries/"*.fa perf_runs/

if [ $? -ne 0 ]; then
    echo "Error: Failed to copy queries to perf_runs directory."
    exit 1
fi
echo "Info: Queries copied successfully."

echo "----------------------------------------"
echo "Info: Running single threaded tests."
echo "----------------------------------------"

# Take user input for the database path 
echo "Info: Please provide the path to the database file (e.g., blast_db):"
read -r DB_PATH
# Ensure that the database file exists with a wildcard check as it might have an extension
shopt -s nullglob
db_files=(${DB_PATH}.*)
if [ ${#db_files[@]} -eq 0 ]; then
    echo "Error: Database file $DB_PATH does not exist."
    exit 1
fi
echo "Info: Using database file: $DB_PATH"
shopt -u nullglob

bash blaze_runner.sh --st "$DB_PATH"

echo "----------------------------------------"
echo "Info: Running multithreaded tests."
echo "----------------------------------------"
bash blaze_runner.sh --mt "$DB_PATH"

echo "----------------------------------------"
echo "Info: Running BLAZE performance tests."
echo "----------------------------------------"
bash blaze_runner.sh --blaze "$DB_PATH" 
echo "----------------------------------------"
echo "Info: Performance tests completed."