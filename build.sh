#!/bin/bash

#Function to reoplace a file with another one
# Throw an error if either file does not exist
replace_file() {
    if [ ! -f "$1" ]; then
        echo "Error: Source file $1 does not exist."
        exit 1
    fi
    if [ ! -f "$2" ]; then
        echo "Error: Destination file $2 does not exist."
        exit 1
    fi

    # Create a backup of the destination file
    cp "$2" "$2.bak"
    cp "$1" "$2"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to replace $2 with $1."
        exit 1
    fi
}

echo "----------------------------------------"
echo "Info: Starting build script for BLAZE"
echo "----------------------------------------"

# BLAZE has been tested with ncbi-blast-2.13.0+-src
NCBI_VERSION="ncbi-blast-2.13.0+-src"
NCBI_SRC="$NCBI_VERSION.zip"
# Get the current working directory
WORK_DIR=$(pwd)
echo "Info: Working directory: $WORK_DIR"
# Get path of NVCC, if not found print error and exit
NVCC=$(which nvcc)
if [ -z "$NVCC" ]; then
    echo "Error: NVCC not found. Please install CUDA toolkit."
    exit 1
fi
echo "Info: NVCC found at: $NVCC"

echo "Info: Done verifying pre-requisites."
echo "----------------------------------------"
echo "Info: Downloading NCBI BLAST src..."
echo "----------------------------------------"

# Check if the NCBI source file exists in the current directory
if [ ! -f "$NCBI_SRC" ]; then
    echo "Info: $NCBI_SRC not found in the current directory."
    echo "Info: Downloading NCBI BLAST source code..."
    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/ncbi-blast-2.13.0+-src.zip
    if [ $? -ne 0 ]; then
        echo "Error: Failed to download $NCBI_SRC."
        exit 1
    fi
    echo "Info: Download complete."
else
    echo "Info: $NCBI_SRC found in the current directory."
fi

# Unzip the NCBI source code
if [ ! -d "ncbi-blast-2.13.0+-src" ]; then
    echo "Info: Unzipping $NCBI_SRC..."
    unzip -q "$NCBI_SRC"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to unzip $NCBI_SRC."
        exit 1
    fi
    echo "Unzip complete."
else
    echo "Info: NCBI source directory already exists."
    echo "Info: This script will delete the existing directory and unzip again."
    # Take user confirmation
    read -p "User Input: Do you want to delete the existing directory and unzip again? (y/n): " confirm
    if [[ "$confirm" == "y" || "$confirm" == "Y" ]]; then
        rm -rf ncbi-blast-2.13.0+-src
        echo "Info: Deleted existing directory."
        unzip -q "$NCBI_SRC"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to unzip $NCBI_SRC."
            exit 1
        fi
        echo "Info: Unzip complete."
    else
        echo "Info: Exiting without changes."
        exit 0
    fi
fi    

echo "Info: Done obtaining NCBI BLAST source code."

echo "----------------------------------------"
echo "Info: Configuring NCBI BLAST build..."
echo "----------------------------------------"

# Change to the NCBI source directory
cd $WORK_DIR/$NCBI_VERSION/c++
echo "Info: Changed directory to $(pwd)"

# Run the configure script with the following options
srcdir=$WORK_DIR/$NCBI_VERSION/c++

echo "Info: Running config with --with-debug --without-boost --with-symbols --with-openmp --with-mt --without-vdb --without-gnutls --without-gcrypt --with-build-root=$srcdir/ReleaseMT --with-experimental=Int8GI "

# Debug print "If options need to be modified, please edit this script and re-run it."
echo "If options need to be modified, please edit this script and re-run it."

echo "Info: Running configure script with specified options..."
# Pipe output to a log file
$srcdir/configure.orig --with-debug --without-boost --with-symbols --with-openmp --with-mt --without-vdb --without-gnutls --without-gcrypt --with-build-root=$srcdir/ReleaseMT --with-experimental=Int8GI > configure.log 2>&1

if [ $? -ne 0 ]; then
    echo "Error: Configuration failed."
    exit 1
fi
echo "Info: Configuration complete."
echo "----------------------------------------"
echo "Info: Setting the DIRs..."
echo "----------------------------------------"

CPP_DIR="$srcdir"
RELEASE_DIR="$srcdir/ReleaseMT"
BIN_DIR="$RELEASE_DIR/bin"
BUILD_DIR="$RELEASE_DIR/build"

# Print the directories
echo "Info: CPP_DIR: $CPP_DIR"
echo "Info: RELEASE_DIR: $RELEASE_DIR"
echo "Info: BUILD_DIR: $BUILD_DIR"
echo "Info: BIN_DIR: $BIN_DIR"

echo "----------------------------------------"
echo "Info: Modifying ncbi-source files to include BLAZE"
echo "----------------------------------------"

echo "Info: Modifying makefiles..." 
replace_file "$WORK_DIR/BLAZE_mk/src_app_blast_Makefile.in" "$CPP_DIR/src/app/blast/Makefile.in"
replace_file "$WORK_DIR/BLAZE_mk/src_app_blastdb_Makefile.in" "$CPP_DIR/src/app/blastdb/Makefile.in"
replace_file "$WORK_DIR/BLAZE_mk/src_app_Makefile.in" "$CPP_DIR/src/app/Makefile.in"
replace_file "$WORK_DIR/BLAZE_mk/src_app_blast_Makefile.blastp.app" "$CPP_DIR/src/app/blast/Makefile.blastp.app"
replace_file "$WORK_DIR/BLAZE_mk/src_app_blastdb_Makefile.makeblastdb.app" "$CPP_DIR/src/app/blastdb/Makefile.makeblastdb.app"

# # Backup Makefile.mk 
# if [ -f "${BUILD_DIR}/Makefile.mk" ]; then
#     cp "${BUILD_DIR}/Makefile.mk" "${BUILD_DIR}/Makefile.mk.bak"
# else
#     echo "Error: ${BUILD_DIR}/Makefile.mk not found."
#     exit 1
# fi

# # Safely append -lblaze -lcudart to CONF_LIBS line
# awk '/^CONF_LIBS/ { $0 = $0 " -lblaze -lcudart" } { print }' "${BUILD_DIR}/Makefile.mk" > "${BUILD_DIR}/Makefile.mk.temp"
# if [ $? -ne 0 ]; then
#     echo "Error: Failed to update CONF_LIBS in Makefile.mk."
#     exit 1
# fi
# mv "${BUILD_DIR}/Makefile.mk.temp" "${BUILD_DIR}/Makefile.mk"

echo "Info: Copying BLAZE source files to NCBI source directory..."

replace_file "$WORK_DIR/BLAZE_ncbi/setup_factory.hpp" "$CPP_DIR/include/algo/blast/api/setup_factory.hpp"
replace_file "$WORK_DIR/BLAZE_ncbi/blast_engine.h" "$CPP_DIR/include/algo/blast/core/blast_engine.h"
replace_file "$WORK_DIR/BLAZE_ncbi/blast_seqsrc.h" "$CPP_DIR/include/algo/blast/core/blast_seqsrc.h"
replace_file "$WORK_DIR/BLAZE_ncbi/prelim_search_runner.hpp" "$CPP_DIR/src/algo/blast/api/prelim_search_runner.hpp"
replace_file "$WORK_DIR/BLAZE_ncbi/prelim_stage.cpp" "$CPP_DIR/src/algo/blast/api/prelim_stage.cpp"
replace_file "$WORK_DIR/BLAZE_ncbi/blast_engine.c" "$CPP_DIR/src/algo/blast/core/blast_engine.c"
replace_file "$WORK_DIR/BLAZE_ncbi/blast_seqsrc.c" "$CPP_DIR/src/algo/blast/core/blast_seqsrc.c"
replace_file "$WORK_DIR/BLAZE_ncbi/seqdbimpl.cpp" "$CPP_DIR/src/objtools/blast/seqdb_reader/seqdbimpl.cpp"

# Copy the header files to the include directory
cp "$WORK_DIR/BLAZE_gpu/gpu_cpu_common.h" "$CPP_DIR/include/algo/blast/core/"

# Make a gpu directory in the src/algo/blast directory
mkdir -p $CPP_DIR/src/algo/blast/gpu/

#Copy the source files to the src directory
cp "$WORK_DIR/BLAZE_gpu/gpuSeedCount.cu" "$CPP_DIR/src/algo/blast/gpu/gpuSeedCount.cu"
cp "$WORK_DIR/BLAZE_gpu/util.hpp" "$CPP_DIR/src/algo/blast/gpu/util.hpp"

echo "----------------------------------------"
echo "Info: Copying libcudart.so..."
echo "----------------------------------------"

# Extract the CUDA directroy from NVCC path
CUDA_DIR=$(dirname $(dirname $NVCC))
echo "Info: CUDA_DIR: $CUDA_DIR"
# Attempting to find libcudart.so depending on the 32/64 bit architecture
if [ "$(uname -m)" == "x86_64" ]; then
    CUDART_LIB="$CUDA_DIR/lib64/libcudart.so"
else
    CUDART_LIB="$CUDA_DIR/lib/libcudart.so"
fi

# Check if libcudart.so exists
if [ ! -f "$CUDART_LIB" ]; then
    echo "Error: libcudart.so not found in $CUDART_LIB. Please check your CUDA installation."
    exit 1
fi  

echo "Info: Copying libcudart.so to $RELEASE_DIR/lib/"
cp "$CUDART_LIB" "$RELEASE_DIR/lib/"

# # Check if the copy was successful
# if [ $? -ne 0 ]; then
#     echo "Error: Failed to copy libcudart.so to $RELEASE_DIR/lib/"
#     exit 1
# fi

echo "----------------------------------------"
echo "Info: Compiling BLAZE into a static library..."
echo "----------------------------------------"

CPP_DIR="$WORK_DIR/$NCBI_VERSION/c++"

# Check the source file exists
if [ ! -f "$CPP_DIR/src/algo/blast/gpu/gpuSeedCount.cu" ]; then
    echo "Error: Source file gpuSeedCount.cu not found at $CPP_DIR/src/algo/blast/gpu/gpuSeedCount.cu"
    exit 1
fi

echo "Info: Compiling for CUDA architecture compute_86, code sm_86"
$NVCC --extra-device-vectorization --use_fast_math -O3 -c -gencode arch=compute_86,code=sm_86 --compiler-options -fPIC --compiler-options -fno-strict-aliasing --compiler-options -fno-inline -DUNIX -O3 --no-align-double -L /lib64  -lusb-1.0 -l pthread -shared  $CPP_DIR/src/algo/blast/gpu/gpuSeedCount.cu -o blazeog.o -I$CPP_DIR/include/ -I$CPP_DIR/ReleaseMT/inc/ > compile_gpuSeedCount.log 2>&1

if [ $? -ne 0 ]; then
  echo -e "\033[31m [!!ERROR!!] Compile failed, check CUDA code\033[0m"
  exit 1
else
    echo -e "\033[32m [SUCCESS] CUDA code compiled!\033[0m"
fi

# Copy the compiled object files to the release directory
cp blazeog.o "$RELEASE_DIR/lib/libblaze.so"

echo "----------------------------------------"
echo "Info: Compiling NCBI BLAST + BLAZE this will take a while..."
echo "----------------------------------------"

cd "$CPP_DIR"
make -j$(nproc) > make_blast.log 2> make_blast.err &
MAKE_PID=$!

# Spinner or periodic message
spin='|/-\\'
i=0
while kill -0 $MAKE_PID 2>/dev/null; do
    i=$(( (i+1) %4 ))
    printf "\r[Info] Compiling... ${spin:$i:1}"
    sleep 1
done
wait $MAKE_PID
printf "\r[Info] Compilation finished.      \n"

if [ $? -ne 0 ]; then
    echo "Error: NCBI BLAST compilation failed. Check make_blast.err for details."
    exit 1
fi

echo "----------------------------------------"
echo "Info: Done compiling NCBI BLAST + BLAZE."
echo "BINs: $BIN_DIR"
echo "Info: You can now run NCBI BLAST with BLAZE support."
echo "----------------------------------------"