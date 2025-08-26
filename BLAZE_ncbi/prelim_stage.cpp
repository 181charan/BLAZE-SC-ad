/* ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Author:  Christiam Camacho
 *
 */

/** @file prelim_stage.cpp
 * NOTE: This file contains work in progress and the APIs are likely to change,
 * please do not rely on them until this notice is removed.
 */

#include <ncbi_pch.hpp>

#include <algo/blast/api/prelim_stage.hpp>
#include <algo/blast/api/uniform_search.hpp>    // for CSearchDatabase
#include <algo/blast/api/blast_mtlock.hpp>
#include <algo/blast/core/blast_hits.h>
#include <algo/blast/core/blast_stat.h>

#include "prelim_search_runner.hpp"
#include "blast_aux_priv.hpp"
#include "psiblast_aux_priv.hpp"
#include "split_query_aux_priv.hpp"
#include "blast_seqalign.hpp"
#include <sstream>

#include <algo/blast/api/blast_dbindex.hpp>
#include <algo/blast/core/gpu_cpu_common.h>
// #include <algo/blast/core/cpu_util.hpp>
#include <time.h>
#include <cstdio> // Add this to use c style file I/O
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <sys/mman.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h> // Add this line
#include <cstring>

int score_matrix[28][28] = {{-16384,-16384,-16384,-16384,-16384,-16384,-16384,-16384,-16384,-16384,-16384,-16384,-16384,-16384,-16384,-16384,-16384,-16384,-16384,-16384,-16384,-16384,-16384,-16384,-16384,-16384,-16384,-16384},
{-16384,4,-2,0,-2,-1,-2,0,-2,-1,-1,-1,-1,-2,-1,-1,-1,1,0,0,-3,-1,-2,-1,0,-4,-1,-1},
{-16384,-2,4,-3,4,1,-3,-1,0,-3,0,-4,-3,4,-2,0,-1,0,-1,-3,-4,-1,-3,0,-3,-4,-1,-3},
{-16384,0,-3,9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-1,-2,-3,9,-4,-1,-1},
{-16384,-2,4,-3,6,2,-3,-1,-1,-3,-1,-4,-3,1,-1,0,-2,0,-1,-3,-4,-1,-3,1,-3,-4,-1,-3},
{-16384,-1,1,-4,2,5,-3,-2,0,-3,1,-3,-2,0,-1,2,0,0,-1,-2,-3,-1,-2,4,-4,-4,-1,-3},
{-16384,-2,-3,-2,-3,-3,6,-3,-1,0,-3,0,0,-3,-4,-3,-3,-2,-2,-1,1,-1,3,-3,-2,-4,-1,0},
{-16384,0,-1,-3,-1,-2,-3,6,-2,-4,-2,-4,-3,0,-2,-2,-2,0,-2,-3,-2,-1,-3,-2,-3,-4,-1,-4},
{-16384,-2,0,-3,-1,0,-1,-2,8,-3,-1,-3,-2,1,-2,0,0,-1,-2,-3,-2,-1,2,0,-3,-4,-1,-3},
{-16384,-1,-3,-1,-3,-3,0,-4,-3,4,-3,2,1,-3,-3,-3,-3,-2,-1,3,-3,-1,-1,-3,-1,-4,-1,3},
{-16384,-1,0,-3,-1,1,-3,-2,-1,-3,5,-2,-1,0,-1,1,2,0,-1,-2,-3,-1,-2,1,-3,-4,-1,-3},
{-16384,-1,-4,-1,-4,-3,0,-4,-3,2,-2,4,2,-3,-3,-2,-2,-2,-1,1,-2,-1,-1,-3,-1,-4,-1,3},
{-16384,-1,-3,-1,-3,-2,0,-3,-2,1,-1,2,5,-2,-2,0,-1,-1,-1,1,-1,-1,-1,-1,-1,-4,-1,2},
{-16384,-2,4,-3,1,0,-3,0,1,-3,0,-3,-2,6,-2,0,0,1,0,-3,-4,-1,-2,0,-3,-4,-1,-3},
{-16384,-1,-2,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-2,7,-1,-2,-1,-1,-2,-4,-1,-3,-1,-3,-4,-1,-3},
{-16384,-1,0,-3,0,2,-3,-2,0,-3,1,-2,0,0,-1,5,1,0,-1,-2,-2,-1,-1,4,-3,-4,-1,-2},
{-16384,-1,-1,-3,-2,0,-3,-2,0,-3,2,-2,-1,0,-2,1,5,-1,-1,-3,-3,-1,-2,0,-3,-4,-1,-2},
{-16384,1,0,-1,0,0,-2,0,-1,-2,0,-2,-1,1,-1,0,-1,4,1,-2,-3,-1,-2,0,-1,-4,-1,-2},
{-16384,0,-1,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1,0,-1,-1,-1,1,5,0,-2,-1,-2,-1,-1,-4,-1,-1},
{-16384,0,-3,-1,-3,-2,-1,-3,-3,3,-2,1,1,-3,-2,-2,-3,-2,0,4,-3,-1,-1,-2,-1,-4,-1,2},
{-16384,-3,-4,-2,-4,-3,1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-2,-3,11,-1,2,-2,-2,-4,-1,-2},
{-16384,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-4,-1,-1},
{-16384,-2,-3,-2,-3,-2,3,-3,2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1,2,-1,7,-2,-2,-4,-1,-1},
{-16384,-1,0,-3,1,4,-3,-2,0,-3,1,-3,-1,0,-1,4,0,0,-1,-2,-2,-1,-2,4,-3,-4,-1,-3},
{-16384,0,-3,9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-1,-2,-3,9,-4,-1,-1},
{-16384,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,1,-4,-4},
{-16384,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-4,-1,-1},
{-16384,-1,-3,-1,-3,-3,0,-4,-3,3,-3,3,2,-3,-3,-2,-2,-2,-1,2,-2,-1,-1,-3,-1,-4,-1,3}};

std::vector<uint64_t> bbins = {
        0, 64, 128, 192, 256, 320, 384, 448, 
        512, 576, 1024, 2048, 1000000
        };  // 12 bins

// Returns the bin index (0-based) using a simple binary search approach
inline size_t findBinIndex(uint64_t value, const std::vector<uint64_t>& bins)
{
    // bins has length N; we want to find i such that bins[i] <= value < bins[i+1]
    // The last bin boundary is bins[N-1], so the maximum index is N-2
    // because we have (N-1) actual bins.
    size_t left = 0;
    size_t right = bins.size() - 2; 
    // Example: bins = [0, 64, 128, ..., 1000000]
    // We only search up to bins.size() - 2 because bins.size() - 1 is the last boundary.

    while (left <= right) {
        size_t mid = (left + right) / 2;
        if (value < bins[mid]) {
            right = (mid == 0) ? 0 : mid - 1;
        } 
        else if (value >= bins[mid + 1]) {
            left = mid + 1;
        }
        else {
            return mid;
        }
    }
    return left;
}


inline uint indx_calc(uint a, uint b, uint c){
    //return(((a << 5*2) | (b << 5*1) | (c << 5*0)));
    return ((a * pow(28,2)) + (b * pow(28,1)) + (c*(pow(28,0))));
}

std::vector<uint> neighbouring_words(uint x, uint y, uint z){

    std::vector<uint> neighbors;
    int score = 0;

    for(int a=0;a<28;a++){
        for(int b=0;b<28;b++){
            for(int c=0;c<28;c++){
                score = int(score_matrix[x][a] + score_matrix[y][b] + score_matrix[z][c]);
                if ( score >= 11  ){
                    if (x == a and y == b and z == c){
                        continue;
                    }
                    neighbors.push_back(indx_calc(a,b,c));
                }
            }
        }
    }

    return neighbors;
}


void write_neighbors_to_file(const std::vector<std::vector<uint>>& neighbors, const std::string& filename) {
    std::ofstream outFile(filename, std::ios::binary);
    if (!outFile) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }

    for (const auto& vec : neighbors) {
        size_t size = vec.size();
        outFile.write(reinterpret_cast<const char*>(&size), sizeof(size));
        outFile.write(reinterpret_cast<const char*>(vec.data()), size * sizeof(uint));
    }

    outFile.close();
}

void read_neighbors_from_file(std::vector<std::vector<uint>>& neighbors, const std::string& filename) {
    std::ifstream inFile(filename, std::ios::binary);
    if (!inFile) {
        std::cerr << "Error opening file for reading: " << filename << std::endl;
        return;
    }

    for (auto& vec : neighbors) {
        size_t size;
        inFile.read(reinterpret_cast<char*>(&size), sizeof(size));
        vec.resize(size);
        inFile.read(reinterpret_cast<char*>(vec.data()), size * sizeof(uint));
    }

    inFile.close();
}


/** @addtogroup AlgoBlast
 *
 * @{
 */

BEGIN_NCBI_SCOPE
USING_SCOPE(objects);
BEGIN_SCOPE(blast)

CBlastPrelimSearch::CBlastPrelimSearch(CRef<IQueryFactory> query_factory,
                                       CRef<CBlastOptions> options,
                                       const CSearchDatabase& dbinfo)
    : m_QueryFactory(query_factory), m_InternalData(new SInternalData),
    m_Options(options), m_DbAdapter(NULL), m_DbInfo(&dbinfo)
{
    BlastSeqSrc* seqsrc = CSetupFactory::CreateBlastSeqSrc(dbinfo);
    CRef<TBlastSeqSrc> wrapped_src(new TBlastSeqSrc(seqsrc, BlastSeqSrcFree));
    x_Init(query_factory, options, CRef<CPssmWithParameters>(), seqsrc);

    m_InternalData->m_SeqSrc = wrapped_src;
}

CBlastPrelimSearch::CBlastPrelimSearch(CRef<IQueryFactory> query_factory,
                                       CRef<CBlastOptions> options,
                                       CRef<CLocalDbAdapter> db,
                                       size_t num_threads)
    : m_QueryFactory(query_factory), m_InternalData(new SInternalData),
    m_Options(options), m_DbAdapter(db), m_DbInfo(NULL)
{
    BlastSeqSrc* seqsrc = db->MakeSeqSrc();
    x_Init(query_factory, options, CRef<CPssmWithParameters>(), seqsrc,
           num_threads);
    m_InternalData->m_SeqSrc.Reset(new TBlastSeqSrc(seqsrc, 0));
    if (num_threads > 1) {
        SetNumberOfThreads(num_threads);
    }
}

CBlastPrelimSearch::CBlastPrelimSearch(CRef<IQueryFactory> query_factory,
                               CRef<CBlastOptions> options,
                               BlastSeqSrc* seqsrc,
                               CConstRef<objects::CPssmWithParameters> pssm)
    : m_QueryFactory(query_factory), m_InternalData(new SInternalData),
    m_Options(options),  m_DbAdapter(NULL), m_DbInfo(NULL)
{
    x_Init(query_factory, options, pssm, seqsrc);
    m_InternalData->m_SeqSrc.Reset(new TBlastSeqSrc(seqsrc, 0));
}

void
CBlastPrelimSearch::SetNumberOfThreads(size_t nthreads)
{
    const bool was_multithreaded = IsMultiThreaded();

    CThreadable::SetNumberOfThreads(nthreads);
    if (was_multithreaded != IsMultiThreaded()) {
        BlastDiagnostics* diags = IsMultiThreaded()
            ? CSetupFactory::CreateDiagnosticsStructureMT()
            : CSetupFactory::CreateDiagnosticsStructure();
        m_InternalData->m_Diagnostics.Reset
            (new TBlastDiagnostics(diags, Blast_DiagnosticsFree));

        CRef<ILocalQueryData> query_data
            (m_QueryFactory->MakeLocalQueryData(&*m_Options));
        unique_ptr<const CBlastOptionsMemento> opts_memento
            (m_Options->CreateSnapshot());
        if (IsMultiThreaded())
            BlastHSPStreamRegisterMTLock(m_InternalData->m_HspStream->GetPointer(),
                                         Blast_CMT_LOCKInit());
    }
}

void
CBlastPrelimSearch::x_Init(CRef<IQueryFactory> query_factory,
                           CRef<CBlastOptions> options,
                           CConstRef<objects::CPssmWithParameters> pssm,
                           BlastSeqSrc* seqsrc,
                           size_t num_threads)
{
    CRef<SBlastSetupData> setup_data =
        BlastSetupPreliminarySearchEx(query_factory, options, pssm, seqsrc,
                                      num_threads);
    m_InternalData = setup_data->m_InternalData;
    copy(setup_data->m_Masks.begin(), setup_data->m_Masks.end(),
         back_inserter(m_MasksForAllQueries));
    m_Messages = setup_data->m_Messages;
}

int hWorldOne() {
    printf("Hello, World one!\n");
    return 0;
}

int hWorldTwo() {
    printf("Hello, World two!\n");
    return 0;
}

int
CBlastPrelimSearch::x_LaunchMultiThreadedSearch(SInternalData& internal_data)
{
    typedef vector< CRef<CPrelimSearchThread> > TBlastThreads;
    TBlastThreads the_threads(GetNumberOfThreads());

    unique_ptr<const CBlastOptionsMemento> opts_memento
        (m_Options->CreateSnapshot());
    _TRACE("Launching BLAST with " << GetNumberOfThreads() << " threads");

    // SCG
    unsigned int tThreads = 0;
    unsigned int gpuReadThreads = API_READ_THREADS_GPU;

    if(gpuReadThreads >= GetNumberOfThreads() - 1) {
        
        // Print an error message
        std::cerr << "The number of threads for reading the database is too high, compared to the threads launched!!" << std::endl;
        exit(1);
    }

    // -RMH- This appears to be a problem right now.  When used...this
    // can cause all the work to go to a single thread!  (-MN- This is fixed in SB-768)
    BlastSeqSrcSetNumberOfThreads(m_InternalData->m_SeqSrc->GetPointer(), GetNumberOfThreads());

    // This is to pass to my code, so I am going to remove the GPU threads from the total number of threads
    BlastSeqSrcSetThreads(m_InternalData->m_SeqSrc->GetPointer(), GetNumberOfThreads() - gpuReadThreads);

    int threadCount = 0; // Initialize the thread count to zero, that way you know that you should be assinging the GPU work to thread with index zero!

    // get the total number of threads
    unsigned totalThreads = GetNumberOfThreads() - gpuReadThreads;
    // create an array of threadStatus
    int * threadStatus = new int[totalThreads];
    // clear the threadStatus array using memset
    memset(threadStatus, 0, totalThreads * sizeof(int));
    
    // Pass the data to the internal data structure
    BlastSeqSrcSetTest(m_InternalData->m_SeqSrc->GetPointer(), threadStatus);

    //  Precisely measure the time taken by the code
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    /////////////////////////////////////START QUERY DATA READ//////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////////////
    // Query Preprocessing //////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////

    const std::string filename = "neighbors.bin";
    std::vector<std::vector<uint>> neighbors(TABLESIZE); // Reserve space for neighbors

    // Check if the file exists
    std::ifstream inFile(filename, std::ios::binary);
    if (inFile.good()) {
        // File exists, read the neighbors array from the file
        read_neighbors_from_file(neighbors, filename);
    } else {
        // File does not exist, calculate the neighbors array and write it to the file
        for (int a = 0; a < 28; ++a) {
            for (int b = 0; b < 28; ++b) {
                for (int c = 0; c < 28; ++c) {
                    neighbors[a * 28 * 28 + b * 28 + c] = neighbouring_words(a, b, c);
                }
            }
        }
        write_neighbors_to_file(neighbors, filename);
    }

    std::vector<std::vector<uint>> kmerIndex(TABLESIZE); // Holds the original kmerIndex
    std::vector<std::vector<uint>> kmerNeighbours(TABLESIZE); // Holds the neighbours too!
    
    uint a = static_cast<uint>(m_InternalData->m_Queries->sequence[0]);
    uint b = static_cast<uint>(m_InternalData->m_Queries->sequence[1]);
    uint c;

    for(int i = 2; i < m_InternalData->m_Queries->length; i++)
    {
        c = static_cast<uint>(m_InternalData->m_Queries->sequence[i]);
        uint index = ((a * pow(28,2)) + (b * pow(28,1)) + (c*(pow(28,0))));
        kmerIndex[index].push_back(i - 2);
        a = b;
        b = c;

    }

    // Build the neighbour table
    for (uint indx = 0; indx < kmerIndex.size(); ++indx) {
        if (!kmerIndex[indx].empty()) {
            kmerNeighbours[indx].insert(kmerNeighbours[indx].end(), kmerIndex[indx].begin(), kmerIndex[indx].end());

            for (const auto& neigh : neighbors[indx]) {
                kmerNeighbours[neigh].insert(kmerNeighbours[neigh].end(), kmerIndex[indx].begin(), kmerIndex[indx].end());
            }
        }
    }

    // Create a bitset of TABLESIZE length this is the presence vector lookip
    std::bitset<TABLESIZE> bitset;
    bitset.reset();

    for (uint indx = 0; indx < kmerNeighbours.size(); ++indx) {
        if (!kmerNeighbours[indx].empty()) {
            bitset.set(indx);
        }
    }

    // Dynamically allocate array to read the querybitmask file
    uint8_t* querybitmaskArray = new uint8_t[sizeof(bitset)];

    // Copy the bitset to the array
    memcpy(querybitmaskArray, (char*)&bitset, sizeof(bitset));

    // Get the size of the query 
    unsigned querySize = m_InternalData->m_Queries->length;

    // the bitmask size is the closest multiple of 32 greater than or equal to the query size
    const unsigned bitmaskSize = ((querySize + 31) / 32) * 32;

    size_t querydbSize = TABLESIZE * bitmaskSize/BYTE_SIZE;

    size_t querybitmaskSize = sizeof(bitset);

    // Allocate array to store the queryDbArray
    uint8_t* querydbArray = new uint8_t[TABLESIZE * bitmaskSize/BYTE_SIZE];

    for (uint kindx = 0; kindx < TABLESIZE; ++kindx) {

        std::vector<uint8_t> qb((bitmaskSize + 7) / 8, 0);
        //std::bitset<bitmaskSize> qb;
        //qb.reset();

        if (!kmerNeighbours[kindx].empty()) {
            // for (const auto& ele : kmerNeighbours[kindx]) {
            //     qb.set(ele);
            // }
            for (const auto& ele : kmerNeighbours[kindx]) {
                // set the bit at position 'ele'
                unsigned byteIndex = ele / 8;
                unsigned bitIndex = ele % 8;
                qb[byteIndex] |= (1 << bitIndex);
            }

        }

        // Copy qb into querydbArray
        size_t bytesPerKmer = qb.size(); 
        memcpy(querydbArray + kindx * bytesPerKmer, qb.data(), bytesPerKmer);

        // Copy the bitset to the array
        //memcpy(querydbArray + kindx * bitmaskSize/BYTE_SIZE, (char*)&qb, sizeof(qb));
    }
    
    uint8_t* gpuQuerydbArray;
    deviceAlloc((void**)&gpuQuerydbArray, querydbSize);

    uint8_t* gpuQuerybitmaskArray;
    deviceAlloc((void**)&gpuQuerybitmaskArray, querybitmaskSize);

    // Copy the querydbArray to the GPU
    deviceMemCpyToDevice((void*)gpuQuerydbArray, (void*)querydbArray, querydbSize);

    // Copy the querybitmaskArray to the GPU
    deviceMemCpyToDevice((void*)gpuQuerybitmaskArray, (void*)querybitmaskArray, querybitmaskSize);

    // Copy data over to the GPU constant memory, this also copies over the matrix to the array only once!
    copyConstantMemory(gpuQuerybitmaskArray, querybitmaskSize);

    uint8_t *gpuQuery;
    deviceAlloc((void**)&gpuQuery, m_InternalData->m_Queries->length);

    // Copy the query to the GPU
    deviceMemCpyToDevice((void*)gpuQuery, (void*)m_InternalData->m_Queries->sequence, m_InternalData->m_Queries->length);

    // Now pass the GPU Query DB Pointer to the BlastSeqSrc
    BlastSeqSrcSetQueryPointers(m_InternalData->m_SeqSrc->GetPointer(), (char*)gpuQuerydbArray, (char*)gpuQuery);

    // Get the total number of sequences and allocate memory for the bitmask
    uint8_t* survivingBitmask = new uint8_t[BlastSeqSrcGetNumSeqs(m_InternalData->m_SeqSrc->GetPointer())];

    // Pass the pointer to the BlastSeqSrc
    BlastSeqSrcSetSequenceBitmaskPointer(m_InternalData->m_SeqSrc->GetPointer(), (char*)survivingBitmask);

    //////////////////////////////////END QUERY DATA READ///////////////////////////

    //////////////////////////////////////////////////////////////////////////////////

    // Get the total size of the db
    uint64_t totalDBBytes = BlastSeqSrcGetTotLen(m_InternalData->m_SeqSrc->GetPointer());

    // Check if an index file exists by appending _index.bin to the totalDBBytes
    std::string indexFileName = std::to_string(totalDBBytes) + "_index.bin";

    // Print the dbName
    //cout << "TotalBytes: " << totalDBBytes << endl;

    size_t totalNumberOfSequences = BlastSeqSrcGetNumSeqs(m_InternalData->m_SeqSrc->GetPointer());

    // Check if the index file exists
    std::ifstream idx_file(indexFileName, std::ios::binary);

    // // Check if the file is valid
    if(!idx_file.is_open())
    {
        // Create a new index file with the name indexFileName
        std::ofstream index_file(indexFileName, std::ios::binary);

        BlastSeqSrcGetSeqArg seq_arg;
        memset((void*) &seq_arg, 0, sizeof(seq_arg));
        seq_arg.encoding = eBlastEncodingProtein;
        
        // This file is not valid so let us create a new one.
        uint64_t prefix_sum = 0;

        for(int i = 0; i < totalNumberOfSequences; i++)
        {
            seq_arg.oid = i;
            if (BlastSeqSrcGetSequence(m_InternalData->m_SeqSrc->GetPointer(), &seq_arg) < 0) continue;

            // Write the prefix_sum to the index file as a uint64_t
            index_file.write((char*)&prefix_sum, sizeof(uint64_t));
            prefix_sum += seq_arg.seq->length;
            BlastSeqSrcReleaseSequence(m_InternalData->m_SeqSrc->GetPointer(), &seq_arg);
        }

        // Write the last prefix_sum to the index file as a uint64_t
        index_file.write((char*)&prefix_sum, sizeof(uint64_t));

        // Close the index file
        index_file.close();

        // // Close the db file
        //db_file.close();

    }

    idx_file.close();

    //exit(1);

    /////////////////////////////////////////////////////////////////////////////////////


    uint16_t* lookupTable = new uint16_t[TABLESIZE];
    // Create a vector to store the flattened kmerNeighbours
    std::vector<uint16_t> flattenedNeighbours;
    uint16_t* lookupHelper = new uint16_t[TABLESIZE];

    // For each index print the number of entries for each index
    // Loop through the kmerNeighbours and create the flattenedNeighbours array
    for (uint indx = 0; indx < kmerNeighbours.size(); ++indx) {
        lookupTable[indx] = kmerNeighbours[indx].size();
        //std::cerr  << indx << "," << kmerNeighbours[indx].size() << std::endl;
        lookupHelper[indx] = flattenedNeighbours.size();
        flattenedNeighbours.insert(flattenedNeighbours.end(), kmerNeighbours[indx].begin(), kmerNeighbours[indx].end());
    }

    // Convert the flattenedNeighbours vector to an array
    uint16_t* flattenedNeighboursArray = new uint16_t[flattenedNeighbours.size()];
    std::copy(flattenedNeighbours.begin(), flattenedNeighbours.end(), flattenedNeighboursArray);

    // Create a GPU pointer to the lookupTable
    uint16_t* gpuLookupTable;
    deviceAlloc((void**)&gpuLookupTable, TABLESIZE * sizeof(uint16_t));
    // Copy the lookupTable to the GPU
    deviceMemCpyToDevice((void*)gpuLookupTable, (void*)lookupTable, TABLESIZE * sizeof(uint16_t));

    // Create a GPU pointer to the lookupHelper
    uint16_t* gpuLookupHelper;
    deviceAlloc((void**)&gpuLookupHelper, TABLESIZE * sizeof(uint16_t));
    // Copy the lookupHelper to the GPU
    deviceMemCpyToDevice((void*)gpuLookupHelper, (void*)lookupHelper, TABLESIZE * sizeof(uint16_t));

    // Create a GPU pointer to the flattenedNeighboursArray
    uint16_t* gpuFlattenedNeighboursArray;
    deviceAlloc((void**)&gpuFlattenedNeighboursArray, flattenedNeighbours.size() * sizeof(uint16_t));
    // Copy the flattenedNeighboursArray to the GPU
    deviceMemCpyToDevice((void*)gpuFlattenedNeighboursArray, (void*)flattenedNeighboursArray, flattenedNeighbours.size() * sizeof(uint16_t));

    // Set the lookupTable, lookupHelper and flattenedNeighboursArray to the BlastSeqSrc
    BlastSeqSrcSetIsolationPointers(m_InternalData->m_SeqSrc->GetPointer(), gpuLookupTable, gpuLookupHelper, gpuFlattenedNeighboursArray);

    int fd = open(indexFileName.c_str(), O_RDWR); // Open for read/write

    if (fd == -1) {
        perror("open");
        return EXIT_FAILURE;
    }

    struct stat file_stat;
    if (fstat(fd, &file_stat) == -1) {
        perror("fstat");
        close(fd);
        return EXIT_FAILURE;
    }

    size_t file_size = file_stat.st_size;
    size_t indexSize = file_size;

    void *mapped_memory = mmap(nullptr, file_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);

    if (mapped_memory == MAP_FAILED) {
        perror("mmap");
        close(fd);
        return EXIT_FAILURE;
    }

    uint8_t* indexArray = static_cast<uint8_t*>(mapped_memory);
    
    // Create a uint64_t pointer to access the indexArray
    uint64_t* indexArray64 = (uint64_t*)indexArray;
    uint64_t indexArraySize64 = indexSize / sizeof(uint64_t);
    uint64_t * indexArrayPointer = &indexArray64[1];


    // Read in a file in which bin boundaries are stored
    std::string binFileName = std::to_string(totalDBBytes) + "_bin_boundaries.bin";

    // Read in the bin boundaries file
    std::ifstream binFile(binFileName, std::ios::binary);

    //std::cout << "Total Number of Sequences: " << totalNumberOfSequences << std::endl;

    // Open the bin list file
    std::string blistFilename = std::to_string(totalDBBytes) + "_blist.bin";

    // Check if the file is valid, if not create it
    if(!binFile.is_open())
    {
        const size_t numElements = indexArraySize64;

        std::vector<uint64_t> original_sequence_lengths(numElements - 1);
        original_sequence_lengths.reserve(numElements - 1);

        for (size_t i = 1; i < numElements; ++i) {
            original_sequence_lengths[i - 1] = indexArray64[i] - indexArray64[i - 1];
        }

        const size_t numBins = bbins.size() - 1;
        std::vector<std::vector<size_t>> index_dict(numBins);
        for (auto &v : index_dict) {
            v.reserve(original_sequence_lengths.size() / numBins + 10); // pre-allocation 
        }
        for (size_t i = 0; i < original_sequence_lengths.size(); ++i) {
            uint64_t val = original_sequence_lengths[i];
            size_t binIndex = findBinIndex(val, bbins);
            index_dict[binIndex].push_back(i);
        }

        std::vector<size_t> starts;
        std::vector<size_t> ends;
        starts.reserve(original_sequence_lengths.size());
        ends.reserve(original_sequence_lengths.size());

        for (size_t b = 0; b < numBins; ++b) {
            const auto& bin_indices = index_dict[b];
            if (bin_indices.empty()) {
                continue;
            }
            size_t prev = bin_indices[0];
            starts.push_back(prev);
            // Detect breaks in consecutive indices
            for (size_t i = 1; i < bin_indices.size(); ++i) {
                size_t idx = bin_indices[i];
                if (idx != prev + 1) {
                    // Close off the previous contiguous segment
                    ends.push_back(prev);
                    // Start a new segment
                    starts.push_back(idx);
                }
                prev = idx;
            }
            ends.push_back(prev);
        }

        std::vector<std::pair<size_t, size_t>> ranges;
        ranges.reserve(starts.size());
        for (size_t i = 0; i < starts.size(); ++i) {
            ranges.emplace_back(starts[i], ends[i]);
        }

        std::sort(ranges.begin(), ranges.end(),
                [](const auto& a, const auto& b) {
                    return a.first < b.first; // descending
                });

        std::vector<uint32_t> nstart(ranges.size()), nend(ranges.size());
        for (size_t i = 0; i < ranges.size(); ++i) {
            nstart[i] = static_cast<uint32_t>(ranges[i].first);
            nend[i]   = static_cast<uint32_t>(ranges[i].second);
        }

        // Figure out bin numbers
        std::vector<uint32_t> blist(ranges.size());     
        for (size_t i = 0; i < ranges.size(); ++i) {
            blist[i] = findBinIndex(original_sequence_lengths[nstart[i]], bbins);
        }

        std::ofstream bFile(blistFilename, std::ios::binary);
        if (!bFile.is_open()) {
            std::cerr << "Error: Cannot open file for writing: " << blistFilename << std::endl;
            return 1;
        }

        bFile.write(reinterpret_cast<const char*>(blist.data()),
                    blist.size() * sizeof(uint32_t));

        bFile.close();

        std::ofstream outFile(binFileName, std::ios::binary);
        if (!outFile.is_open()) {
            std::cerr << "Error: Cannot open file for writing: " << binFileName << std::endl;
            return 1;
        }
        
        // Write starts
        outFile.write(reinterpret_cast<const char*>(nstart.data()),
                    nstart.size() * sizeof(uint32_t));
        // Write ends
        outFile.write(reinterpret_cast<const char*>(nend.data()),
                    nend.size() * sizeof(uint32_t));
        outFile.close();

    }

    /////////////////////////////////////////////////////////////////////////////////////
    // Attempt to read the binFile again ///////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////

    if(!binFile.is_open())
    {
        // Open the bin boundaries file
        binFile.open(binFileName, std::ios::binary);

        // Check if the file is valid
        if(!binFile.is_open())
        {
            cerr << "Error: Unable to open the bin file, after just creating it!" << endl;
            exit(1);
        }
    }

    // Get the size of the bin boundaries file
    binFile.seekg(0, binFile.end);
    size_t binSize = binFile.tellg();
    binFile.seekg(0, binFile.beg);

    // Create an array to store the bin boundaries
    uint8_t* binArray = new uint8_t[binSize];
    // Read the bin boundaries file into the array
    binFile.read((char*)binArray, binSize);
    // Close the file
    binFile.close() ;

    // Read in the bin list file
    std::ifstream blistFile(blistFilename, std::ios::binary);

    // Check if the file is valid
    if(!blistFile.is_open())
    {
        cerr << "Error: Unable to open the bin list file" << endl;
        exit(1);
    }

    // Get the size of the bin list file
    blistFile.seekg(0, blistFile.end);
    size_t blistSize = blistFile.tellg();
    blistFile.seekg(0, blistFile.beg);

    // Create an array to store the bin list
    uint8_t* blistArray = new uint8_t[blistSize];
    // Read the bin list file into the array
    blistFile.read((char*)blistArray, blistSize);
    // Close the file
    blistFile.close() ;

    // Create a uint32_t pointer to access the blistArray
    uint32_t* blistArray32 = (uint32_t*)blistArray;

    /////////////////////////////////////////////////////////////////////////////////////
    // Get the total size of GPU memory and make sure to use only portion of it, it defaults to 90%
    uint64_t spaceAvailable = (size_t)(GPU_MEM_PERCENT * getFreeGPUMemory());
    //uint64_t spaceAvailable = (uint64_t)(0.025 * getFreeGPUMemory());
    //uint64_t spaceAvailable = (uint64_t)(0.02 * getFreeGPUMemory());

    // Create a uint32_t pointer to access the binArray
    uint32_t* binArray32 = (uint32_t*)binArray;
    uint32_t binArraySize32 = binSize / sizeof(uint32_t);

    uint32_t* halfbinArray32 = binArray32 + binArraySize32/2;

    uint64_t totSeqs = 0;

    uint64_t delta = 0;

    uint64_t numChunks = 0;
    
    //for(int i = binArraySize32/2 - 1; i >= 0; i--)
    for(int i = 0; i <= binArraySize32/2 - 1; i++)
    {
        if(binArray32[i] == 0){
            delta = 0;
        }else{
            delta = indexArrayPointer[binArray32[i] - 1];
        }

        // Get the total number of bytes in the bin chunk
        uint64_t totalBytes = indexArrayPointer[halfbinArray32[i]] - delta;

        if(blistArray32[i] < CPU_OFFLOAD_BIN){
            continue;
        }

        // Print the start and the end
        // cout << "Start: " << binArray32[i] << " End: " << halfbinArray32[i] << " TotalBytes: " << totalBytes << endl;

        // Find the number of chunks
        numChunks += (totalBytes + spaceAvailable - 1) / spaceAvailable;
    }

    // Print the total number of chunks
    // cout << "TotalChunks: " << numChunks << endl;


    // Now that we have the total number of chunks, we allocate the final perThreadChunk Arrays
    uint64_t* allChunkStart = new uint64_t[numChunks*totalThreads];
    uint64_t* allChunkEnd = new uint64_t[numChunks*totalThreads];

    uint64_t chunkOffset = 0;

    uint64_t maxBytesPerThread = 0;

    uint64_t accumBytes = 0;

    // Print the ranges of the bins from the end
    for(int i = 0; i <= binArraySize32/2 - 1; i++)
    {

        if(binArray32[i] == 0){
            delta = 0;
        }else{
            delta = indexArrayPointer[binArray32[i] - 1];
        }

        if(blistArray32[i] < CPU_OFFLOAD_BIN){
            continue;
        }

        // Get the total number of bytes in the bin chunk
        uint64_t totalBytes = indexArrayPointer[halfbinArray32[i]] - delta;

        // Find the number of chunks
        uint64_t nChunks = (totalBytes + spaceAvailable - 1) / spaceAvailable;

        // Create an array to store the start and end of each chunk
        uint64_t* chunkArray = new uint64_t[nChunks];

        // Temp array to store the chunks for chunkID per thread
        uint64_t* cPT = new uint64_t[totalThreads];

        // Get address to the perThreadChunkStart
        uint64_t* perThreadChunkStartPointer = &allChunkStart[chunkOffset];

        // Get address to the perThreadChunkEnd
        uint64_t* perThreadChunkEndPointer = &allChunkEnd[chunkOffset];

        // Print the start and the end
        //cout << "Start: " << binArray32[i] << " End: " << halfbinArray32[i] << " TotalBytes: " << totalBytes << " nChunks: " << nChunks << endl;

        find_balanced_indices(&indexArrayPointer[binArray32[i]], halfbinArray32[i] - binArray32[i] + 1, nChunks, chunkArray, delta);

        uint32_t chunkVerify = 0;

        uint32_t cStart = binArray32[i];

        for(int j = 0 ; j < nChunks; j++)
        {
            uint64_t tStart;
            uint64_t tEnd;

            // Clear the chunkPerThread array
            memset(cPT, 0, totalThreads * sizeof(uint64_t));
       
            if(nChunks == 1){
                tStart = binArray32[i];
                tEnd = halfbinArray32[i];
            }
            else if(j == 0){
                tStart = binArray32[i];
                tEnd = cStart + chunkArray[0] - 1;
            }else if(j == nChunks - 1){
                tStart = cStart + chunkArray[nChunks-2];
                tEnd = halfbinArray32[i];
            }else{
                tStart = cStart + chunkArray[j-1];
                tEnd = cStart + chunkArray[j] - 1;
            }

            // Print the chunk, start and end
            //cout << "Chunk: " << j << " Start: " << tStart << " End: " << tEnd << endl;

            // Calculate the chunkVerify
            chunkVerify += tEnd - tStart + 1;

            uint64_t d1 = 0;

            if(tStart == 0){
            d1 = 0;
            }else{
                d1 = indexArrayPointer[tStart - 1];
            }

            // // // Calculate the chunkPerThread
            find_balanced_indices(&indexArrayPointer[tStart], (uint64_t)(tEnd - tStart + 1), (uint64_t)(totalThreads), cPT, d1);

            for(int k = 0; k<totalThreads; k++){

                uint64_t kS;
                uint64_t kE;

                if(totalThreads == 1){
                    kS = tStart;
                    kE = tEnd;
                }
                else if(k == 0){
                    kS = tStart;
                    kE = tStart + cPT[0] - 1;
                }else if(k == totalThreads - 1){
                    kS = tStart + cPT[totalThreads-2];
                    kE = tEnd;
                }else{
                    kS = tStart + cPT[k-1];
                    kE = tStart + cPT[k] - 1;
                }

                // Print the thread, start and end
                //cout << "Thread: " << k << " Start: " << kS << " End: " << kE << endl;

                // Assign the start and end to the perThreadChunkStart and perThreadChunkEnd
                perThreadChunkStartPointer[j*totalThreads+k] = kS;
                perThreadChunkEndPointer[j*totalThreads+k] = kE;

                uint64_t kd = 0;

                if(kS == 0){
                    kd = 0;
                }else{
                    kd = indexArrayPointer[kS - 1];
                }

                // Find the bytes per thread
                uint64_t bpt = indexArrayPointer[kE] - kd;

                maxBytesPerThread = max(maxBytesPerThread, bpt);

            }

        }

        // Update the chunkOffset
        chunkOffset += nChunks*totalThreads;

        totSeqs += halfbinArray32[i] - binArray32[i] + 1;

        if(chunkVerify != halfbinArray32[i] - binArray32[i] + 1){
            cout << "Error: Chunk Verify Failed" << endl;
            // Print the chunkVerify and the total number of sequences
            cout << "Chunk Verify: " << chunkVerify << " Total Sequences: " << halfbinArray32[i] - binArray32[i] + 1 << endl;
            exit(1);
        }

        accumBytes += totalBytes;

        // Free the chunkArray
        delete[] chunkArray;

        // Free the cPT
        delete[] cPT;
    
    }

    // Print the total number of sequences
    cout << "Total Sequences: " << totSeqs << endl;

    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    
    uint64_t numChunks_gpu = 0;
    uint64_t totSeqs_gpu = 0;


    for(int i = 0; i <= binArraySize32/2 - 1; i++)
    {
        if(binArray32[i] == 0){
            delta = 0;
        }else{
            delta = indexArrayPointer[binArray32[i] - 1];
        }

        // Get the total number of bytes in the bin chunk
        uint64_t totalBytes = indexArrayPointer[halfbinArray32[i]] - delta;

        if(blistArray32[i] >= CPU_OFFLOAD_BIN){
            continue;
        }

        // Print the start and the end
        // cout << "Start: " << binArray32[i] << " End: " << halfbinArray32[i] << " TotalBytes: " << totalBytes << endl;

        // Find the number of chunks
        numChunks_gpu += (totalBytes + spaceAvailable - 1) / spaceAvailable;
    }

    // Print the total number of chunks
    cout << "TotalChunks_gpu: " << numChunks_gpu << endl;

    //exit(1);

    // Now that we have the total number of chunks, we allocate the final perThreadChunk Arrays
    uint64_t* allChunkStart_gpu = new uint64_t[numChunks_gpu*gpuReadThreads];
    uint64_t* allChunkEnd_gpu = new uint64_t[numChunks_gpu*gpuReadThreads];

    uint64_t chunkOffset_gpu = 0;

    uint64_t maxBytesPerThread_gpu = 0;

    uint64_t accumBytes_gpu = 0;

    // Print the ranges of the bins from the end
    //for(int i = binArraySize32/2 - 1; i >= 0; i--)
    for(int i = 0; i <= binArraySize32/2 - 1; i++)
    {

        if(binArray32[i] == 0){
            delta = 0;
        }else{
            delta = indexArrayPointer[binArray32[i] - 1];
        }

        if(blistArray32[i] >= CPU_OFFLOAD_BIN){
            continue;
        }

        // Get the total number of bytes in the bin chunk
        uint64_t totalBytes = indexArrayPointer[halfbinArray32[i]] - delta;

        // Find the number of chunks
        uint64_t nChunks = (totalBytes + spaceAvailable - 1) / spaceAvailable;

        // Create an array to store the start and end of each chunk
        uint64_t* chunkArray = new uint64_t[nChunks];

        // Temp array to store the chunks for chunkID per thread
        uint64_t* cPT = new uint64_t[gpuReadThreads];

        // Get address to the perThreadChunkStart
        uint64_t* perThreadChunkStartPointer = &allChunkStart_gpu[chunkOffset_gpu];

        // Get address to the perThreadChunkEnd
        uint64_t* perThreadChunkEndPointer = &allChunkEnd_gpu[chunkOffset_gpu];

        find_balanced_indices(&indexArrayPointer[binArray32[i]], halfbinArray32[i] - binArray32[i] + 1, nChunks, chunkArray, delta);

        uint32_t chunkVerify = 0;

        uint32_t cStart = binArray32[i];

        for(int j = 0 ; j < nChunks; j++)
        {
            uint64_t tStart;
            uint64_t tEnd;

            // Clear the chunkPerThread array
            memset(cPT, 0, gpuReadThreads * sizeof(uint64_t));
       
            if(nChunks == 1){
                tStart = binArray32[i];
                tEnd = halfbinArray32[i];
            }
            else if(j == 0){
                tStart = binArray32[i];
                tEnd = cStart + chunkArray[0] - 1;
            }else if(j == nChunks - 1){
                tStart = cStart + chunkArray[nChunks-2];
                tEnd = halfbinArray32[i];
            }else{
                tStart = cStart + chunkArray[j-1];
                tEnd = cStart + chunkArray[j] - 1;
            }

            // Print the chunk, start and end
            //cout << "Chunk: " << j << " Start: " << tStart << " End: " << tEnd << endl;

            // Calculate the chunkVerify
            chunkVerify += tEnd - tStart + 1;

            uint64_t d1 = 0;

            if(tStart == 0){
            d1 = 0;
            }else{
                d1 = indexArrayPointer[tStart - 1];
            }

            // // // Calculate the chunkPerThread
            find_balanced_indices(&indexArrayPointer[tStart], (uint64_t)(tEnd - tStart + 1), (uint64_t)(gpuReadThreads), cPT, d1);

            for(int k = 0; k<gpuReadThreads; k++){

                uint64_t kS;
                uint64_t kE;

                if(gpuReadThreads == 1){
                    kS = tStart;
                    kE = tEnd;
                }
                else if(k == 0){
                    kS = tStart;
                    kE = tStart + cPT[0] - 1;
                }else if(k == gpuReadThreads - 1){
                    kS = tStart + cPT[gpuReadThreads-2];
                    kE = tEnd;
                }else{
                    kS = tStart + cPT[k-1];
                    kE = tStart + cPT[k] - 1;
                }

                // Print the thread, start and end
                //cout << "Thread: " << k << " Start: " << kS << " End: " << kE << endl;

                // Assign the start and end to the perThreadChunkStart and perThreadChunkEnd
                perThreadChunkStartPointer[j*gpuReadThreads+k] = kS;
                perThreadChunkEndPointer[j*gpuReadThreads+k] = kE;

                uint64_t kd = 0;

                if(kS == 0){
                    kd = 0;
                }else{
                    kd = indexArrayPointer[kS - 1];
                }

                // Find the bytes per thread
                uint64_t bpt = indexArrayPointer[kE] - kd;

                maxBytesPerThread_gpu = max(maxBytesPerThread_gpu, bpt);

            }

        }

        // Update the chunkOffset
        chunkOffset_gpu += nChunks*gpuReadThreads;

        totSeqs_gpu += halfbinArray32[i] - binArray32[i] + 1;

        if(chunkVerify != halfbinArray32[i] - binArray32[i] + 1){
            cout << "Error: Chunk Verify Failed" << endl;
            // Print the chunkVerify and the total number of sequences
            cout << "Chunk Verify: " << chunkVerify << " Total Sequences: " << halfbinArray32[i] - binArray32[i] + 1 << endl;
            exit(1);
        }

        accumBytes_gpu += totalBytes;

        // Free the chunkArray
        delete[] chunkArray;

        // Free the cPT
        delete[] cPT;
    
    }

    // Print the total number of sequences
    cout << "Total Sequences GPU: " << totSeqs_gpu << endl;

    uint64_t * dbOffsets_gpu = new uint64_t[gpuReadThreads];

    // Use the maxBytes to assign the dbOffsets
    for(int j=0;j<gpuReadThreads;j++){
        dbOffsets_gpu[j] = j * maxBytesPerThread_gpu;
    }

    BlastSeqSrcPassDBWorkInfoGPU(m_InternalData->m_SeqSrc->GetPointer(), &indexArray64[0], indexArraySize64, numChunks_gpu, allChunkStart_gpu, allChunkEnd_gpu, dbOffsets_gpu);

    ////////////////////////////////////////////////////////////////////////////////////////////////    
    //////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Update the spaceAvaliable if the maxBytesPerThread*totalThreads is greater than the spaceAvailable
    if(maxBytesPerThread_gpu*gpuReadThreads > spaceAvailable){
        spaceAvailable = maxBytesPerThread_gpu*gpuReadThreads;
    }

    uint64_t * dbOffsets = new uint64_t[gpuReadThreads];

    // Calculate the time taken to allocate the GPU memory
    clock_gettime(CLOCK_MONOTONIC, &end);

    // Use the maxBytes to assign the dbOffsets
    for(int j=0;j<gpuReadThreads;j++){
        dbOffsets[j] = j * maxBytesPerThread_gpu;
    }

    // Make sure that maxBytes*totalThreads is less than the spaceAvailable
    if(maxBytesPerThread_gpu*gpuReadThreads > spaceAvailable){
        cerr << "Error: The total memory required is greater than the available memory" << endl;
        exit(1);
    }
    

    // Now allocate the space on the GPU using hostAlloc so that you can perform Async Memcpy
    uint8_t* gpuDBArray; 
    deviceAlloc((void**)&gpuDBArray, maxBytesPerThread_gpu*gpuReadThreads);

    // // lets allocate the same amount of storage on the host.
    uint8_t* hostDBArray; 
    hostAlloc((void**)&hostDBArray, maxBytesPerThread_gpu*gpuReadThreads);

    // Now pass the GPU DB Pointer to the BlastSeqSrc
    BlastSeqSrcSetDBPointer(m_InternalData->m_SeqSrc->GetPointer(), (char*)gpuDBArray);

    // Now pass the Host DB Pointer to the BlastSeqSrc
    BlastSeqSrcSetHostDBPointer(m_InternalData->m_SeqSrc->GetPointer(), (char*)hostDBArray);

    // Print the pointers for debugging
    //cout << "GPU DB Pointer: " << (void*)gpuDBArray << endl;
    //cout << "Host DB Pointer: " << (void*)hostDBArray << endl;

    
    // Calculate the time taken in seconds
    // double tt = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1000000000.0;

    // Print the time taken in seconds
    //cout << "Time to allocate GPU memory: " << tt << " seconds" << endl;

    // Once read I can pass the indexArray64, indexArraySize64, chunkArray, numChunks to the GPU
    BlastSeqSrcPassDBWorkInfo(m_InternalData->m_SeqSrc->GetPointer(), &indexArray64[0], indexArraySize64, numChunks, allChunkStart, allChunkEnd, dbOffsets);

    /////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////

    // Stop the timer
    clock_gettime(CLOCK_MONOTONIC, &end);

    // Calculate the time taken
    double time_taken = (end.tv_sec - start.tv_sec) * 1e9;

    time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;

    // // Print the time taken
    // cout << "Time taken for this setup: " << time_taken << " seconds" << endl;

    // Create the threads ...
    NON_CONST_ITERATE(TBlastThreads, thread, the_threads) {

        // Assign thread count to internal data
        internal_data.threadIndex = threadCount;

        //internal_data.ptrToSeqsrcArray = &seqsrcArray;

        // Assign the current thread to the BlastSeqSrc
        BlastSeqSrcSetCurrentThread(m_InternalData->m_SeqSrc->GetPointer(), threadCount);

        // Set the thread Type.

        // If the current thread count is less than 2, then set the pFlag to 1
        if(tThreads < gpuReadThreads){
            internal_data.pflag = 1;
            internal_data.threadCount = gpuReadThreads;
            // Print the thread count
            // std::cerr << "GPU Thread No: " << threadCount << std::endl;
            BlastSeqSrcSetThreadType(m_InternalData->m_SeqSrc->GetPointer(), 1);
        }else{
            internal_data.pflag = 0;
            internal_data.threadCount = GetNumberOfThreads() - gpuReadThreads;
            // Print the thread count
            // std::cerr << "CPU Thread No: " << threadCount << std::endl;
            BlastSeqSrcSetThreadType(m_InternalData->m_SeqSrc->GetPointer(), 0);
        }

        thread->Reset(new CPrelimSearchThread(internal_data,
                                              opts_memento.get()));
        if (thread->Empty()) {
            NCBI_THROW(CBlastSystemException, eOutOfMemory,
                       "Failed to create preliminary search thread");
        }

        threadCount++;
        // Increment the thread count
        if(tThreads == gpuReadThreads - 1)
            threadCount = 0;
        
        
        tThreads++;
    }

    // Inform indexing library about the number of concurrent
    // search threads.
    // GetDbIndexSetNumThreadsFn()( GetNumberOfThreads() );
    GetDbIndexSetNumThreadsFn()( GetNumberOfThreads() ); //[SCG]

    clock_gettime(CLOCK_MONOTONIC, &start);

    // ... launch the threads ...
    NON_CONST_ITERATE(TBlastThreads, thread, the_threads) {
        (*thread)->Run();
    }

    // ... and wait for the threads to finish
    Uint8 retv(0);
    for (size_t i = 0; i < the_threads.size(); i++) {
        void* result(0);

        the_threads[i]->Join(&result);
        
        if (result) {
            // Thread is not really returning a pointer, it's actually returning an int
            retv = reinterpret_cast<Uint8>(result);
        }
    }


    // Print the thread status to stderr
    for(int i = 0; i < totalThreads; i++)
    {
        cerr << "Thread Status: " << threadStatus[i] << endl;
    }

    // Free the arrays
    delete [] threadStatus;
    delete [] allChunkStart;
    delete [] allChunkEnd;
    delete [] dbOffsets;
    delete [] survivingBitmask;

    BlastSeqSrcSetNumberOfThreads(m_InternalData->m_SeqSrc->GetPointer(), 0);

    if (retv) {
          NCBI_THROW(CBlastException, eCoreBlastError,
                                   BlastErrorCode2String((Int2)retv));
    }
    return 0;
}

CRef<SInternalData>
CBlastPrelimSearch::Run()
{
    if (! BlastSeqSrcGetNumSeqs(m_InternalData->m_SeqSrc->GetPointer())) {
        string msg =
            "Filtering resulted in an empty database.";

        m_Messages.AddMessageAllQueries(eBlastSevWarning,
                                        kBlastMessageNoContext,
                                        msg);
    }

    BlastSeqSrcResetChunkIterator(m_InternalData->m_SeqSrc->GetPointer());

    CEffectiveSearchSpacesMemento eff_memento(m_Options);
    SplitQuery_SetEffectiveSearchSpace(m_Options, m_QueryFactory,
                                       m_InternalData);
    int retval = 0;

    unique_ptr<const CBlastOptionsMemento> opts_memento
        (m_Options->CreateSnapshot());
    BLAST_SequenceBlk* queries = m_InternalData->m_Queries;
    LookupTableOptions * lut_options = opts_memento->m_LutOpts;
    BlastInitialWordOptions * word_options = opts_memento->m_InitWordOpts;

    // Query splitting data structure (used only if applicable)
    CRef<SBlastSetupData> setup_data(new SBlastSetupData(m_QueryFactory, m_Options));
    CRef<CQuerySplitter> query_splitter = setup_data->m_QuerySplitter;
    if (query_splitter->IsQuerySplit()) {

        CRef<CSplitQueryBlk> split_query_blk = query_splitter->Split();

        for (Uint4 i = 0; i < query_splitter->GetNumberOfChunks(); i++) {
            try {
                CRef<IQueryFactory> chunk_qf =
                    query_splitter->GetQueryFactoryForChunk(i);
                _TRACE("Query chunk " << i << "/" <<
                       query_splitter->GetNumberOfChunks());
                CRef<SInternalData> chunk_data =
                    SplitQuery_CreateChunkData(chunk_qf, m_Options,
                                               m_InternalData,
                                               GetNumberOfThreads());

                CRef<ILocalQueryData> query_data(
                        chunk_qf->MakeLocalQueryData( &*m_Options ) );
                BLAST_SequenceBlk * chunk_queries =
                    query_data->GetSequenceBlk();
                GetDbIndexSetUsingThreadsFn()( IsMultiThreaded() );
                GetDbIndexRunSearchFn()(
                        chunk_queries, lut_options, word_options );

                if (IsMultiThreaded()) {
                     x_LaunchMultiThreadedSearch(*chunk_data);
                } else {
                    retval =
                        CPrelimSearchRunner(*chunk_data, opts_memento.get())();
                    if (retval) {
                        NCBI_THROW(CBlastException, eCoreBlastError,
                                   BlastErrorCode2String(retval));
                    }
                }


                _ASSERT(chunk_data->m_HspStream->GetPointer());
                BlastHSPStreamMerge(split_query_blk->GetCStruct(), i,
                                chunk_data->m_HspStream->GetPointer(),
                                m_InternalData->m_HspStream->GetPointer());
                _ASSERT(m_InternalData->m_HspStream->GetPointer());
                // free this as the query_splitter keeps a reference to the
                // chunk factories, which in turn keep a reference to the local
                // query data.
                query_data->FlushSequenceData();
            } catch (const CBlastException& e) {
                // This error message is safe to ignore for a given chunk,
                // because the chunks might end up producing a region of
                // the query for which ungapped Karlin-Altschul blocks
                // cannot be calculated
                const string err_msg1("search cannot proceed due to errors "
                                     "in all contexts/frames of query "
                                     "sequences");
                const string err_msg2(kBlastErrMsg_CantCalculateUngappedKAParams);
                if (e.GetMsg().find(err_msg1) == NPOS && e.GetMsg().find(err_msg2) == NPOS) {
                    throw;
                }
            }
        }

        // Restore the full query sequence for the traceback stage!
        if (m_InternalData->m_Queries == NULL) {
            CRef<ILocalQueryData> query_data
                (m_QueryFactory->MakeLocalQueryData(&*m_Options));
            // Query masking info is calculated as a side-effect
            CBlastScoreBlk sbp
                (CSetupFactory::CreateScoreBlock(opts_memento.get(), query_data,
                                        NULL, m_Messages, NULL, NULL));
            m_InternalData->m_Queries = query_data->GetSequenceBlk();
        }
    } else {

        GetDbIndexSetUsingThreadsFn()( IsMultiThreaded() );
        GetDbIndexRunSearchFn()( queries, lut_options, word_options );

        if (IsMultiThreaded()) {
             x_LaunchMultiThreadedSearch(*m_InternalData);
        } else {
            retval = CPrelimSearchRunner(*m_InternalData, opts_memento.get())();
            if (retval) {
                NCBI_THROW(CBlastException, eCoreBlastError,
                           BlastErrorCode2String(retval));
            }
        }
    }

    return m_InternalData;
}

int
CBlastPrelimSearch::CheckInternalData()
{
    int retval = 0;
    retval = BlastScoreBlkCheck(m_InternalData->m_ScoreBlk->GetPointer());
    return retval;
}


BlastHSPResults*
CBlastPrelimSearch::ComputeBlastHSPResults(BlastHSPStream* stream,
                                           Uint4 max_num_hsps,
					   bool* rm_hsps,
					   vector<bool> *rm_hsps_info) const
{
    bool any_query_hsp_limited = false;
    unique_ptr<const CBlastOptionsMemento> opts_memento
        (m_Options->CreateSnapshot());

    _ASSERT(m_InternalData->m_QueryInfo->num_queries > 0);
    Boolean *removed_hsps = new Boolean [ m_InternalData->m_QueryInfo->num_queries ];
    SBlastHitsParameters* hit_param = NULL;
    SBlastHitsParametersNew(opts_memento->m_HitSaveOpts,
                            opts_memento->m_ExtnOpts,
                            opts_memento->m_ScoringOpts,
                            &hit_param);
    BlastHSPResults* retval =
        Blast_HSPResultsFromHSPStreamWithLimitEx(stream,
           (Uint4) m_InternalData->m_QueryInfo->num_queries,
           hit_param,
           max_num_hsps,
           removed_hsps);
    if( rm_hsps_info){
	rm_hsps_info->reserve(m_InternalData->m_QueryInfo->num_queries );
        for( int query_index = 0 ; query_index < m_InternalData->m_QueryInfo->num_queries ; query_index ++ ){
	  (*rm_hsps_info)[ query_index ] = removed_hsps[query_index] == FALSE ? false : true;
	  if( (*rm_hsps_info)[ query_index ] ) any_query_hsp_limited = true;
        }
    }
    delete [] removed_hsps;
    if( rm_hsps ) *rm_hsps = any_query_hsp_limited ;
    // applications assume the HSPLists in the HSPResults are
    // sorted in order of worsening best e-value
    Blast_HSPResultsSortByEvalue(retval);
    return retval;
}

bool CBlastPrelimSearch::Run( vector<list<CRef<CStd_seg> > >  & l )
{
	Run();
	return x_BuildStdSegList(l);
}

void s_FixNumIdent(BlastHSPList *hsp_list, bool gapped_calculation)
{
	   BlastHSP* hsp;
	   int i;

	   for (i=0; i < hsp_list->hspcnt; i++)
	   {
	      hsp = hsp_list->hsp_array[i];
	      if (gapped_calculation)
	    	  hsp->num_ident = -1;
	   }
}

void s_GetBitScores(BlastHitList * hit_list, bool gapped_calculation, const BlastScoreBlk * sbp)
{

	 for (int i = 0; i < hit_list->hsplist_count; i++)
	 {
	    	BlastHSPList* hsp_list = hit_list->hsplist_array[i];
		    if (!hsp_list)
		    	continue;

		    Blast_HSPListGetBitScores(hsp_list, gapped_calculation, sbp);
		    s_FixNumIdent(hsp_list, gapped_calculation);
	 }
}

// Results is trimmed by Blast Hits Save Options if set
bool CBlastPrelimSearch::x_BuildStdSegList( vector<list<CRef<CStd_seg> > >  & l )
{
	if(m_InternalData->m_HspStream.Empty())
	{
		_TRACE("HSP Stream is empty");
		return false;
	}

	if(NULL != m_DbInfo)
	{
		m_DbAdapter.Reset(new CLocalDbAdapter(*m_DbInfo));
	}

	if(m_DbAdapter.Empty())
	{
		_TRACE("This method does not support CBlastPrelimSearch constructed with BlastSeqSrc");
		return false;
	}

	BlastHSPStream * hsp_stream  = m_InternalData->m_HspStream->GetPointer();
	if(NULL == hsp_stream)
	{
		_TRACE("NULL HSP Stream Pointer");
		return false;
	}

	IBlastSeqInfoSrc * s_seqInfoSrc = m_DbAdapter->MakeSeqInfoSrc();
	EBlastProgramType program = hsp_stream->program;

	CStructWrapper<BlastHSPResults> results
            (ComputeBlastHSPResults(hsp_stream), Blast_HSPResultsFree);

	if(NULL == results.GetPointer())
		return false;

	int num_queries = results->num_queries;

	BlastHitList ** q_list_ptr = results->hitlist_array;
	CRef<ILocalQueryData> local_query_data = m_QueryFactory->MakeLocalQueryData(m_Options.GetPointer());
	l.resize(num_queries);
	const BlastScoreBlk * sbp = m_InternalData->m_ScoreBlk->GetPointer();
	bool gapped_cal = m_Options->GetGappedMode();
	for(int i=0; i < num_queries; i++)
	{
		CConstRef<CSeq_loc> query_loc = local_query_data->GetSeq_loc(i);
		TSeqPos	query_length = local_query_data->GetSeqLength(i);
		BlastHitList * hit_list = q_list_ptr[i];
		if(NULL != hit_list)
		{

			s_GetBitScores(hit_list, gapped_cal, sbp);
			BLASTPrelminSearchHitListToStdSeg(program, hit_list, *query_loc, query_length, s_seqInfoSrc, l[i]);
		}
	}

	return true;
}
END_SCOPE(blast)
END_NCBI_SCOPE

/* @} */

