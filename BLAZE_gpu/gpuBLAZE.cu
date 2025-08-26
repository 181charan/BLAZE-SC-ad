#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <omp.h>
#include <cmath>
#include <numeric>
#include <bitset>
#include <chrono>
#include <unistd.h>
#include <iostream>
#include "util.hpp"

#include <algo/blast/core/gpu_cpu_common.h>

typedef struct {
    int best;
    int best_gap;
} BlastGapDP;

#define X_DROPOFF 38
#define GAPPED_CUTOFF 67

#define WARP_SIZE 32
#define BYTE_SIZE 8

#define kRestrictedMult 0.68

#define QUERY_BITSET_SIZE 256
#define SUBJECT_LEN 256
#define THREADS_PER_PROBLEM 64 
#define MAX_SEEDS 8 
#define MAX_SEEDS_GAPPED 8 
#define WARPS_PER_BLOCK_FLOAT (static_cast<float>(THREADS_PER_PROBLEM) / static_cast<float>(WARP_SIZE)) 
#define WARPS_PER_BLOCK THREADS_PER_PROBLEM/WARP_SIZE
#define BITSET_BYTES QUERY_BITSET_SIZE/BYTE_SIZE
#define DIAGS QUERY_BITSET_SIZE+SUBJECT_LEN-1
#define MAX_UG_SEEDS 128

// Defaults from BLASTP
//////////////////////////////////////////////
//////////////////////////////////////////////
#define ndropoff -16
#define dropoff 16
#define window 40
#define wordsize 3
#define cutoff 41
#define HSP_MAX_WINDOW 11
#define RESTRICT_SIZE 10

#define gapOpen -11
#define gapExtend -1

#define GAP_OPEN 11
#define GAP_EXTEND 1
#define GAP_OPEN_EXTEND 12
#define MATCH_SCORE 1
#define MISMATCH_SCORE -1
#define MININT (-2147483647-1)/2
//////////////////////////////////////////////
//////////////////////////////////////////////

using namespace std;

__constant__ int score_matrix[28][28];
__constant__ uint8_t query_bitmask[2744];

__device__ inline unsigned generateMask(int laneID) {
    // If laneID is 31, return all 1s (32 bits)
    if (laneID == 31) {
        return 0xFFFFFFFF;
    }
    // Otherwise, set bits up to laneID
    return (1U << (laneID + 1)) - 1U;
}

__inline__ __device__
unsigned semiGappedKernelLeft(const char *A, const char *B, int M, int N, BlastGapDP * score_array, bool reverse_sequence){

    int sizeA = M;
    int sizeB = N;

    // Get the thread id
    int tid = threadIdx.x;

    // laneID within the warp
    int laneID = tid % WARP_SIZE;

    int warpID = tid / WARP_SIZE;

    char * A_block = (char *)A;
    char * B_block = (char *)B;
    char * tempBlockPtr;

    // min of A and B
    int minAB = min(sizeA, sizeB);
    int maxAB = max(sizeA, sizeB);

    if(minAB == sizeA){
        tempBlockPtr = A_block;
        A_block = B_block;
        B_block = tempBlockPtr;
        // swap the sizes
        int tempSize = sizeA;
        sizeA = sizeB;
        sizeB = tempSize;
    }

    int sweepcnt = ceil(((double)maxAB/(double)WARP_SIZE));

    int temp1 = 0;
    int temp2 = 0;

    int globalMax = 0;          // globalMax is the maximum score encountered in the DP calculation.

    int score;
    int score_gap_row;
    int score_gap_col;
    int next_score;
    int best_score = MININT;

    score = -GAP_OPEN_EXTEND;
    score_array[0].best = 0;
    score_array[0].best_gap = -GAP_OPEN_EXTEND;

    if (laneID == 0) {
        for (int i = 1; i <= minAB; i++) {
            if (score < -X_DROPOFF) {
                score_array[i].best = MININT;
                score_array[i].best_gap = MININT;
            } else {
                score_array[i].best = score;
                score_array[i].best_gap = score - GAP_OPEN_EXTEND;
                score -= GAP_EXTEND;
            }
        }
    }

    __syncwarp();

    int best_prev = 0, bestgap_prev = 0;
    int b_first_index = 0, b_last_index = minAB;
    int prev_b_first_index = 0, prev_b_last_index = minAB;
    int b_started = 0, b_ended = 0;

    for (int e = 0; e < sweepcnt; e++) {
        score = MININT;
        int score_gap_row = MININT;

        best_score = __reduce_max_sync(0xffffffff, best_score);
        b_first_index = __reduce_max_sync(0xffffffff, b_first_index);
        b_last_index = __reduce_max_sync(0xffffffff, b_last_index);

        prev_b_first_index = b_first_index;
        prev_b_last_index = b_last_index;

        if (b_first_index >= minAB || b_last_index < b_first_index) break;

        b_started = 0;
        b_ended = 0;

        for (int x = (e * WARP_SIZE); x < (e * WARP_SIZE + WARP_SIZE + minAB); x++) {
            int i = e * WARP_SIZE + laneID;
            int j = x - i;

            int temp1 = __shfl_up_sync(0xffffffff, best_prev, 1, 32);
            int temp2 = __shfl_up_sync(0xffffffff, bestgap_prev, 1, 32);

            int tmp_best_score = __shfl_up_sync(0xffffffff, best_score, 1, 32);
            if (laneID != 0) best_score = max(tmp_best_score, best_score);

            unsigned int term_mask = __ballot_sync(0xffffffff, b_ended);
            if (term_mask == 0xFFFFFFFF) break;

            if (j >= prev_b_first_index) {
                if (!(i < 0 || (i + 1) > maxAB || j < 0 || (j + 1) > minAB)) {
                    int i_temp = i + 1;
                    int b_index = j;

                    int next_score, score_gap_col;
                    if (laneID == 0) {
                        if (b_index <= prev_b_last_index) {
                            next_score = score_array[b_index].best;
                            score_gap_col = score_array[b_index].best_gap;
                        } else {
                            next_score = MININT;
                            score_gap_col = MININT;
                        }
                    } else {
                        next_score = temp1;
                        score_gap_col = temp2;
                    }

                    next_score += score_matrix[A_block[maxAB - i_temp]][B_block[minAB - (b_index + 1)]];

                    if (score < score_gap_col) score = score_gap_col;
                    if (score < score_gap_row) score = score_gap_row;

                    if (best_score - score > X_DROPOFF) {
                        if (!b_started) {
                            b_first_index++;
                        } else if (!b_ended) {
                            b_last_index = b_index;
                            b_ended = 1;
                        }
                        if(b_started){
                            best_prev = MININT;
                        }
                    } else {
                        b_started = 1;
                        if (score > best_score) best_score = score;

                        score_gap_row -= GAP_EXTEND;
                        score_gap_col -= GAP_EXTEND;
                        bestgap_prev = max(score - GAP_OPEN_EXTEND, score_gap_col);
                        score_gap_row = max(score - GAP_OPEN_EXTEND, score_gap_row);
                        best_prev = score;
                    }

                    score = next_score;

                    if (laneID == 31) {
                        score_array[b_index].best = best_prev;
                        score_array[b_index].best_gap = bestgap_prev;
                    }
                }
            }
        }
    }

    // holds the global maximum across all threads
    globalMax = __reduce_max_sync( 0xFFFFFFFF, best_score);

    return globalMax;
    
}

__inline__ __device__
unsigned semiGappedKernelRight(const char *A, const char *B, int M, int N, BlastGapDP * score_array, bool reverse_sequence){

    int sizeA = M;
    int sizeB = N;

    // Get the thread id
    int tid = threadIdx.x;

    // laneID within the warp
    int laneID = tid % WARP_SIZE;

    int warpID = tid / WARP_SIZE;

    char * A_block = (char *)A;
    char * B_block = (char *)B;
    char * tempBlockPtr;

    // min of A and B
    int minAB = min(sizeA, sizeB);
    int maxAB = max(sizeA, sizeB);

    if(minAB == sizeA){
        tempBlockPtr = A_block;
        A_block = B_block;
        B_block = tempBlockPtr;
        // swap the sizes
        int tempSize = sizeA;
        sizeA = sizeB;
        sizeB = tempSize;
    }

    int sweepcnt = ceil(((double)maxAB/(double)WARP_SIZE));

    int temp1 = 0;
    int temp2 = 0;

    int globalMax = 0;          // globalMax is the maximum score encountered in the DP calculation.

    int score;
    int score_gap_row;
    int score_gap_col;
    int next_score;
    int best_score = MININT;

    score = -GAP_OPEN_EXTEND;
    score_array[0].best = 0;
    score_array[0].best_gap = -GAP_OPEN_EXTEND;

    if (laneID == 0) {
        for (int i = 1; i <= minAB; i++) {
            if (score < -X_DROPOFF) {
                score_array[i].best = MININT;
                score_array[i].best_gap = MININT;
            } else {
                score_array[i].best = score;
                score_array[i].best_gap = score - GAP_OPEN_EXTEND;
                score -= GAP_EXTEND;
            }
        }
    }

    __syncwarp();

    int best_prev = 0, bestgap_prev = 0;
    int b_first_index = 0, b_last_index = minAB;
    int prev_b_first_index = 0, prev_b_last_index = minAB;
    int b_started = 0, b_ended = 0;

    for (int e = 0; e < sweepcnt; e++) {
        score = MININT;
        int score_gap_row = MININT;

        best_score = __reduce_max_sync(0xffffffff, best_score);
        b_first_index = __reduce_max_sync(0xffffffff, b_first_index);
        b_last_index = __reduce_max_sync(0xffffffff, b_last_index);

        prev_b_first_index = b_first_index;
        prev_b_last_index = b_last_index;

        if (b_first_index >= minAB || b_last_index < b_first_index) break;

        b_started = 0;
        b_ended = 0;

        for (int x = (e * WARP_SIZE); x < (e * WARP_SIZE + WARP_SIZE + minAB); x++) {
            int i = e * WARP_SIZE + laneID;
            int j = x - i;

            int temp1 = __shfl_up_sync(0xffffffff, best_prev, 1, 32);
            int temp2 = __shfl_up_sync(0xffffffff, bestgap_prev, 1, 32);

            int tmp_best_score = __shfl_up_sync(0xffffffff, best_score, 1, 32);
            if (laneID != 0) best_score = max(tmp_best_score, best_score);

            unsigned int term_mask = __ballot_sync(0xffffffff, b_ended);
            if (term_mask == 0xFFFFFFFF) break;

            if (j >= prev_b_first_index) {
                if (!(i < 0 || (i + 1) > maxAB || j < 0 || (j + 1) > minAB)) {
                    int i_temp = i + 1;
                    int b_index = j;

                    int next_score, score_gap_col;
                    if (laneID == 0) {
                        if (b_index <= prev_b_last_index) {
                            next_score = score_array[b_index].best;
                            score_gap_col = score_array[b_index].best_gap;
                        } else {
                            next_score = MININT;
                            score_gap_col = MININT;
                        }
                    } else {
                        next_score = temp1;
                        score_gap_col = temp2;
                    }

                    next_score += score_matrix[A_block[i_temp]][B_block[b_index+1]];

                    if (score < score_gap_col) score = score_gap_col;
                    if (score < score_gap_row) score = score_gap_row;

                    if (best_score - score > X_DROPOFF) {
                        if (!b_started) {
                            b_first_index++;
                        } else if (!b_ended) {
                            b_last_index = b_index;
                            b_ended = 1;
                        }
                        if(b_started){
                            best_prev = MININT;
                        }
                        
                    } else {
                        b_started = 1;
                        if (score > best_score) best_score = score;

                        score_gap_row -= GAP_EXTEND;
                        score_gap_col -= GAP_EXTEND;
                        bestgap_prev = max(score - GAP_OPEN_EXTEND, score_gap_col);
                        score_gap_row = max(score - GAP_OPEN_EXTEND, score_gap_row);
                        best_prev = score;
                    }

                    score = next_score;

                    if (laneID == 31) {
                        score_array[b_index].best = best_prev;
                        score_array[b_index].best_gap = bestgap_prev;
                    }
                }
            }
        }
    }

    // holds the global maximum across all threads
    globalMax = __reduce_max_sync( 0xFFFFFFFF, best_score);

    return globalMax;
    
}

__inline__ __device__
void s_BlastAaExtendLeft( const char * s, const char * q, uint16_t s_off, uint16_t q_off , int n , int * left_score_o, int * left_d_o, int laneID, uint32_t problemID){

    s += s_off - n;
    q += q_off - n;

    int left_d = 0;
    int left_score = 0;

    uint loopcnt = ceil(((double)(n+1)/(double)WARP_SIZE)); 

    int globalMax = 0;
    int globalMaxD = n+1;
    int personalMax = 0;

    int prevScore = 0;

    for(int l=0;l<loopcnt;l++){
        
            int j = n - laneID - l*WARP_SIZE;

            bool valid = j >= 0;

            // Get a mask of all threads which are valid.
            unsigned valid_mask = __ballot_sync(0xFFFFFFFF, valid);

            int localScore = valid ? score_matrix[q[j]][s[j]] : MININT;

            for (int offset = 1; offset < warpSize; offset *= 2) { 
                int value = __shfl_up_sync(0xFFFFFFFF, localScore, offset);
                if (laneID >= offset) {
                    localScore += value;
                }
            }

            localScore += prevScore;
            
            if(laneID == 0)
            personalMax = max(personalMax, localScore); 
            else
                personalMax =  localScore;

            for (int offset = 1; offset < warpSize; offset *= 2) {
                int candidate = __shfl_up_sync(0xFFFFFFFF, personalMax, offset);
                if (laneID >= offset) {
                    personalMax = max(personalMax, candidate);
                }
            }

            unsigned drop_mask = __ballot_sync(valid_mask, ((personalMax - localScore) >= dropoff) ); // This is the dropflag array from the "cuBLASTP Paper"
            unsigned non_zero_dm = __ffs(drop_mask) - 1;

            if(drop_mask != 0){
                
                int theAbsoluteMax = __shfl_sync(0xFFFFFFFF, personalMax, non_zero_dm); 

                if(theAbsoluteMax == globalMax){
                    left_d = globalMaxD + 1;
                    left_score = globalMax;
                    break;
                }else{

                    unsigned tmask = generateMask(non_zero_dm);
                    unsigned smask = __ballot_sync(0xFFFFFFFF, localScore == theAbsoluteMax); 
                    smask = smask & tmask; 
                    uint8_t tempID = __ffs(smask) - 1; 
                    left_score = theAbsoluteMax; 
                    left_d = l*WARP_SIZE + tempID + 1;
                    break;
                }
        
            }else{

                int thread31_max = __shfl_sync(0xFFFFFFFF, personalMax, 31); 
                if(thread31_max != globalMax){
                    unsigned smask = __ballot_sync(0xFFFFFFFF, localScore == thread31_max);
                    uint8_t maxID = __ffs(smask) - 1;
                    globalMax = thread31_max;
                    globalMaxD = l*WARP_SIZE + maxID; 
                } 

                if(l == loopcnt-1){
                    left_d = globalMaxD + 1;
                    left_score = globalMax;
                    break;
                }

                prevScore = __shfl_sync(0xFFFFFFFF, localScore, (warpSize - 1)); // this is the cumulative SUM
                personalMax = __shfl_sync(0xFFFFFFFF, personalMax, (warpSize - 1)); // this is the cumulative SUM
                
            }    
        
    }

    __syncwarp();

    *left_d_o = __shfl_sync(0xFFFFFFFF, left_d, 0);
    *left_score_o = __shfl_sync(0xFFFFFFFF, left_score, 0); 

}

__inline__ __device__
void s_BlastAaExtendRight( const char * st, const char * qt, uint16_t s_off, uint16_t q_off, int n , int score, int * right_score_o, int * right_d_o, int laneID, uint32_t problemID){

    st += s_off + 1;
    qt += q_off + 1;
    
    int right_d = 0;
    int right_score = 0;

    int prevScore = score;

    int personalMax = score;
    int personalMaxD = -1;

    uint loopcnt = ceil(((double)(n+1)/(double)WARP_SIZE)); 
    int globalMax = score;
    int globalMaxD = -1;

    int currentLocalMax = score;

    for(int l=0;l<loopcnt;l++){
        
        int j = l*WARP_SIZE + laneID;
        
            bool valid = j < n;

            unsigned valid_mask = __ballot_sync(0xFFFFFFFF, valid);

            int localScore = valid ? score_matrix[qt[j]][st[j]] : MININT;
            
            for (int offset = 1; offset < warpSize; offset *= 2) { 
                int value = __shfl_up_sync(0xFFFFFFFF, localScore, offset);
                if (laneID >= offset) {
                    localScore += value;
                }
            }

            localScore += prevScore;
            
            if(laneID == 0)
            personalMax = max(personalMax, localScore);  
            else
                personalMax =  localScore; 
            
            for (int offset = 1; offset < warpSize; offset *= 2) {
                int candidate = __shfl_up_sync(0xFFFFFFFF, personalMax, offset);
                if (laneID >= offset) {
                    personalMax = max(personalMax, candidate);
                }
            }

            unsigned drop_mask = __ballot_sync(valid_mask, ((localScore <= 0) || (personalMax - localScore) >= dropoff) ); 
            unsigned non_zero_dm = __ffs(drop_mask) - 1;

            if(drop_mask != 0){
                
                int theAbsoluteMax = __shfl_sync(0xFFFFFFFF, personalMax, non_zero_dm); 
                
                if(theAbsoluteMax == globalMax){
                    right_d = globalMaxD + 1;
                    right_score = globalMax;
                    break;
                }else{
                    unsigned tmask = generateMask(non_zero_dm); 
                    unsigned smask = __ballot_sync(0xFFFFFFFF, localScore == theAbsoluteMax); 
                    smask = smask & tmask; 
                    uint8_t tempID = __ffs(smask) - 1; 

                    right_score = theAbsoluteMax; 
                    right_d = l*WARP_SIZE + tempID + 1;

                    break;
                }
        
            }else{

                int thread31_max = __shfl_sync(0xFFFFFFFF, personalMax, 31); 

                if(thread31_max != globalMax){
                    unsigned smask = __ballot_sync(0xFFFFFFFF, localScore == thread31_max);
                    uint8_t maxID = __ffs(smask) - 1;
                    globalMax = thread31_max;
                    globalMaxD = l*WARP_SIZE + maxID;                     
                } 

                if(l == loopcnt-1){
                    right_d = globalMaxD + 1;
                    right_score = globalMax;
                    break;
                }

                prevScore = __shfl_sync(0xFFFFFFFF, localScore, (warpSize - 1)); // this is the cumulative SUM
                personalMax = __shfl_sync(0xFFFFFFFF, personalMax, (warpSize - 1)); // this is the cumulative SUM

            }    
        
    }

    __syncwarp();

    *right_d_o = __shfl_sync(0xFFFFFFFF, right_d, 0);
    *right_score_o = __shfl_sync(0xFFFFFFFF, right_score, 0); 

}

template<unsigned int T_QUERY_BITSET_SIZE, unsigned int T_SUBJECT_LEN, unsigned int T_THREADS_PER_PROBLEM>
__global__ void gpuBLAZE_full(uint8_t * gpuDbArray, uint32_t * gpuSequenceLength, uint8_t * gpuQuerydbArray, char * gpuQuery, size_t gpuQueryLength, uint32_t streamID, uint8_t * gpuNumSeeds, int ungapped_cutoff, int gapped_cutoff, uint32_t * gappedSeeds, uint32_t * gpuHSPLen, uint16_t * query_lookup, uint16_t *query_prefix, uint16_t *seed_lookup_table)
{
    // Calculate the constants during compile time
    constexpr float T_WARPS_PER_BLOCK_FLOAT = static_cast<float>(T_THREADS_PER_PROBLEM) / WARP_SIZE;
    constexpr int T_WARPS_PER_BLOCK = T_THREADS_PER_PROBLEM / WARP_SIZE;
    constexpr int T_BITSET_BYTES = (T_QUERY_BITSET_SIZE + 7) / 8;
    constexpr int T_DIAGS = T_QUERY_BITSET_SIZE + T_SUBJECT_LEN - 1;

    // Hold the diagonal private information about the lasthit and the flag
    __shared__ uint32_t diagInfo[T_DIAGS]; // M+N-1 = T_DIAGS

    __shared__ uint32_t finalSeeds[T_WARPS_PER_BLOCK*MAX_UG_SEEDS]; // Holds the seeds that each thread finds.
    __shared__ uint16_t lastHits[T_WARPS_PER_BLOCK*MAX_UG_SEEDS]; // Holds the seeds that each thread finds.
        
    __shared__ uint16_t warpTotals[T_WARPS_PER_BLOCK]; // Holds the total number of seeds found by each warp.

    __shared__ uint8_t moreSeedsFlag[T_WARPS_PER_BLOCK]; // Holds information on seeds overflow
    
    uint8_t seedCnt = 0; // Local seed count for each thread

    int warpID = threadIdx.x / WARP_SIZE;

    uint8_t laneID = threadIdx.x & 0x1F; // Gets the lane ID.

    int tid = threadIdx.x; // Get the thread ID

    // Calculate the problem ID
    size_t problemID = blockIdx.x; 

    // Create a 32-bit pointer to the constant memory so that I can do 32-bit operations to load the query bitmask
    uint32_t * gpuQuerydbArray32 = (uint32_t *)gpuQuerydbArray;

    // Get the length of the db sequence
    uint32_t seqLen = gpuSequenceLength[problemID+1] - gpuSequenceLength[problemID];
    
    if(seqLen > CPU_OFFLOAD_SIZE_GT_EQ){
        return;
    }

    // Get the address of the db sequence
    size_t dbAddr = gpuSequenceLength[problemID];
    // Get the db sequence
    uint8_t* dbSeq = gpuDbArray + dbAddr;

    // Clear the finalSeeds array
    for(uint i=tid;i<T_WARPS_PER_BLOCK*MAX_UG_SEEDS;i+=blockDim.x){
        finalSeeds[i] = 0;
        lastHits[i] = 0;    
    }

    // Clear the diagInfo
    for(uint i=tid;i<T_DIAGS;i+=blockDim.x){
        diagInfo[i] = 0;
    }

    if(laneID == 0){
        warpTotals[warpID] = 0;
        moreSeedsFlag[warpID] = 0;
    }

    uint8_t a;
    uint8_t b;
    uint8_t c;

    int diag_offset = window;

    int seedHitCounter = 0;
    bool moreSeeds = false;

    uint32_t seedOffset = 0;

    uint32_t * seedPointer = &finalSeeds[warpID*MAX_UG_SEEDS];

    uint32_t seed = 0;
    int seedlast = 0;

    for(int j=0; j<seqLen-2; j++){

        // Synchronize threads in the block, ensures that all writes to the diagInfo are visible to all threads and also
        __syncthreads();

        // Account for the seeds that are found in the last iteration
        // Using a ballot to get the mask of all threads in the warp that have a seed
        unsigned valid_seed_mask = __ballot_sync(0xFFFFFFFF, (seed != 0));
        
        if(valid_seed_mask){
            int pos = __popc(valid_seed_mask & ((1 << laneID) - 1));
            uint32_t threadIndex = pos + seedOffset;
            if(seed != 0 && threadIndex < T_WARPS_PER_BLOCK*MAX_UG_SEEDS){
                seedPointer[threadIndex] = seed;
                lastHits[threadIndex] = seedlast;
            }
            if(threadIndex >= T_WARPS_PER_BLOCK*MAX_UG_SEEDS){
                moreSeeds = true;
            }

            seedOffset += __popc(valid_seed_mask);
        }

        // Read the sequence from dbSequence
        a = dbSeq[j];
        b = dbSeq[j+1];
        c = dbSeq[j+2];

        // Calculate the index
        uint index = ((a * 784) + (b * 28) + (c)); // uint index = ((a * 28*28) + (b * 28) + (c));

        // Convert the index to the byte and bit index in the gpuQuerybitmaskArray
        uint byteIndex = index/BYTE_SIZE;
        uint bitIndex = index%BYTE_SIZE;

        // Get the byte from constant memory
        uint8_t bitsetByte = query_bitmask[byteIndex];

        // Check if the bit is set then copy the 32 bytes from the gpuQuerydbArray to the dbArray
        size_t queryAddr = index*T_BITSET_BYTES;

        seed = 0;

        if(bitsetByte & (1 << bitIndex)){

            uint16_t query_prefix_val = query_prefix[index];
            uint16_t query_lookup_val = query_prefix[index+1] - query_prefix_val;


            if(query_lookup_val > WARP_SIZE){
               moreSeeds = true; // If the query_lookup_val is greater than WARP_SIZE then we have more seeds than we can handle 
            }

            // If the laneID is less than the query_lookup_val then read the seedhit from the seed_lookup_table
            if(tid < query_lookup_val){
                uint16_t seedHit = seed_lookup_table[query_prefix_val + tid];
                //uint16_t seedHit = seed_holder[query_prefix_val + tid];
                uint16_t query_offset = seedHit;
                uint16_t subject_offset = j;

                // Calculate the diagonal number
                int diagNo = T_QUERY_BITSET_SIZE - 1 + ((int)subject_offset - (int)query_offset);

                // The lower 16 bits of diagInfo is the last hit
                volatile int lastHit = diagInfo[diagNo];

                // increment the seedHitCounter
                seedHitCounter++;

                int diff = subject_offset - lastHit + diag_offset; // Find the difference between the current hit and the last hit

                if(diff >= window){  
                    lastHit = subject_offset + diag_offset;
                    diagInfo[diagNo] = lastHit;
                    continue;
                }

                if(diff < wordsize){
                    continue;
                }        

                seed = (query_offset << 16) | subject_offset;
                seedlast = lastHit;

                seedCnt++;

                lastHit = subject_offset + diag_offset;
                diagInfo[diagNo] = lastHit;
    
            }
        }
        
    }

    // Using a ballot to get the mask of all threads in the warp that have a seed
    unsigned valid_seed_mask = __ballot_sync(0xFFFFFFFF, (seed != 0));
        
    if(valid_seed_mask){
        int pos = __popc(valid_seed_mask & ((1 << laneID) - 1));
        uint32_t threadIndex = pos + seedOffset;
        if(seed != 0 && threadIndex < T_WARPS_PER_BLOCK*MAX_UG_SEEDS){
            seedPointer[threadIndex] = seed;
            lastHits[threadIndex] = seedlast;
        }
        if(threadIndex >= T_WARPS_PER_BLOCK*MAX_UG_SEEDS){
            moreSeeds = true;
        }

        seedOffset += __popc(valid_seed_mask);
    }

    uint32_t ballot_more = __ballot_sync(0xFFFFFFFF, moreSeeds);   

    if(laneID == 0){
        warpTotals[warpID] = static_cast<uint8_t>(seedOffset);
        if(ballot_more != 0){
            moreSeedsFlag[warpID] = 1;
        }
    }

    __syncthreads(); 

    #pragma unroll
    for(uint i=0;i<T_WARPS_PER_BLOCK;i++){
        if(moreSeedsFlag[i]){
            if(warpID == 0 && laneID == 0){
                gpuNumSeeds[blockIdx.x] = 255;
            }
            return;
        }
    }

    int gappedExtension = 0;

    for(uint w=0; w < T_WARPS_PER_BLOCK; w++){

        uint8_t warpSeedCnt = warpTotals[w];

        uint warpResponsibility = ceil(warpSeedCnt/(float)T_WARPS_PER_BLOCK_FLOAT);

        for(uint wr = 0; wr < warpResponsibility; wr++){

            uint16_t addr = w*MAX_UG_SEEDS + wr*T_WARPS_PER_BLOCK + warpID;


            if(wr*T_WARPS_PER_BLOCK + warpID >= warpSeedCnt)
                break;
                    
            uint query_offset = (finalSeeds[addr] >> 16);
            uint subject_offset = (finalSeeds[addr] & 0xFFFF);
            int lastHit = lastHits[addr];
        
            int left_d = 0;
            int right_d = 0;
            int left_score = 0;
            int right_score = 0;
            int score = 0;

            #pragma unroll
            for(int j=0;j< wordsize; j++){
                
                score += score_matrix[gpuQuery[query_offset+j]][dbSeq[subject_offset+j]];
        
                // update the score
                if(score > left_score){
                    left_score = score;
                    right_d = j + 1;
                }
            }
        
            uint16_t q_off = query_offset + right_d - 1;
            uint16_t s_off = subject_offset + right_d -1;
        
            right_d = 0;
                
            int n = min((uint)s_off,(uint)q_off);
        
            s_BlastAaExtendLeft( (char *)dbSeq, (char *)gpuQuery, s_off, q_off, n, &left_score, &left_d, laneID, streamID+blockIdx.x);
        
            n = min((uint)(seqLen - s_off - 1),(uint)(gpuQueryLength - q_off - 1));

            if(left_d >= ( (s_off+1) - (lastHit+wordsize))){
                s_BlastAaExtendRight( (char *)dbSeq, (char *)gpuQuery, s_off, q_off, n, left_score , &right_score, &right_d, laneID, streamID+blockIdx.x);
            }

            if(right_score <= left_score){
                right_score = left_score;
                right_d = 0;
            }

            score = max(left_score, right_score);

            uint32_t gappedSeed = 0;

            if(score >= ungapped_cutoff){
                gappedExtension = 1; 

                int hsp_q = q_off + 1 - left_d;
                int hsp_s = s_off + 1 - left_d;
                int hsp_len = left_d + right_d;

                int q_start = hsp_q;
                int q_length = hsp_len;
                int s_start = hsp_s;
                int s_length = hsp_len;
                int max_offset = 0;

                volatile int score = 0;

                char * s = (char *) dbSeq;
                char * q = (char *) gpuQuery;
                
                if (q_length <= HSP_MAX_WINDOW) {
                    max_offset = q_start + q_length/2;
                    s_length = subject_offset + max_offset - query_offset + 1;
                    q_length = max_offset + 1;
                }else{

                    int hsp_end = q_start + HSP_MAX_WINDOW;
                    score = 0;
                    for(int index=laneID;index<HSP_MAX_WINDOW;index+=WARP_SIZE){
                        score = score_matrix[q[q_start+index]][s[s_start+index]];
                    }
                    score = __reduce_add_sync(0xFFFFFFFF, score);

                    int max_score = score;
                    int max_offset = hsp_end - 1;

                    hsp_end = q_start + min(q_length, s_length);

                    for(int j=0;j<((hsp_end-(q_start+HSP_MAX_WINDOW)));j++){
                        score -= score_matrix[q[q_start+j]][s[s_start+j]];
                        score += score_matrix[q[q_start+HSP_MAX_WINDOW+j]][s[s_start+HSP_MAX_WINDOW+j]];
                        if(score > max_score){
                            max_score = score;
                            max_offset = q_start + HSP_MAX_WINDOW + j;
                        }
                    }

                    if(max_score > 0){
                        max_offset -= HSP_MAX_WINDOW/2;
                    }else{
                        max_offset = q_start;
                    }
                    s_length = subject_offset + max_offset - query_offset + 1;
                    q_length = max_offset + 1;
                }

                gappedSeed = (q_length << 16) | s_length;
            }

            if(laneID == 0){
                finalSeeds[addr] = gappedSeed;
            }
        }
        
    }

    // sync the threads in the block
    __syncthreads();

    // Using lane0 of warp0 to write the gapped seeds to the gappedSeeds array
    if(warpID == 0 && laneID == 0){

        int seedCount = 0;

        for(int i=0;i<T_WARPS_PER_BLOCK;i++){
            // Get the warp seed count
            int warpSeedCount = warpTotals[i];
            for(int j=0;j<warpSeedCount;j++){
                if(finalSeeds[i*MAX_UG_SEEDS + j] != 0){
                    if(seedCount < MAX_SEEDS_GAPPED){
                        gappedSeeds[blockIdx.x*MAX_SEEDS_GAPPED + seedCount] = finalSeeds[i*MAX_UG_SEEDS + j];
                    }
                    seedCount++;
                }
            }
        }

        // If the laneZeroSeedCnt is greater than the MAX_SEEDS_GAPPED then set the gpuNumSeeds to 2
        if(seedCount > MAX_SEEDS_GAPPED){
            gpuNumSeeds[blockIdx.x] = 255;
        }
        else{
            if(seedCount > 0)
                gpuNumSeeds[blockIdx.x] = seedCount;
        }

    }

}

template<unsigned int T_QUERY_BITSET_SIZE, unsigned int T_SUBJECT_LEN, unsigned int T_THREADS_PER_PROBLEM>
__global__ void gpuBLAZE_everything(uint8_t * gpuDbArray, uint32_t * gpuSequenceLength, uint8_t * gpuQuerydbArray, char * gpuQuery, size_t gpuQueryLength, uint32_t streamID, uint8_t * gpuNumSeeds, int ungapped_cutoff, int gapped_cutoff, uint32_t * gappedSeeds, uint32_t * gpuHSPLen)
{
    // Calculate the constants during compile time
    constexpr float T_WARPS_PER_BLOCK_FLOAT = static_cast<float>(T_THREADS_PER_PROBLEM) / WARP_SIZE;
    constexpr int T_WARPS_PER_BLOCK = T_THREADS_PER_PROBLEM / WARP_SIZE;
    constexpr int T_BITSET_BYTES = (T_QUERY_BITSET_SIZE + 7) / 8;
    constexpr int T_DIAGS = T_QUERY_BITSET_SIZE + T_SUBJECT_LEN - 1;
    constexpr int T_BUFFER_SIZE = (T_QUERY_BITSET_SIZE < T_SUBJECT_LEN) ? T_QUERY_BITSET_SIZE : T_SUBJECT_LEN;

    // Return if the gpuNumSeeds is 255 or 0
    if(gpuNumSeeds[blockIdx.x] == 255 || gpuNumSeeds[blockIdx.x] == 0){
        return;
    }

    __shared__ uint8_t warpTotals[T_WARPS_PER_BLOCK]; // Holds the total number of seeds found by each warp.

    //get the warp id
    int warpID = threadIdx.x / WARP_SIZE;

    // Get the laneID
    uint8_t laneID = threadIdx.x & 0x1F; // Gets the lane ID.

    // Calculate the problem ID
    size_t problemID = blockIdx.x;

    // Get the length of the db sequence
    uint32_t seqLen = gpuSequenceLength[problemID+1] - gpuSequenceLength[problemID];

    if(seqLen > CPU_OFFLOAD_SIZE_GT_EQ){
        return;
    }

    // Get the address of the db sequence
    size_t dbAddr = gpuSequenceLength[problemID];
    // Get the db sequence
    uint8_t* dbSeq = gpuDbArray + dbAddr;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////// SW_OWN/////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    // allocate shared memory for the scoring row of size MAX_QRY_LEN and type BlastGapDP
    __shared__ BlastGapDP score_array[T_WARPS_PER_BLOCK][T_BUFFER_SIZE+1];

    if(laneID == 0){
        warpTotals[warpID] = 0;
    }

    // Get the total number of problems and find the warp responsibility for each warp
    uint8_t warpSeedCnt = gpuNumSeeds[blockIdx.x];

    uint warpResponsibility = ceil(warpSeedCnt/(float)T_WARPS_PER_BLOCK_FLOAT);

        // Warp Seed Offset
    uint32_t warpSeedOffset = warpID*warpResponsibility;

    int gappedExtension = 0;

    char * s = (char *) dbSeq;
    char * q = (char *) gpuQuery;

    for(uint wr = 0; wr < warpResponsibility; wr++){

        // Loop and check if any of the warpTotals are greater than 0, if so then break out of the loop
        for(int i=0;i<T_WARPS_PER_BLOCK;i++){
            if(warpTotals[i] > 0){
                break;
            }
        }

        uint32_t seedIndex = warpSeedOffset + wr;

        if(seedIndex >= warpSeedCnt)
            break;

        // Get the seed from the gappedSeeds array
        uint32_t seed = gappedSeeds[blockIdx.x*MAX_SEEDS_GAPPED + seedIndex];

        int q_length = seed >> 16;
        int s_length = seed & 0xFFFF;
        int ungapped_score = 0;

        int gapped_right = 0;
        int gapped_left = 0;

        gapped_left = semiGappedKernelLeft(q, s, q_length, s_length, &score_array[warpID][0], 1);
        if (q_length < gpuQueryLength && s_length < seqLen) {
            gapped_right = semiGappedKernelRight(q+q_length-1, s+s_length-1, gpuQueryLength-q_length, seqLen-s_length, &score_array[warpID][0], 0);
        }

        int totalScore = gapped_left + gapped_right;

        if(totalScore >= gapped_cutoff){
            gappedExtension = 1;
            if(laneID == 0){
                warpTotals[warpID] = 1;
            }
            break;
        }

        __syncwarp();       

    }

    // sync the threads in the block
    __syncthreads();

    int warpGappedExtension = __ballot_sync(0xFFFFFFFF, gappedExtension);
    
    if(laneID == 0){
        if(warpGappedExtension != 0)
            warpTotals[warpID] = 1;
        else
            warpTotals[warpID] = 0;
    }

    // Synchronize the threads in the block
    __syncthreads();

    if(warpID == 0 and laneID == 0){

        int tmp = 0;
        for(int i=0;i<T_WARPS_PER_BLOCK;i++){
            if(warpTotals[i] != 0){
                tmp = 1;
                break;
            }
        }

        gpuNumSeeds[blockIdx.x] = tmp;
    }

    __syncthreads();

}

#define CHECK_CUDA_ERROR(val) check((val), #val, __FILE__, __LINE__)
void check(cudaError_t err, const char* const func, const char* const file,
           const int line)
{
    if (err != cudaSuccess)
    {
        std::cerr << "CUDA Runtime Error at: " << file << ":" << line
                  << std::endl;
        std::cerr << cudaGetErrorString(err) << " " << func << std::endl;
        std::exit(EXIT_FAILURE);
    }
}


#define CHECK_LAST_CUDA_ERROR() checkLast(__FILE__, __LINE__)
void checkLast(const char* const file, const int line)
{
    cudaError_t const err{cudaGetLastError()};
    if (err != cudaSuccess)
    {
        std::cerr << "CUDA Runtime Error at: " << file << ":" << line
                  << std::endl;
        std::cerr << cudaGetErrorString(err) << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

extern "C" {

    int BLAZE_GPU(BlastSeqSrc* seq_src, uint32_t threadIndex, uint32_t totalThreads, uint32_t queryLen, uint32_t * survivingSeqs, uint32_t * num_survivors, int ungapped_cutoff, int gapped_cutoff){

        // Measure the time it takes to execute the kernel using c++ chrono library
        auto start = std::chrono::high_resolution_clock::now();

        BlastSeqSrcGetSeqArg seq_arg;
        memset((void *)&seq_arg, 0, sizeof(seq_arg));
        seq_arg.encoding = eBlastEncodingProtein;

        uint32_t thread_id = threadIndex;

        // Get DBWorkInfo
        uint64_t * indexArray;
        uint64_t * prefixIndexArray;

        uint64_t indexArraySize;
        uint64_t *chunkIndexArray;
        uint64_t num_chunks;

        uint64_t *chunkStart;
        uint64_t *chunkEnd;

        uint64_t * dbOffset;
        char * dbHostData;

        char * gpuDBAlloc;

        char * queryDBPointer;
        char * queryPointer;


        uint32_t seqsProcessed = 0;

        BlastSeqSrcGetDBPointer(seq_src, &gpuDBAlloc);
        BlastSeqSrcGetQueryPointers(seq_src, &queryDBPointer, &queryPointer);
        BlastSeqSrcGetDBWorkInfoGPU(seq_src, &prefixIndexArray, &indexArraySize, &num_chunks, &chunkStart, &chunkEnd, &dbOffset);

        // Get the Isolation query pointers.
        uint16_t *gpuQueryLookup;
        uint16_t *gpuQueryOffsets;
        uint16_t *gpuSeedTable;

        BlastSeqSrcGetIsolationPointers(seq_src, &gpuQueryLookup, &gpuQueryOffsets, &gpuSeedTable);

        // Get the indexArray, the 0th entry 
        indexArray = &prefixIndexArray[1];

        BlastSeqSrcGetHostDBPointer(seq_src, &dbHostData);

        size_t totalNumSeqs = BlastSeqSrcGetNumSeqs(seq_src);

        uint32_t seqCounter = 0;

        unsigned qSize = queryLen;

        qSize =  ((qSize + 31) / 32) * 32;

        // So at this point you are going to replace the value of j with the thread_id and are going to find the maximum required bytes per thread that is going to be allocated.
        uint64_t maxBytes = 0;
        size_t max_num_oids = 0;

        for(int i = 0; i < num_chunks; i++){
                
            uint64_t oidStart = chunkStart[i*totalThreads+thread_id];
            uint64_t oidEnd = chunkEnd[i*totalThreads+thread_id];

            // Calculate the number of bytes required 
            uint64_t requiredBytes = indexArray[oidEnd] - indexArray[oidStart - 1];

            // Update the maxBytes
            maxBytes = MAX(maxBytes, requiredBytes);

            // Update the max_num_oids
            max_num_oids = MAX(max_num_oids, oidEnd - oidStart + 1);
        }

        cudaEvent_t events[NUM_BUFFERS];
        uint32_t * sequenceLengths[NUM_BUFFERS]; 
        uint32_t * gpuSequenceLengths[NUM_BUFFERS];
        uint8_t *bitmaskArray[num_chunks];
        uint8_t * gpuBitmaskArray[NUM_BUFFERS];
        uint8_t *threadSequenceData[NUM_BUFFERS];
        uint8_t *threadSequenceDataGPU[NUM_BUFFERS];
        uint32_t * gpuSeeds[NUM_BUFFERS];
        uint32_t * gpuHSPLen[NUM_BUFFERS];

        uint32_t * sequenceLengthsPtr; 
        uint32_t * gpuSequenceLengthsPtr;
        uint8_t *bitmaskArrayPtr;
        uint8_t * gpuBitmaskArrayPtr;
        uint8_t *threadSequenceDataPtr;
        uint8_t *threadSequenceDataGPUPtr;
        uint32_t * gpuSeedsPtr;
        uint32_t * gpuHSPLenPtr;


        maxBytes = dbOffset[1] - dbOffset[0]; 

        // CudaAllocHost the bitmask array wuing the max_num_oids for the total number of chunks
        // Allocate the bitmask array
        for(uint c = 0; c < num_chunks; c++){
            // Allocate the bitmask array
            CHECK_CUDA_ERROR(cudaHostAlloc(&bitmaskArray[c], (max_num_oids) * sizeof(uint8_t),  cudaHostAllocDefault));
        }


        for(uint b = 0; b < NUM_BUFFERS; b++){
            CHECK_CUDA_ERROR(cudaHostAlloc(&sequenceLengths[b], (max_num_oids + 1) * sizeof(uint32_t),  cudaHostAllocDefault));
            CHECK_CUDA_ERROR(cudaMalloc(&gpuSequenceLengths[b], (max_num_oids + 1) * sizeof(uint32_t)));
            CHECK_CUDA_ERROR(cudaMalloc(&gpuBitmaskArray[b], (max_num_oids) * sizeof(uint8_t)));
            CHECK_CUDA_ERROR(cudaHostAlloc(&threadSequenceData[b], maxBytes * sizeof(uint8_t),  cudaHostAllocDefault));
            CHECK_CUDA_ERROR(cudaMalloc(&threadSequenceDataGPU[b], maxBytes * sizeof(uint8_t)));        
            CHECK_CUDA_ERROR(cudaEventCreate(&events[b]));
            CHECK_CUDA_ERROR(cudaMalloc(&gpuSeeds[b], (max_num_oids) * MAX_SEEDS_GAPPED * sizeof(uint32_t)));
            CHECK_CUDA_ERROR(cudaMalloc(&gpuHSPLen[b], (max_num_oids) * MAX_SEEDS_GAPPED * sizeof(uint32_t)));

            // CLear the gpuBitmaskArray
            CHECK_CUDA_ERROR(cudaMemset(gpuBitmaskArray[b], 0, (max_num_oids) * sizeof(uint8_t)));
            // Clear the gpuSeeds
            CHECK_CUDA_ERROR(cudaMemset(gpuSeeds[b], 0, (max_num_oids) * MAX_SEEDS_GAPPED * sizeof(uint32_t)));
            // Clear the gpuHSPLen
            CHECK_CUDA_ERROR(cudaMemset(gpuHSPLen[b], 0, (max_num_oids) * MAX_SEEDS_GAPPED * sizeof(uint32_t)));

        }

        // Create a cuda stream.
        cudaStream_t stream;
        CHECK_CUDA_ERROR(cudaStreamCreate(&stream));

        char * survivingSequences;
        BlastSeqSrcGetSequenceBitmaskPointer(seq_src, &survivingSequences);

        for(int i = 0; i < num_chunks; i++){
            
            uint64_t oidStart = chunkStart[i*totalThreads+thread_id];
            uint64_t oidEnd = chunkEnd[i*totalThreads+thread_id];

            uint32_t subject_seq_idx = 0;

            unsigned BufferID = i % NUM_BUFFERS;

            bitmaskArrayPtr = bitmaskArray[i]; 
            threadSequenceDataPtr = threadSequenceData[BufferID];
            gpuSequenceLengthsPtr = gpuSequenceLengths[BufferID];
            gpuBitmaskArrayPtr = gpuBitmaskArray[BufferID];
            sequenceLengthsPtr = sequenceLengths[BufferID];
            threadSequenceDataGPUPtr = threadSequenceDataGPU[BufferID];
            gpuSeedsPtr = gpuSeeds[BufferID];
            gpuHSPLenPtr = gpuHSPLen[BufferID];

            size_t seq_length = 0;

            int OIDs_processed = oidEnd - oidStart + 1;

            // Synchronize the event
            if(i >= NUM_BUFFERS){
                CHECK_CUDA_ERROR(cudaEventSynchronize(events[BufferID]));
            }

            for(int oid_search = oidStart; oid_search <= oidEnd; oid_search++){

                seqsProcessed++;

                // Update the seq_arg.oid
                seq_arg.oid = oid_search;

                // Get the sequence
                BlastSeqSrcGetSequence(seq_src, &seq_arg);
                
                // Get the length of the sequence
                seq_length = seq_arg.seq->length;

                // Copy the sequence to the threadSequenceData
                memcpy(threadSequenceDataPtr + subject_seq_idx, seq_arg.seq->sequence, seq_length);

                // Copy the sequence length to the sequenceLengths array
                sequenceLengthsPtr[oid_search - oidStart] = subject_seq_idx;

                // Release the sequence
                BlastSeqSrcReleaseSequence(seq_src, &seq_arg);

                // Update the subject_seq_idx
                subject_seq_idx += seq_length;
            }

            //Add the last index to the sequenceLengths array
            sequenceLengthsPtr[oidEnd - oidStart + 1] = subject_seq_idx;


            CHECK_CUDA_ERROR(cudaMemcpyAsync(threadSequenceDataGPUPtr , threadSequenceDataPtr, maxBytes * sizeof(uint8_t), cudaMemcpyHostToDevice, stream));
            CHECK_CUDA_ERROR(cudaMemcpyAsync(gpuSequenceLengthsPtr, sequenceLengthsPtr, (OIDs_processed + 1) * sizeof(uint32_t), cudaMemcpyHostToDevice, stream));

            // Get the query length
            int sLen = 256;

            sLen = (seq_length <= 64)   ? 64   :
                (seq_length <= 128)  ? 128  :
                (seq_length <= 192)  ? 192  :
                (seq_length <= 256)  ? 256  :
                (seq_length <= 320)  ? 320  :
                (seq_length <= 384)  ? 384  :
                (seq_length <= 448)  ? 448  :
                (seq_length <= 512)  ? 512  :
                (seq_length <= 576)  ? 576  :
                (seq_length <= 1024) ? 1024 :
                (seq_length <= 2048) ? 2048 :
                                            2048;
            

            auto it = kernelMap.find({qSize, sLen});
            if (it != kernelMap.end()) {
                dim3 block(T_PER_PROB), grid(OIDs_processed);
                it->second<<<grid, block, 0, stream>>>(threadSequenceDataGPUPtr, gpuSequenceLengthsPtr, (uint8_t*)queryDBPointer, queryPointer, queryLen, oidStart, gpuBitmaskArrayPtr, ungapped_cutoff, gapped_cutoff, gpuSeedsPtr, gpuHSPLenPtr, gpuQueryLookup, gpuQueryOffsets, gpuSeedTable);
            } else {
                printf("No matching kernel for Q=%d, S=%d.\n", qSize, sLen);
            }

            auto itE = kernelMapEverything.find({qSize, sLen});
            if (itE != kernelMapEverything.end()) {
                dim3 block(T_PER_PROB), grid(OIDs_processed);
                itE->second<<<grid, block, 0, stream>>>(threadSequenceDataGPUPtr, gpuSequenceLengthsPtr, (uint8_t*)queryDBPointer, queryPointer, queryLen, oidStart, gpuBitmaskArrayPtr, ungapped_cutoff, gapped_cutoff, gpuSeedsPtr, gpuHSPLenPtr);
            } else {
                printf("No matching kernel for Q=%d, S=%d.\n", qSize, sLen);
            }

            CHECK_CUDA_ERROR(cudaMemcpyAsync(bitmaskArrayPtr, gpuBitmaskArrayPtr, OIDs_processed * sizeof(uint8_t), cudaMemcpyDeviceToHost, stream));

            // Record the event after the stream operations
            CHECK_CUDA_ERROR(cudaEventRecord(events[BufferID], stream));

            // Asynchronously clear the gpuBitmaskArray
            CHECK_CUDA_ERROR(cudaMemsetAsync(gpuBitmaskArrayPtr, 0, OIDs_processed * sizeof(uint8_t), stream));

            // Asynchronously clear the gpuSeeds
            CHECK_CUDA_ERROR(cudaMemsetAsync(gpuSeedsPtr, 0, OIDs_processed * MAX_SEEDS_GAPPED * sizeof(uint32_t), stream));

            // Asynchronously clear the gpuHSPLen
            CHECK_CUDA_ERROR(cudaMemsetAsync(gpuHSPLenPtr, 0, OIDs_processed * MAX_SEEDS_GAPPED * sizeof(uint32_t), stream));
        }

        // // Make sure that all the events are recorded
        for(uint b = 0; b < NUM_BUFFERS; b++){
            CHECK_CUDA_ERROR(cudaEventSynchronize(events[b]));
        }

        for(int i = 0; i < num_chunks; i++){
            
            uint64_t oidStart = chunkStart[i*totalThreads+thread_id];
            uint64_t oidEnd = chunkEnd[i*totalThreads+thread_id];

            bitmaskArrayPtr = bitmaskArray[i]; 

            for(int k=0; k < oidEnd - oidStart + 1; k++){
                // Read each byte of the bitmaskArray
                if(bitmaskArrayPtr[k] != 0){
                    survivingSeqs[seqCounter] = oidStart + k;
                    seqCounter++;
                }
            }

        }

        // Update the number of survivors
        *num_survivors = seqCounter;

        // Stop the timer
        auto end = std::chrono::high_resolution_clock::now();

        // Calculate the duration
        std::chrono::duration<double> duration = end - start;

        // // Print the seqs processed and the execution time
        // fprintf(stdout, "Seqs Processed: %d, Execution Time: %f seconds\n", seqsProcessed, duration.count());

        return 0;

    }

    void copySurvivals(uint8_t * survivals_cpu, uint8_t *gpuSurvivals, size_t OIDs_processed, void* stream){

        // Get the stream
        cudaStream_t* cudaStream = (cudaStream_t*)stream;  

        // Copy the numSeeds from the GPU to the CPU
        cudaMemcpy(survivals_cpu, gpuSurvivals, OIDs_processed * sizeof(uint8_t), cudaMemcpyDeviceToHost);

    }

    void createCudaStream(void** streamPtr) {
        cudaStream_t* stream = (cudaStream_t*)malloc(sizeof(cudaStream_t));
        cudaStreamCreate(stream);
        *streamPtr = (void*)stream;
    }

    void deleteCudaStream(void* stream) {
        cudaStream_t* cudaStream = (cudaStream_t*)stream;
        cudaStreamDestroy(*cudaStream);
        free(cudaStream);
    }


    void hostAlloc(void** pointer, size_t size) {
        cudaHostAlloc(pointer, size, cudaHostAllocDefault);
        // CHECK FOR ERRORS
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess) {
            std::cerr << "CUDA error in hostAlloc: " << cudaGetErrorString(error) << std::endl;
        }
    }

    void deviceAlloc(void** pointer, size_t size) {
        cudaMalloc(pointer, size);

        // CHECK FOR ERRORS
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess) {
            std::cerr << "CUDA error in deviceAlloc: " << cudaGetErrorString(error) << std::endl;
        }
    }

    void deviceMemCpyToDevice(void*dst, void* src, size_t size){
        cudaMemcpy(dst, src, size, cudaMemcpyHostToDevice);

        // CHECK FOR ERRORS
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess) {
            std::cerr << "CUDA error in deviceMemCpyToDevice: " << cudaGetErrorString(error) << std::endl;
        }
    }

    void deviceMemCpyToDeviceAsync(void *dst, void *src, size_t size, void *stream){
        // Get the stream
        cudaStream_t* cudaStream = (cudaStream_t*)stream;  

        // Perform an asynchronous copy
        cudaMemcpyAsync(dst, src, size, cudaMemcpyHostToDevice, *cudaStream);

        // CHECK FOR ERRORS
        cudaError_t error = cudaGetLastError();

        if (error != cudaSuccess) {
            std::cerr << "CUDA error in deviceMemCpyToDeviceAsync: " << cudaGetErrorString(error) << std::endl;
        }

    }

    void streamSynchronize(void *stream){
        // Get the stream
        cudaStream_t* cudaStream = (cudaStream_t*)stream;  

        // Synchronize the stream
        CHECK_CUDA_ERROR(cudaStreamSynchronize(*cudaStream));

        // CHECK FOR ERRORS
        cudaError_t error = cudaGetLastError();

        if (error != cudaSuccess) {
            std::cerr << "CUDA error in streamSynchronize: " << cudaGetErrorString(error) << std::endl;
        }
    }

    void hostFree(void* pointer) {
        cudaFreeHost(pointer);

        // CHECK FOR ERRORS
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess) {
            std::cerr << "CUDA error in hostFree: " << cudaGetErrorString(error) << std::endl;
        }
    }

    void deviceFree(void* pointer) {
        cudaFree(pointer);

        // CHECK FOR ERRORS
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess) {
            std::cerr << "CUDA error in deviceFree: " << cudaGetErrorString(error) << std::endl;
        }
    }


    size_t getTotalGPUMemory() {
        size_t free, total;
        cudaMemGetInfo(&free, &total);
        return total;
    }

    size_t getFreeGPUMemory() {
        size_t free, total;
        cudaMemGetInfo(&free, &total);
        return free;
    }

 
    void copyConstantMemory(uint8_t * gpuQuerybitmaskArray, size_t querybitmaskarray_gpu_size){
        // Copy the bitmaskarraay to the constant memory
        cudaMemcpyToSymbol( query_bitmask, gpuQuerybitmaskArray , querybitmaskarray_gpu_size*sizeof(uint8_t), 0, cudaMemcpyHostToDevice );

        // Check for errors
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess) {
            std::cerr << "CUDA error in copyConstantMemory1: " << cudaGetErrorString(error) << std::endl;
        }

        // copy the score matrix to the constant memory
        cudaMemcpyToSymbol( score_matrix, matrix, (28*28)*sizeof(int), 0, cudaMemcpyHostToDevice );

        // Check for errors
        error = cudaGetLastError();
        if (error != cudaSuccess) {
            std::cerr << "CUDA error in copyConstantMemory2: " << cudaGetErrorString(error) << std::endl;
        }
    
    }

    void copyConstantMemoryLUT(uint16_t * seedLUT, size_t seedLUTSize){
        // Copy the bitmaskarraay to the constant memory
        cudaMemcpyToSymbol( seed_holder, seedLUT , seedLUTSize*sizeof(uint16_t), 0, cudaMemcpyHostToDevice );

        // Check for errors
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess) {
            std::cerr << "CUDA error in copyConstantMemory3: " << cudaGetErrorString(error) << std::endl;
        }
    
    }

}