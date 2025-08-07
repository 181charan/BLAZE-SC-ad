/*
Base code setup done, have a bunch of debug code too!
*/

#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/execution_policy.h>

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

// Include the header to use sleep
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


//#define DEBUG_PROBLEM 196686927

#define DEBUG_PROBLEM 241450368
#define S_OFF 397
#define Q_OFF 190

#define QUERY_BITSET_SIZE 256
#define SUBJECT_LEN 256
#define THREADS_PER_PROBLEM 64 // This is the number of threads that will be assigned to each problem
#define MAX_SEEDS 8 // This is the max number of seeds that I can perform ungapped extension on. You need to get this number after removing the filterting from ungapped extension.
#define MAX_SEEDS_GAPPED 8 // This is the max number of seeds that I can perform gapped extension on. You need to get this number after removing the filterting from gapped extension.
#define WARPS_PER_BLOCK_FLOAT (static_cast<float>(THREADS_PER_PROBLEM) / static_cast<float>(WARP_SIZE)) //#define WARPS_PER_BLOCK_FLOAT 2.0
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
// __constant__ uint16_t seed_holder[27648];

__device__ inline unsigned generateMask(int laneID) {
    // If laneID is 31, return all 1s (32 bits)
    if (laneID == 31) {
        return 0xFFFFFFFF;
    }
    // Otherwise, set bits up to laneID
    return (1U << (laneID + 1)) - 1U;
}

__inline__ __device__
unsigned restrictKernel(const char *A, const char *B, int M, int N, BlastGapDP * score_array, bool reverse_sequence){

    // Return if either of the sequence is greater than 256
    // if(M > 256 || N > 256){
    //     return 0;
    // }


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


    // // Print the sequences using thread thread zero of warp zero
    // if(tid == 0 ){
    //     printf("A(%d): ", sizeA);
    //     for(int i=0;i<=sizeA;i++){
    //         printf("%d,", (int)(A_block[i]));
    //     }
    //     printf("\n");

    //     printf("B(%d): ", sizeB);
    //     for(int i=0;i<=sizeB;i++){
    //         printf("%d,", (int)(B_block[i]));
    //     }
    //     printf("\n");
    // }

    int sweepcnt = ceil(((double)maxAB/(double)WARP_SIZE));

    // Print the sweep count
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

    //TODO: Can definitely parallelize this one!
    for (int i = 1; i <= minAB; i++) {

        if (score < -X_DROPOFF){
            // break;
            score_array[i].best = MININT;
            score_array[i].best_gap = MININT;
        }else{
            score_array[i].best = score;
            score_array[i].best_gap = score - GAP_OPEN_EXTEND;
            score -= GAP_EXTEND;
        }
    }

    __syncwarp();

    int b_gap = 0;

    int best_prev = 0;
    int bestgap_prev = 0;

    for(int e=0;e<sweepcnt;e++){

        // Reset the score and the score_gap_row
        score = MININT;
        score_gap_row = MININT;

        // Get the best score from thread 31 and assign it to all
        best_score = __shfl_sync(0xffffffff, best_score, 31);

        for(int x=(e*WARP_SIZE); x < (e*WARP_SIZE + WARP_SIZE + minAB); x++){
        
            int i = e*WARP_SIZE + laneID;
            int j = x - i;

            temp1 = __shfl_up_sync(0xffffffff, best_prev, 1, 32);
            temp2 = __shfl_up_sync(0xffffffff, bestgap_prev, 1, 32);

            int tmp_best_score = __shfl_up_sync(0xffffffff, best_score, 1, 32);

            if(laneID != 0){
                best_score = max(tmp_best_score, best_score); // only update the bes score from previous threads if it is better than what you found!
            }

            if(!(i < 0 || ( (i + 1) > maxAB) || j < 0 || (j + 1) > minAB)){

                int i_temp = i + 1;
                int b_index = j;

                if(laneID == 0){
                    next_score = score_array[b_index].best;
                    score_gap_col = score_array[b_index].best_gap;
                }
                else{
                    next_score = temp1;
                    score_gap_col = temp2;
                }

                if(reverse_sequence){
                    next_score += score_matrix[A_block[maxAB - i_temp]][B_block[ minAB - (b_index+1)]];

                }else{
                    next_score += score_matrix[A_block[i_temp]][B_block[b_index+1]];
                }

                //next_score += score_matrix[A_block[i_temp]][B_block[b_index+1]];

                if (i_temp % RESTRICT_SIZE != 0) {

                    if (b_index % RESTRICT_SIZE != 0) {

                        /* the majority of cases fall here; a gap
                        may not start in either A or B */

                        if (best_score - score > X_DROPOFF) {
                            best_prev = MININT;
                        }
                        else {

                            if (score > best_score) {
                                best_score = score;
                            }
                            best_prev = score;
                        }
                    }
                    else {
                        b_gap += RESTRICT_SIZE;

                        if (score < score_gap_col)
                            score = score_gap_col;

                        if (best_score - score > X_DROPOFF) {
                            best_prev = MININT;
                        }
                        else {
                            if (score > best_score) {
                                best_score = score;
                            }

                            score_gap_col -= GAP_EXTEND;
                            bestgap_prev = max(score - GAP_OPEN_EXTEND, score_gap_col);
                            best_prev = score;

                        }
                    }



                }
                else{

                    //if (b_index != b_gap) {
                    if (b_index % RESTRICT_SIZE != 0) {

                        /* gap may not start in B. Compute
                            the resulting two-term recurrence */

                        if (score < score_gap_row)
                            score = score_gap_row;

                        if (best_score - score > X_DROPOFF) {
                            best_prev = MININT;
                        }
                        else {

                            if (score > best_score) {
                                best_score = score;
                            }

                            score_gap_row -= GAP_EXTEND;
                            score_gap_row = max(score - GAP_OPEN_EXTEND, score_gap_row);
                            best_prev = score;
                        }
                    }
                    else {

                        /* the ordinary case: a gap may start in A
                            or B and the full three-term recurrence
                            must be computed. This happens once every
                            RESTRICT_SIZE*RESTRICT_SIZE cells */

                        b_gap += RESTRICT_SIZE;
                        //score_gap_col = score_array[b_index].best_gap;

                        if (score < score_gap_col)
                            score = score_gap_col;
                        if (score < score_gap_row)
                            score = score_gap_row;

                        if (best_score - score > X_DROPOFF) {
                            best_prev = MININT;
                        }
                        else {
                            if (score > best_score) {
                                best_score = score;
                            }

                            score_gap_row -= GAP_EXTEND;
                            score_gap_col -= GAP_EXTEND;
                            bestgap_prev = max(score - GAP_OPEN_EXTEND, score_gap_col);
                            score_gap_row = max(score - GAP_OPEN_EXTEND,
                                                score_gap_row);
                            best_prev = score;
                        }

                    }

                }

                score = next_score;

                if(laneID == 31){
                    score_array[b_index].best = best_prev;
                    score_array[b_index].best_gap = bestgap_prev;
                }

            }


        }

        
    }

    // holds the global maximum across all threads
    globalMax = __reduce_max_sync( 0xFFFFFFFF, best_score);

    return globalMax;
    
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

    // Print the sweep count
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

    // Print the sweep count
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

    uint loopcnt = ceil(((double)(n+1)/(double)WARP_SIZE)); // fixed a corner case bug here, as the loop for left extension starts at i=n and does down till i>=0, so the loopcount should be n+1 and not n.

    // max score is definitely zero, however the best_i is initialized to n+1
    int globalMax = 0;
    int globalMaxD = n+1;
    int personalMax = 0;

    int prevScore = 0;

    for(int l=0;l<loopcnt;l++){
        
            int j = n - laneID - l*WARP_SIZE;

            bool valid = j >= 0;

            // Get a mask of all threads which are valid.
            unsigned valid_mask = __ballot_sync(0xFFFFFFFF, valid);

            // Calcscore
            int localScore = valid ? score_matrix[q[j]][s[j]] : MININT;

            for (int offset = 1; offset < warpSize; offset *= 2) { // This essentailly calculates the prefix sum by staggering the computation into log number of lops.
                int value = __shfl_up_sync(0xFFFFFFFF, localScore, offset);
                if (laneID >= offset) {
                    localScore += value;
                }
            }

            localScore += prevScore;

            // // For the debug problems print the local score and the j and the laneID in this format "score: 11, matrix[q[i]][s[i]]: 11, i: 114"
            // if(problemID == DEBUG_PROBLEM && s_off == S_OFF && q_off == Q_OFF){
            //     printf("score: %d, matrix[q[i]][s[i]]: %d, i: %d\n", localScore, score_matrix[q[j]][s[j]], j);
            // }

            // Print the calc score score matrix vlaue and the j and the laneID in this format "score: 11, matrix[q[i]][s[i]]: 11, i: 114"
            //printf("Localscore: %d, matrix[q[i]][s[i]]: %d, i: %d\n", localScore, (valid ? score_matrix[qt[j]][st[j]] : MININT), j);
            
            if(laneID == 0)
            personalMax = max(personalMax, localScore);  // initial value for each thread
            else
                personalMax =  localScore;// this is the max score that is broadcast to everyone in the warp.
            
            for (int offset = 1; offset < warpSize; offset *= 2) {
                // Use personalMax here so that the running maximum is updated
                int candidate = __shfl_up_sync(0xFFFFFFFF, personalMax, offset);
                if (laneID >= offset) {
                    personalMax = max(personalMax, candidate);
                }
            }

            // Print the personal max, personal max D and j
            //printf("personalMax: %d, personalMaxD: %d, j: %d\n", personalMax, personalMaxD, j);

            // The current score is the same as personalMax
            unsigned drop_mask = __ballot_sync(valid_mask, ((personalMax - localScore) >= dropoff) ); // This is the dropflag array from the "cuBLASTP Paper"
            unsigned non_zero_dm = __ffs(drop_mask) - 1;

            if(drop_mask != 0){
                
                int theAbsoluteMax = __shfl_sync(0xFFFFFFFF, personalMax, non_zero_dm); // this is the score that is broarcast to everyone in the warp.
                
                // Check if the absolute max is different than the global Max,
                if(theAbsoluteMax == globalMax){
                    // The previous max score is still valid so get all threads to update their local variables and break out of the thread.
                    left_d = globalMaxD + 1;
                    left_score = globalMax;

                    break;
                }else{

                    unsigned tmask = generateMask(non_zero_dm); // This is the mask of all threads less than or equal to the calling thread.
                    unsigned smask = __ballot_sync(0xFFFFFFFF, localScore == theAbsoluteMax); // Get the mask of threads, which have scores equal to the glibal max.
                    smask = smask & tmask; // This is the mask of all threads less than or equal to the calling thread that have a score equal to the absolute max.
                    uint8_t tempID = __ffs(smask) - 1; // Find the first non-zero bit in the mask // You get the ID of the first thread for which this is true.
                    left_score = theAbsoluteMax; // This is the score that is broadcast to all threads in the warp.
                    left_d = l*WARP_SIZE + tempID + 1;
                    break;
                }
        
            }else{

                int thread31_max = __shfl_sync(0xFFFFFFFF, personalMax, 31); // this is the max score that is broadcast to all threads in the warp.
                if(thread31_max != globalMax){
                    // Do not update the global max score and the index.
                    unsigned smask = __ballot_sync(0xFFFFFFFF, localScore == thread31_max);
                    uint8_t maxID = __ffs(smask) - 1;
                    // Update the global max score and the index.
                    globalMax = thread31_max;
                    globalMaxD = l*WARP_SIZE + maxID; // This is the index of the first thread in the warp with the same score as thread 31.
                    // If this is the last iteration then you just return the global max score and the index.
                } //else when they are equal do not have to update the max value

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

    // Broadcast the left_d and the left_score from the first thread of each warp to all the threads of the warp
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

    uint loopcnt = ceil(((double)(n+1)/(double)WARP_SIZE)); // fixed a corner case bug here, as the loop for left extension starts at i=n and does down till i>=0, so the loopcount should be n+1 and not n.

    // max score is definitely zero, however the best_i is initialized to n+1
    int globalMax = score;
    int globalMaxD = -1;

    int currentLocalMax = score;

    for(int l=0;l<loopcnt;l++){
        
        int j = l*WARP_SIZE + laneID;

        //if(problemID == DEBUG_PROBLEM && s_off == S_OFF && q_off == Q_OFF){
            //printf("personalMax: %d, personalMaxD: %d, j: %d\n", personalMax, personalMaxD, j);
        
            bool valid = j < n;

            // Get a mask of all threads which are valid.
            unsigned valid_mask = __ballot_sync(0xFFFFFFFF, valid);

            // Calcscore
            int localScore = valid ? score_matrix[qt[j]][st[j]] : MININT;
            
            for (int offset = 1; offset < warpSize; offset *= 2) { // This essentailly calculates the prefix sum by staggering the computation into log number of lops.
                int value = __shfl_up_sync(0xFFFFFFFF, localScore, offset);
                if (laneID >= offset) {
                    localScore += value;
                }
            }

            localScore += prevScore;

            // Print the calc score score matrix vlaue and the j and the laneID in this format "score: 11, matrix[q[i]][s[i]]: 11, i: 114"
            //printf("Localscore: %d, matrix[q[i]][s[i]]: %d, i: %d\n", localScore, (valid ? score_matrix[qt[j]][st[j]] : MININT), j);
            
            if(laneID == 0)
            personalMax = max(personalMax, localScore);  // initial value for each thread
            else
                personalMax =  localScore;// this is the max score that is broadcast to everyone in the warp.
            
            for (int offset = 1; offset < warpSize; offset *= 2) {
                // Use personalMax here so that the running maximum is updated
                int candidate = __shfl_up_sync(0xFFFFFFFF, personalMax, offset);
                if (laneID >= offset) {
                    personalMax = max(personalMax, candidate);
                }
            }

            // Print the personal max, personal max D and j
            //printf("personalMax: %d, personalMaxD: %d, j: %d\n", personalMax, personalMaxD, j);

            // The current score is the same as personalMax
            unsigned drop_mask = __ballot_sync(valid_mask, ((localScore <= 0) || (personalMax - localScore) >= dropoff) ); // This is the dropflag array from the "cuBLASTP Paper"
            unsigned non_zero_dm = __ffs(drop_mask) - 1;

            if(drop_mask != 0){
                // This means that some thread has called terminate. So you have to break out of the loop
                // However the only think you need to terminate is the previouis best score and its index position.
                // Broadcast the currentLocalMax of the termination calling thread to all threads in the warp.
                // This is the max score, check if this score matches with previous max score, if it does then you can report the previous max score and its index,
                // Else check which scores in the current iteration match the currentLocalMax of the terminating thread and report its offset to all the threads in the warp.
                // USe the thread 0 to update the final pointer variable.
                
                // Step 1. Find the index of the first thread that called it quits.
                // Print non_zero_dm and j and laneID
                //printf("dropoff was %d and maxVal %d and non_zero_dm: %d, j: %d, laneID: %d\n", dropoff, globalMax, non_zero_dm, j, laneID);
                
                // Broadcast the currentLocalMax from the thread that called it quits to all threads in the warp.
                int theAbsoluteMax = __shfl_sync(0xFFFFFFFF, personalMax, non_zero_dm); // this is the score that is broarcast to everyone in the warp.
                // Print the absolute max
                //printf("theAbsoluteMax: %d and laneID: %d\n", theAbsoluteMax, laneID);
                
                // Check if the absolute max is different than the global Max,
                if(theAbsoluteMax == globalMax){
                    // The previous max score is still valid so get all threads to update their local variables and break out of the thread.
                    right_d = globalMaxD + 1;
                    right_score = globalMax;
                    
                    // Print that the score was the old score so you do not have to update it.
                    //printf("The score was the old score so you do not have to update it\n");
                    break;
                }else{

                    // Print that the score was the new score so you have to update it.
                    //printf("The score was the new score so you have to update it\n");

                    // Get a mask of all threads lesser than or equal to the thread calling it quits that have a score equal to the absolute max.
                    // Generate a mask of all threads less than or equal to the calling thread.
                    unsigned tmask = generateMask(non_zero_dm); // This is the mask of all threads less than or equal to the calling thread.
                    // Using this mask get the first thread that has a score equal to the absolute max.
                    unsigned smask = __ballot_sync(0xFFFFFFFF, localScore == theAbsoluteMax); // Get the mask of threads, which have scores equal to the glibal max.
                    // And and the two masks together
                    smask = smask & tmask; // This is the mask of all threads less than or equal to the calling thread that have a score equal to the absolute max.
                    // Get the first non-zero bit in the mask
                    uint8_t tempID = __ffs(smask) - 1; // Find the first non-zero bit in the mask // You get the ID of the first thread for which this is true.

                    right_score = theAbsoluteMax; // This is the score that is broadcast to all threads in the warp.
                    right_d = l*WARP_SIZE + tempID + 1;

                    // Print the tempID, the absolute max score and the right_d
                    //printf("tempID: %d, theAbsoluteMax: %d and right_d: %d\n", tempID, theAbsoluteMax, right_d);

                    break;
                }
        
            }else{
                // SO in this iteration the drop was not called. So you will have to get the max score from thread 31 and broadcast it to all threads.
                // If this is equal to the previous max score then you do not update the max score and the index, 
                // Else check which is the first thread in the warp with the same score as the thread 31's amx score and update the global max for all threads and the index for all threads.
                // There can be anotehr corner case situation where this might be the last iteration, in that case the last thread might nit be thread 31. In that case you are going to get the Index of the last
                // Active thread in the warp and not the thread 31. So you have to check if the last thread is the thread 31 or not.
                // Get the maxscore from thread 31
                int thread31_max = __shfl_sync(0xFFFFFFFF, personalMax, 31); // this is the max score that is broadcast to all threads in the warp.
    
                // Check if the max score is equal to the global max score, if it is then you do not update the global max score and the index.
                if(thread31_max != globalMax){
                    // Do not update the global max score and the index.
                    unsigned smask = __ballot_sync(0xFFFFFFFF, localScore == thread31_max);
                    uint8_t maxID = __ffs(smask) - 1;
                    // Update the global max score and the index.
                    globalMax = thread31_max;
                    globalMaxD = l*WARP_SIZE + maxID; // This is the index of the first thread in the warp with the same score as thread 31.
                    // If this is the last iteration then you just return the global max score and the index.
                    

                } //else when they are equal do not have to update the max value

                if(l == loopcnt-1){
                    right_d = globalMaxD + 1;
                    right_score = globalMax;
                    break;
                }

                // But you still have to update the previous score and the personal max score.
                //prevScore = __shfl_sync(active, calcScore, (31 - __clz(active))); // this is the cumulative SUM
                // Broadcast the calcScore from the 31st thread to all other threads, and if the last thread is not 31
                // you are in the last iteration anyway, so whatever you do here is not going to matter anyways.
                prevScore = __shfl_sync(0xFFFFFFFF, localScore, (warpSize - 1)); // this is the cumulative SUM
                personalMax = __shfl_sync(0xFFFFFFFF, personalMax, (warpSize - 1)); // this is the cumulative SUM
                
                // Print the prevScore
                //printf("prevScore: %d and laneID: %d\n", prevScore, laneID);

            }    
        
        //}

    //     if(j >= n){

    //         // Maybe here I have to check for the active mask, if all the threads are there only then you update it, else someone else will do it for you!
    //         right_d = globalMaxD + 1;
    //         right_score = globalMax;

    //         if(right_score < 0)
    //             right_score = 0;

    //         break;
    //     }

    //     int calcScore = score_matrix[qt[j]][st[j]]; // This gets the current score for the protein alphabet

    //     int active = __activemask();        // Gets the mask of threads that are active
        
    //     for (int offset = 1; offset < warpSize; offset *= 2) { // This essentailly calculates the prefix sum by staggering the computation into log number of lops.
    //         int value = __shfl_up_sync(active, calcScore, offset);
    //         if (laneID >= offset) {
    //             calcScore += value;
    //         }
    //     }

    //     calcScore += prevScore;

    //     int maxVal = __reduce_max_sync(active, calcScore); // Gets the maximum score, which includes the cumulative sum.

    //     globalMax = max(globalMax, maxVal); // check if the max you just found is greater than the current max, this score is retained across the various threads as technically the max value is and it is broadcast across the entire warp.

    //     // Broadcast the calcScore from last active thread to all other threads
    //     prevScore = __shfl_sync(active, calcScore, (31 - __clz(active))); // this is the cumulative SUM

    //     unsigned mask = __ballot_sync(active, calcScore == globalMax); // Get the mask of threads, which have scores equal to the glibal max.

    //     uint8_t maxID = __ffs(mask) - 1; // Find the first non-zero bit in the mask // You get the ID of the first thread for which this is true.

    //     if(maxVal == globalMax){
    //         globalMaxD = l*WARP_SIZE + maxID; // Even this is this guaranteed that everyone has the same score???? Makes sense as everyone calculates the same maxID.
    //     }

    //     mask = __ballot_sync(active, laneID > maxID); // Find all the threads that have a higher thread ID than the thread with the max value, this is very important. Do I really need to check ONLY after the maxID? Essentially I need this to be zero. TODO: Might need to modify this. This is because thread zero is going to be skipped even if the maxID is zero. 
    //     // So there should be some form of checking to make sure that all threads are participaint in this dropoff thresholding in the circumstance in which you do not update the max score? Maybe not, maybe you do not have to? But I will fix this first and only after I KNOW that I have everything working is when I will try to remove the probably redundant checks.
        
    //     // [TODO:] One plaguing question I have is, in the situation where none of them are equal, Essentially in the situation in which you do not improve the best score you have seen in a particular iteration.
    //     // I strongly belive that this should be fixed. I belive that I am lucking out in this situation. 
    //     // This should address that situation!
    //     if(maxID == 255){
    //         mask = active;
    //     }

    //     uint32_t activeMask = mask;

    //     //mask = __ballot_sync(mask, ((calcScore <= 0) || (calcScore - maxVal) < dropoff) ); // This is the dropflag array from the "cuBLASTP Paper"

    //     mask = __ballot_sync(mask, ((calcScore <= 0) || (globalMax - calcScore) >= dropoff) ); // This is the dropflag array from the "cuBLASTP Paper"

    //     mask = mask & activeMask;

    //     if (mask != 0) {

    //         right_d = globalMaxD + 1;
    //         right_score = globalMax;

    //         // If the problem number is the debug problem then print the query_offset, subject_offset, left_d, left_score.
    //         if(right_score < 0)
    //             right_score = 0;

    //         break; // Break out of the loop if you meet the dropoff condition  
    //     }

    //     // If you are on the last iteration of the loop that means that you reached the end of the loop and you still havent met the dropoff condition
    //     // So report the right_d as the l*WARP_SIZE + (31 - __clz(mask)) and the right_score as the maxVal
    //     if(l == loopcnt-1){

    //         right_d = globalMaxD + 1;
    //         right_score = globalMax;

    //         if(right_score < 0)
    //             right_score = 0;

    //     }

    }

    __syncwarp();

    *right_d_o = __shfl_sync(0xFFFFFFFF, right_d, 0);
    *right_score_o = __shfl_sync(0xFFFFFFFF, right_score, 0); 

}

// This is the function for right extension
__forceinline__ __device__ int o_BlastAaExtendRight( int bmatrix[][28], uint8_t * subject, char * query, uint32_t subject_length , uint32_t query_length , uint32_t s_off, uint32_t q_off, int bdropoff, int * blen, int maxscore, uint32_t * s_last_off, uint tid)
{

    int i, n, best_i = -1;
    int score = maxscore;

    n = min(subject_length - s_off, query_length - q_off);

    #ifdef UG_DEBUG
    int searchSpace = n;
    #endif

    char * s = (char *) subject;
    char * q = (char *) query;

    s += s_off;
    q += q_off;

    for (i = 0; i < n; i++) {
        score += bmatrix[q[i]][s[i]];

        if (score > maxscore) {
            maxscore = score;
            best_i = i;
        }

        /* The comparison below is really >= and is different than the old
           code (e.g., blast.c:BlastWordExtend_prelim). In the old code the
           loop continued as long as sum > X (X being negative).  The loop
           control here is different and we *break out* when the if statement 
           below is true. */
        if (score <= 0 || (maxscore - score) >= bdropoff){

            #ifdef UG_DEBUG
            searchSpace = i;
            #endif

            break;
        }
        
    }

    // only if #UG_DEBUG is defined print the following
    #ifdef UG_DEBUG
    if(tid == probDebug)
        printf("ugextr:%d,%d(%d)\n",best_i + 1, searchSpace + 1, tid);
    #endif

    *blen = best_i + 1;
    *s_last_off = s_off + i;
    return maxscore;


}

// This is the function for left extension
__forceinline__ __device__ int o_BlastAaExtendLeft( int matrix[][28], uint8_t * subject, char * query, uint32_t s_off, uint32_t q_off, int bdropoff, int * length, int maxscore, uint tid){

        int i, n, best_i;
        int score = maxscore;
        
        n = min((uint)s_off,(uint)q_off);
        best_i = n + 1;

        #ifdef UG_DEBUG
        int searchSpace = 0;
        #endif

        char * s = (char *) subject;
        char * q = (char *) query;

        s += s_off - n;
        q += q_off - n;
        
        for (i = n; i >= 0; i--){

            score += matrix[q[i]][s[i]];

            if (score > maxscore){
                maxscore = score;
                best_i = i;
            }
            /* The comparison below is really >= and is different than the old
           code (e.g., blast.c:BlastWordExtend_prelim). In the old code the
           loop continued as long as sum > X (X being negative).  The loop
           control here is different and we *break out* when the if statement 
           below is true. */
           if ((maxscore - score) >= bdropoff)
           {
            #ifdef UG_DEBUG
            searchSpace = i; // [SCG]
            #endif
            
            break;
           }
            
        }

        #ifdef UG_DEBUG
        if(tid == probDebug)
            printf("ugextl:%d,%d(%d)\n", n - best_i + 1, n - searchSpace + 1, tid);
        #endif

        *length = n - best_i + 1;
        return maxscore;

}


// This is the two hit extension function
__forceinline__ __device__ int s_BlastAaExtendTwoHit(int matrix[][28], uint8_t * subject, char * query, uint32_t subject_length, uint32_t query_length , uint32_t s_left_off, uint32_t s_right_off, uint32_t q_right_off, int bdropoff, bool * right_extend, uint32_t * s_last_off, uint tid, int * right_d, int * left_d){

    // These are char because maximum length I can have is 255, will have to change this!
    //int left_d = 0;
    *left_d = 0;
    *right_d = 0;
    //score = 0;
    //int right_d = 0;
    int left_score = 0;
    int right_score = 0;
    int score = 0;
    int i = 0;

    for(i=0;i< wordsize; i++){
        // Calculate the score
        score += matrix[query[q_right_off+i]][subject[s_right_off+i]];

        // update the score
        if(score > left_score){
            left_score = score;
            *right_d = i + 1;
        }

    }

    q_right_off += *right_d;
    s_right_off += *right_d;

    *right_d = 0;
    *s_last_off = s_right_off;

    left_score = o_BlastAaExtendLeft(matrix, subject, query, s_right_off - 1, q_right_off - 1, bdropoff, left_d, 0, tid);

    /* Extend to the right only if left extension reached the first hit. */
    if (*left_d >= (s_right_off - s_left_off)) {
        *right_extend = true;
        right_score = o_BlastAaExtendRight(matrix, subject, query, subject_length, query_length, s_right_off, q_right_off, bdropoff, right_d, left_score, s_last_off, tid);
    }

    return (max(left_score,right_score));

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

    //get the warp id
    int warpID = threadIdx.x / WARP_SIZE;

    // Get the laneID
    uint8_t laneID = threadIdx.x & 0x1F; // Gets the lane ID.

    int tid = threadIdx.x; // Get the thread ID

    // Calculate the problem ID
    size_t problemID = blockIdx.x; // Currently I launch all the problems in a single block, hence the problemID is the same as the blockID

    // Create a 32-bit pointer to the constant memory so that I can do 32-bit operations to load the query bitmask
    uint32_t * gpuQuerydbArray32 = (uint32_t *)gpuQuerydbArray;

    // Get the length of the db sequence
    uint32_t seqLen = gpuSequenceLength[problemID+1] - gpuSequenceLength[problemID];
    
    // Return if either of the sequence is greater than 256
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
        
        // If there are valid seeds then for the debug problem print the query and subject offset
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

        // YOu know that you do not have to go into the global memory if not!
        if(bitsetByte & (1 << bitIndex)){

            // Each thread looks up the index
            //uint16_t query_lookup_val = query_lookup[index];
            //uint16_t query_prefix_val = query_prefix[index];
            //uint16_t query_lookup_val = query_prefix[index+1] - query_prefix_val;

            uint16_t query_prefix_val = query_prefix[index];
            uint16_t query_lookup_val = query_prefix[index+1] - query_prefix_val;

            // uint16_t query_prefix_val = seed_holder[index];
            // uint16_t query_lookup_val = seed_holder[index+1] - query_prefix_val;



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


    // Account for the seeds that are found in the last iteration
    // Using a ballot to get the mask of all threads in the warp that have a seed
    unsigned valid_seed_mask = __ballot_sync(0xFFFFFFFF, (seed != 0));
        
    // If there are valid seeds then for the debug problem print the query and subject offset
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

    __syncthreads(); // Makes sure that all the seeds and the warpTotals are written to shared memory
    
    // if the more seeds flag is set for any of the warps then use warp zero thread zero to set the gpuNumSeeds to 1
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

    // // Each warp handles warpTotals[warpID]/WARP_SIZE seeds for each warp, so loop over the T_WARPS_PER_BLOCK
    for(uint w=0; w < T_WARPS_PER_BLOCK; w++){

        uint8_t warpSeedCnt = warpTotals[w];

        //////////////////////////////// IMP /////////////////////////////////////////////////////////////
        /////////////////////////////// HARDCOADED VALUE ////////////////////////////////////////////////
        uint warpResponsibility = ceil(warpSeedCnt/(float)T_WARPS_PER_BLOCK_FLOAT);
        ////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////

        // Loop over the warpResponsibility for each warp
        for(uint wr = 0; wr < warpResponsibility; wr++){

            //uint16_t addr = w*WARP_SIZE*MAX_SEEDS + wr*T_WARPS_PER_BLOCK + warpID;
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
                
                // Calculate the score current all the threads of the warp are doing this, is this good?
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
                right_score = left_score; // Not sure if I wanna hardcode this here, maybe there is a better way to do this?
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

                    // This is totally redundant work that is being done by each thread. You can use some of the techniques that cuBLASTP proposed to parallelize this across threads of the warp.
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

            // Write whatever is the gappedseed value is back to the finalSeeds array using lane 0
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


        // For all the warps in the block iterate over the final seeds array and if
        // the element is not zero then write it to the gappedSeeds array
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
    size_t problemID = blockIdx.x; // Currently I launch all the problems in a single block, hence the problemID is the same as the blockID

    // Get the length of the db sequence
    uint32_t seqLen = gpuSequenceLength[problemID+1] - gpuSequenceLength[problemID];

    // Return if either of the sequence is greater than 256
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
        //int hsp_len = gpuHSPLen[blockIdx.x*MAX_SEEDS_GAPPED + seedIndex];
        int ungapped_score = 0;

        int gapped_right = 0;
        int gapped_left = 0;
        // bool isRestricted = false;
        // int restricted_cutoff = kRestrictedMult * gapped_cutoff;

        // if(ungapped_score < restricted_cutoff){
        //     isRestricted = true;
        // }

        // // Restricted
        // if( isRestricted ){
        //     gapped_left = restrictKernel(q, s, q_length, s_length, &score_array[warpID][0], 1);
        //     if (q_length < gpuQueryLength && s_length < seqLen) {
        //         gapped_right = restrictKernel(q+q_length-1, s+s_length-1, gpuQueryLength-q_length, seqLen-s_length, &score_array[warpID][0], 0);
        //     }

        // }else{

        gapped_left = semiGappedKernelLeft(q, s, q_length, s_length, &score_array[warpID][0], 1);
        if (q_length < gpuQueryLength && s_length < seqLen) {
            gapped_right = semiGappedKernelRight(q+q_length-1, s+s_length-1, gpuQueryLength-q_length, seqLen-s_length, &score_array[warpID][0], 0);
        }

        //}

        int totalScore = gapped_left + gapped_right;

        // // This results in some ambiguity so you have to mak this for recomputation!
        // if(isRestricted &&
        // totalScore < gapped_cutoff &&
        // totalScore >= restricted_cutoff){

        //     gapped_right = 0;
        //     gapped_left = semiGappedKernel(q, s, q_length, s_length, &score_array[warpID][0], 1);
        //     if (q_length < gpuQueryLength && s_length < seqLen) {
        //         gapped_right = semiGappedKernel(q+q_length-1, s+s_length-1, gpuQueryLength-q_length, seqLen-s_length, &score_array[warpID][0], 0);
        //     }
        //     totalScore = gapped_left + gapped_right;
        //     //gappedExtension = 1;
        // }

        if(totalScore >= gapped_cutoff){
            gappedExtension = 1;
            // using lane0  update the warpTotals to 1 if the gappedExtension is found
            if(laneID == 0){
                warpTotals[warpID] = 1;
            }
            break;

            //break;
        }

        // // fOR THE DEBUG PROBLEM PRINT THE gapped_left and gapped_right score and the q_length and s_length
        // if(streamID+blockIdx.x == DEBUG_PROBLEM && laneID == 0){
        //     printf("gapped_left: %d, gapped_right: %d, q_length: %d, s_length: %d and gapped extension %d\n", gapped_left, gapped_right, q_length, s_length, gappedExtension);
        // }

        __syncwarp();       

    }

    // sync the threads in the block
    __syncthreads();

    int warpGappedExtension = __ballot_sync(0xFFFFFFFF, gappedExtension);
    
    // If the warp has a gapped extension then set then get lane 0 to set warpTotals[warpID] to 1 alse set it to 0
    if(laneID == 0){
        if(warpGappedExtension != 0)
            warpTotals[warpID] = 1;
        else
            warpTotals[warpID] = 0;
    }

    // Synchronize the threads in the block
    __syncthreads();

    // use warpID 0 to check if either of the warps have a gapped extension if so then set gpuNumSeeds[blockIdx.x] = 1; else set it to 0
    if(warpID == 0 and laneID == 0){

        // loop over the warps and check if any of the warps have a gapped extension
        int tmp = 0;
        for(int i=0;i<T_WARPS_PER_BLOCK;i++){
            if(warpTotals[i] != 0){
                tmp = 1;
                break;
            }
        }

        // if(streamID+blockIdx.x == DEBUG_PROBLEM && laneID == 0){
        //     // Print that you found a gapped extesion for this seed iuncluding the debug problem and tmp values
        //     printf("Found a gapped extension for this seed, streamID: %d, blockIdx.x: %d, tmp: %d\n", streamID+blockIdx.x, blockIdx.x, tmp);
        //     // Print the warp gapped extension and warp totals
        //     for(int i=0;i<T_WARPS_PER_BLOCK;i++){
        //         printf("warpGappedExtension: %d, warpTotals: %d\n", warpGappedExtension, warpTotals[i]);
        //     }
        // }

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
        // We don't exit when we encounter CUDA errors in this example.
        // std::exit(EXIT_FAILURE);
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
        // We don't exit when we encounter CUDA errors in this example.
        // std::exit(EXIT_FAILURE);
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


        // MaxBytes is actually constant so replacing, but clean up TODO: [SCG]
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
            //CHECK_CUDA_ERROR(cudaHostAlloc(&bitmaskArray[b], (max_num_oids) * sizeof(uint8_t),  cudaHostAllocDefault));
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

            // Overwite the address of the bitmask array with the pointer of the surviving sequences
            //bitmaskArrayPtr = reinterpret_cast<uint8_t*>(&survivingSequences[oidStart]);
            bitmaskArrayPtr = bitmaskArray[i]; // I am going to use the bitmask array for the current chunk
            threadSequenceDataPtr = threadSequenceData[BufferID];
            gpuSequenceLengthsPtr = gpuSequenceLengths[BufferID];
            gpuBitmaskArrayPtr = gpuBitmaskArray[BufferID];
            sequenceLengthsPtr = sequenceLengths[BufferID];
            threadSequenceDataGPUPtr = threadSequenceDataGPU[BufferID];
            gpuSeedsPtr = gpuSeeds[BufferID];
            gpuHSPLenPtr = gpuHSPLen[BufferID];

            size_t seq_length = 0;

            // Get the number of OIDs processed
            int OIDs_processed = oidEnd - oidStart + 1;

            // CHECK_CUDA_ERROR(cudaStreamSynchronize(stream));

            // // //If BufferID is 0, then synchronize the event
            // if(BufferID == 0 && i != 0){
            //     CHECK_CUDA_ERROR(cudaEventSynchronize(events[BufferID]));
            // }
            // As soon as the i is greater than the number of buffers, you are in a situation where you have to reuse buffers
            // So make sure that the event on the buffer is synchronized before you use it.
        
            // Synchronize the event
            if(i >= NUM_BUFFERS){
                CHECK_CUDA_ERROR(cudaEventSynchronize(events[BufferID]));
            }

            for(int oid_search = oidStart; oid_search <= oidEnd; oid_search++){

                seqsProcessed++;

                // Update the seq_arg.oid
                seq_arg.oid = oid_search;

                // if chunk == 0 and thread_id == 0 and oid_search == 0
                if(i == 0 && thread_id == 0 && oid_search == oidStart){
                    seq_arg.oid = DEBUG_PROBLEM;
                }

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
            
            //sLen = 1024;

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

            //bitmaskArrayPtr = reinterpret_cast<uint8_t*>(&survivingSequences[oidStart]);
            bitmaskArrayPtr = bitmaskArray[i]; // I am going to use the bitmask array for the current chunk

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

        // Print the seqs processed and the execution time
        fprintf(stdout, "Seqs Processed: %d, Execution Time: %f seconds\n", seqsProcessed, duration.count());

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
        //CHECK_CUDA_ERROR(cudaStreamSynchronize(*cudaStream));
        cudaStreamSynchronize(*cudaStream);

        // // CHECK FOR ERRORS
        // cudaError_t error = cudaGetLastError();

        // if (error != cudaSuccess) {
        //     std::cerr << "CUDA error in streamSynchronize: " << cudaGetErrorString(error) << std::endl;
        // }
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
        // cudaMemcpyToSymbol( seed_holder, seedLUT , seedLUTSize*sizeof(uint16_t), 0, cudaMemcpyHostToDevice );

        // Check for errors
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess) {
            std::cerr << "CUDA error in copyConstantMemory3: " << cudaGetErrorString(error) << std::endl;
        }
    
    }

}