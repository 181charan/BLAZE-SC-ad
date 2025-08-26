/*  $Id: blast_seqsrc.c 604796 2020-04-02 12:32:43Z fongah2 $
 * ===========================================================================
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
 *
 */

/** @file blast_seqsrc.c
 * Definition of ADT to retrieve sequences for the BLAST engine and
 * low level details of the implementation of the BlastSeqSrc framework.
 */

#include <algo/blast/core/blast_seqsrc.h>
#include <algo/blast/core/blast_seqsrc_impl.h>
#include <limits.h>


/** Complete type definition of Blast Sequence Source ADT.
 * The members of this structure should only be accessed by BlastSeqSrc
 * implementations using the _BlastSeqSrcImpl_* functions.
 */
struct BlastSeqSrc {

    uint32_t * survivorsPointer;
    uint32_t survivorsCount;
    int isGPU;

    // Creating test pointers
    int *test;
    int *test2;
    int num_threads;
    int curr_thread;

    int gpu_read_thread;

    int *inputA;
    int *inputB;
    int *result;

    uint8_t *survivingSeedsPtr;

    char *DBPointer;
    FILE *fullDBPointer;

    uint64_t *indexPointer;
    uint64_t sizeOfPointerArray;
    uint64_t *startChunkArray;
    uint64_t *endChunkArray;
    uint64_t sizeOfChunkArray;

    uint64_t *indexPointerGPU;
    uint64_t sizeOfPointerArrayGPU;
    uint64_t *startChunkArrayGPU;
    uint64_t *endChunkArrayGPU;
    uint64_t sizeOfChunkArrayGPU;
    uint64_t *dbOffsetGPU;

    uint16_t* qlookupTable;
    uint16_t* qlookupPrefix;
    uint16_t* actualLookupTable;

    uint64_t *dbOffset;
    uint8_t *hostDBPtr;

    uint8_t *gpuQueryDB;
    uint8_t *gpuQueryBitmask;

    BlastSeqSrcConstructor NewFnPtr;       /**< Constructor */
    BlastSeqSrcDestructor  DeleteFnPtr;    /**< Destructor */
    BlastSeqSrcCopier      CopyFnPtr;      /**< Copier */

   /* Functions to set the number of threads in MT mode */
    SetInt4FnPtr      SetNumberOfThreads;  /**< Set number of threads */

   /* Functions to get information about database as a whole */
    GetInt4FnPtr      GetNumSeqs;     /**< Get number of sequences in set */
    GetInt4FnPtr      GetNumSeqsStats; /**< Number of sequences for statistical purposes. */
    GetInt4FnPtr      GetMaxSeqLen;   /**< Get length of longest seq in set */
    GetInt4FnPtr      GetMinSeqLen;   /**< Get length of longest seq in set */
    GetInt4FnPtr      GetAvgSeqLen;   /**< Get average length of sequences in 
                                         the set */
    GetInt8FnPtr      GetTotLen;      /**< Get tot length of all seqs in set */
    GetInt8FnPtr      GetTotLenStats; /**< Total length of all seqs for statistical purposes. */
    GetStrFnPtr       GetName;        /**< Get the name of the database */
    GetBoolFnPtr      GetIsProt;      /**< Find if database is a protein or 
                                         nucleotide */

   /* Functions that supports partial sequence fetching */
    GetBoolFnPtr      GetSupportsPartialFetching; /**< Find if database supports partial fetching */
    SetSeqRangeFnPtr  SetSeqRange;    /**< Setting ranges for partial fetching */

   /* Functions that deal with individual sequences */
    GetSeqBlkFnPtr    GetSequence;    /**< Retrieve individual sequence */
    GetInt4FnPtr      GetSeqLen;      /**< Retrieve given sequence length */
    ReleaseSeqBlkFnPtr ReleaseSequence; /**< Deallocate individual sequence 
                                         (if applicable) */

   /* Functions to iterate over sequences in the database */
    AdvanceIteratorFnPtr IterNext;    /**< Gets next oid from the iterator */

    ResetChunkIteratorFnPtr ResetChunkIterator; /**< Reset the implementation's
                                                  chunk "bookmark"
                                                  */
   
    void*             DataStructure;  /**< ADT holding the sequence data */

    char*             InitErrorStr;   /**< initialization error string */
#ifdef KAPPA_PRINT_DIAGNOSTICS
    GetGisFnPtr       GetGis;         /**< Retrieve a sequence's gi(s) */
#endif /* KAPPA_PRINT_DIAGNOSTICS */
};

BlastSeqSrc* BlastSeqSrcNew(const BlastSeqSrcNewInfo* bssn_info)
{
    BlastSeqSrc* retval = NULL;

    if (!bssn_info) {
        return (BlastSeqSrc*) NULL;
    }

    if ( !(retval = (BlastSeqSrc*) calloc(1, sizeof(BlastSeqSrc)))) {
        return (BlastSeqSrc*) NULL;
    }

    /* Save the constructor and invoke it */
    if ((retval->NewFnPtr = bssn_info->constructor)) {
        retval = (*retval->NewFnPtr)(retval, bssn_info->ctor_argument);
    } else {
        sfree(retval);
    }

    return retval;
}

BlastSeqSrc* BlastSeqSrcFree(BlastSeqSrc* seq_src)
{
    BlastSeqSrcDestructor destructor_fnptr = NULL;
    BlastSeqSrc* retval;

    if (!seq_src) {
        return (BlastSeqSrc*) NULL;
    }

    if (seq_src->InitErrorStr) {
        sfree(seq_src->InitErrorStr);
    }

    /* This could leave a memory leak if destructor function pointer is not
     * initialized! It is the implementation's resposibility to provide this */
    if ( !(destructor_fnptr = (*seq_src->DeleteFnPtr))) {
        sfree(seq_src);
        return (BlastSeqSrc*) NULL;
    }

    retval = (BlastSeqSrc*) (*destructor_fnptr)(seq_src);
    ASSERT(retval == NULL);
    sfree(seq_src);
    return retval;
}

BlastSeqSrc* BlastSeqSrcCopy(const BlastSeqSrc* seq_src)
{
    BlastSeqSrcCopier copy_fnptr = NULL;
    BlastSeqSrc* retval = NULL;

    if (!seq_src) {
        return (BlastSeqSrc*) NULL;
    }

    if ( !(retval = (BlastSeqSrc*) BlastMemDup(seq_src, sizeof(BlastSeqSrc)))) {
        return (BlastSeqSrc*) NULL;
    }

    /* If copy function is not provided, just return a copy of the structure */
    if ( !(copy_fnptr = (*seq_src->CopyFnPtr))) {
        return retval;
    }

    return (BlastSeqSrc*) (*copy_fnptr)(retval);
}

char* BlastSeqSrcGetInitError(const BlastSeqSrc* seq_src)
{
    if ( !seq_src || !seq_src->InitErrorStr ) {
        return NULL;
    }

    return strdup(seq_src->InitErrorStr);
}

void BlastSeqSrcSetNumberOfThreads(BlastSeqSrc* seq_src, int n_threads)
{
    
    if (!seq_src || !seq_src->SetNumberOfThreads) {
        return;
    }

    (*seq_src->SetNumberOfThreads)(seq_src->DataStructure, n_threads);
}

void BlastSeqSrcSetThreads(BlastSeqSrc* seq_src, int n_threads)
{

    seq_src->num_threads = n_threads; 

}

Int4
BlastSeqSrcGetNumSeqs(const BlastSeqSrc* seq_src)
{
    ASSERT(seq_src);
    ASSERT(seq_src->GetNumSeqs);
    return (*seq_src->GetNumSeqs)(seq_src->DataStructure, NULL);
}

Int4
BlastSeqSrcGetNumSeqsStats(const BlastSeqSrc* seq_src)
{
    ASSERT(seq_src);
    ASSERT(seq_src->GetNumSeqsStats);
    return (*seq_src->GetNumSeqsStats)(seq_src->DataStructure, NULL);
}

Int4
BlastSeqSrcGetMaxSeqLen(const BlastSeqSrc* seq_src)
{
    ASSERT(seq_src);
    ASSERT(seq_src->GetMaxSeqLen);
    return (*seq_src->GetMaxSeqLen)(seq_src->DataStructure, NULL);
}

Int4
BlastSeqSrcGetMinSeqLen(const BlastSeqSrc* seq_src)
{
    ASSERT(seq_src);
    /* TODO this function may not be available for all seq_src */
    return (seq_src->GetMinSeqLen) ?
        (*seq_src->GetMinSeqLen)(seq_src->DataStructure, NULL) 
        : BLAST_SEQSRC_MINLENGTH;
}

Int4
BlastSeqSrcGetAvgSeqLen(const BlastSeqSrc* seq_src)
{
    ASSERT(seq_src);
    ASSERT(seq_src->GetAvgSeqLen);
    return (*seq_src->GetAvgSeqLen)(seq_src->DataStructure, NULL);
}

Int8
BlastSeqSrcGetTotLen(const BlastSeqSrc* seq_src)
{
    ASSERT(seq_src);
    ASSERT(seq_src->GetTotLen);
    return (*seq_src->GetTotLen)(seq_src->DataStructure, NULL);
}

Int8
BlastSeqSrcGetTotLenStats(const BlastSeqSrc* seq_src)
{
    ASSERT(seq_src);
    ASSERT(seq_src->GetTotLenStats);
    return (*seq_src->GetTotLenStats)(seq_src->DataStructure, NULL);
}

const char*
BlastSeqSrcGetName(const BlastSeqSrc* seq_src)
{
    ASSERT(seq_src);
    ASSERT(seq_src->GetName);
    return (*seq_src->GetName)(seq_src->DataStructure, NULL);
}

Boolean
BlastSeqSrcGetIsProt(const BlastSeqSrc* seq_src)
{
    ASSERT(seq_src);
    ASSERT(seq_src->GetIsProt);
    return (*seq_src->GetIsProt)(seq_src->DataStructure, NULL);
}

Boolean
BlastSeqSrcGetSupportsPartialFetching(const BlastSeqSrc* seq_src)
{
    ASSERT(seq_src);
    if (seq_src->GetSupportsPartialFetching) {
        return (*seq_src->GetSupportsPartialFetching)(seq_src->DataStructure, NULL);
    }
    return FALSE;
}

void
BlastSeqSrcSetSeqRanges(const BlastSeqSrc* seq_src,
                        BlastSeqSrcSetRangesArg* arg)
{
    ASSERT(seq_src);
    if (seq_src->SetSeqRange) {
        (*seq_src->SetSeqRange)(seq_src->DataStructure, arg);
    }
}

Int2
BlastSeqSrcGetSequence(const BlastSeqSrc* seq_src, 
                       BlastSeqSrcGetSeqArg* getseq_arg)
{
    ASSERT(seq_src);
    ASSERT(seq_src->GetSequence);
    ASSERT(getseq_arg);
    return (*seq_src->GetSequence)(seq_src->DataStructure, getseq_arg);
}

Int4
BlastSeqSrcGetSeqLen(const BlastSeqSrc* seq_src, void* oid)
{
    ASSERT(seq_src);
    ASSERT(seq_src->GetSeqLen);
    return (*seq_src->GetSeqLen)(seq_src->DataStructure, oid);
}

void
BlastSeqSrcReleaseSequence(const BlastSeqSrc* seq_src,
                           BlastSeqSrcGetSeqArg* getseq_arg)
{
    ASSERT(seq_src);
    ASSERT(seq_src->ReleaseSequence);
    ASSERT(getseq_arg);
    (*seq_src->ReleaseSequence)(seq_src->DataStructure, getseq_arg);
}

int BlastSeqSrcGetNumberOfThreads(BlastSeqSrc* seq_src)
{
    if (!seq_src) {
        return 0; // This is an error condition, ideally check for this
    }

    return seq_src->num_threads;
}

void BlastSeqSrcSetCurrentThread(BlastSeqSrc* seq_src, int thread)
{
    if (!seq_src) {
        return;
    }

    seq_src->curr_thread = thread;
}

void BlastSeqSrcSetThreadType(BlastSeqSrc* seq_src, int is_gpu){
   
    if (!seq_src) {
        return;
    }

    seq_src->gpu_read_thread = is_gpu;
}

int BlastSeqSrcGetIsGPU(BlastSeqSrc* seq_src){

    if (!seq_src) {
        return 0;
    }

    return seq_src->isGPU;
}

void BlastSeqSrcSetIsGPU(BlastSeqSrc* seq_src, int is_gpu){

    if (!seq_src) {
        return 0;
    }

    seq_src->isGPU = is_gpu;
}

void BlastSeqSrcSetSurvivors(BlastSeqSrc* seq_src, uint32_t* survivorsPointer, uint32_t survivorsCount){

    if (!seq_src) {
        return;
    }

    seq_src->survivorsPointer = survivorsPointer;
    seq_src->survivorsCount = survivorsCount;

}

void BlastSeqSrcGetSurvivors(BlastSeqSrc* seq_src, uint32_t** survivorsPointer, uint32_t* survivorsCount){

    if (!seq_src) {
        return;
    }

    *survivorsPointer = seq_src->survivorsPointer;
    *survivorsCount = seq_src->survivorsCount;
}



int BlastSeqSrcGetThreadType(BlastSeqSrc* seq_src){

    if (!seq_src) {
        return 0;
    }

    return seq_src->gpu_read_thread;
}

int BlastSeqSrcGetCurrentThread(BlastSeqSrc* seq_src)
{
    if (!seq_src) {
        return 0; // This is an error condition, ideally check for this
    }

    return seq_src->curr_thread;
}

void BlastSeqSrcSetTest(BlastSeqSrc* seq_src, int* test)
{
    seq_src->test = test;
}

void BlastSeqSrcSetTest2(BlastSeqSrc* seq_src, int* test)
{
    seq_src->test2 = test;
}

int BlastSeqSrcGetTest(const BlastSeqSrc* seq_src, unsigned indx)
{
    // Get the value using the index to the test pointer
    int *addr = seq_src->test;
    return addr[indx];
}

int BlastSeqSrcGetTest2(const BlastSeqSrc* seq_src, unsigned indx)
{
    int *addr = seq_src->test2;
    return addr[indx];
}

void BlastSeqSrcEnableTest(const BlastSeqSrc* seq_src, unsigned indx)
{
    // Get the value using the index to the test pointer
    int *addr = seq_src->test;
    addr[indx] = 1;

}

void BlastSeqSrcEnableTest2(const BlastSeqSrc* seq_src, unsigned indx)
{
    int *addr = seq_src->test2;
    addr[indx] = 1;
}

void BlastSeqSrcWait(const BlastSeqSrc* seq_src)
{
    // Get the value using the index to the test pointer
    int *addr = seq_src->test;
    // Wait for all threads to finish
    while (1) {
        int count = 0;
        for (int i = 0; i < seq_src->num_threads; i++) {
            if (addr[i] == 1) {
                count++;
            }
        }
        if (count == seq_src->num_threads) {
            break;
        }
    }

}


void BlastSeqSrcSetFullDBPointer(BlastSeqSrc* seq_src, FILE* db_pointer)
{
    seq_src->fullDBPointer = db_pointer;
}

void BlastSeqSrcGetFullDBPointer(const BlastSeqSrc* seq_src, FILE** db_pointer)
{
    *db_pointer = seq_src->fullDBPointer;
}


void BlastSeqSrcSetDBPointer(BlastSeqSrc* seq_src, char* db_pointer)
{
    seq_src->DBPointer = db_pointer;
}

void BlastSeqSrcGetDBPointer(const BlastSeqSrc* seq_src, char** db_pointer)
{
    *db_pointer = seq_src->DBPointer;
}

void BlastSeqSrcSetDBIndexPointer(BlastSeqSrc* seq_src, uint64_t* db_pointer)
{
    seq_src->indexPointer = db_pointer;
}

void BlastSeqSrcGetDBIndexPointer(const BlastSeqSrc* seq_src, uint64_t** db_pointer)
{
    *db_pointer = seq_src->indexPointer;
}

void BlastSeqSrcSetPointers(BlastSeqSrc* seq_src, int* inputA, int* inputB, int* result)
{
    seq_src->inputA = inputA;
    seq_src->inputB = inputB;
    seq_src->result = result;
}

void BlastSeqSrcGetPointers(const BlastSeqSrc* seq_src, int** inputA, int** inputB, int** result)
{
    *inputA = seq_src->inputA;
    *inputB = seq_src->inputB;
    *result = seq_src->result;
}

void find_balanced_indices(uint64_t* prefix_sum, uint64_t size, uint64_t num_parts, uint64_t* indices, uint64_t prefix_diff) {
    uint64_t total_sum = prefix_sum[size - 1] - prefix_diff;
    double part_sum = (double)total_sum / num_parts;

    for (uint64_t part = 1; part < num_parts; ++part) {
        double target = part * part_sum;
        uint64_t left = 0, right = size - 1;
        uint64_t best_index = -1;
        double best_diff = (double)INT_MAX;

        while (left <= right) {
            uint64_t mid = (left + right) / 2;
            uint64_t current_sum = prefix_sum[mid] - prefix_diff;
            double diff = fabs((double)current_sum - target);

            if (diff < best_diff) {
                best_diff = diff;
                best_index = mid;
            }

            if (current_sum < target) {
                left = mid + 1;
            } else {
                right = mid - 1;
            }
        }

        indices[part - 1] = best_index;
    }
}

void BlastSeqSrcSetHostDBPointer(BlastSeqSrc* seq_src, char* db_pointer){
    seq_src->hostDBPtr = db_pointer;
}

void BlastSeqSrcGetHostDBPointer(const BlastSeqSrc* seq_src, char** db_pointer){
    *db_pointer = seq_src->hostDBPtr;
}

void BlastSeqSrcSetSequenceBitmaskPointer(BlastSeqSrc* seq_src, char* ptr){
    seq_src->survivingSeedsPtr = ptr;
}

void BlastSeqSrcGetSequenceBitmaskPointer(const BlastSeqSrc* seq_src, char** ptr){
    *ptr = seq_src->survivingSeedsPtr;
}



void BlastSeqSrcSetQueryPointers(BlastSeqSrc* seq_src, char* query_pointer, char* query_bitmask_pointer){
    seq_src->gpuQueryDB = query_pointer;
    seq_src->gpuQueryBitmask = query_bitmask_pointer;
}

void BlastSeqSrcGetQueryPointers(const BlastSeqSrc* seq_src, char** db_pointer, char** db_bitmask_pointer){
    *db_pointer = seq_src->gpuQueryDB;
    *db_bitmask_pointer = seq_src->gpuQueryBitmask;
}

void BlastSeqSrcSetIsolationPointers(BlastSeqSrc* seq_src, uint16_t* lookuptable, uint16_t* prefixTable, uint16_t* seedInformationTable){
    seq_src->qlookupTable = lookuptable;
    seq_src->qlookupPrefix = prefixTable;
    seq_src->actualLookupTable = seedInformationTable;
}

void BlastSeqSrcGetIsolationPointers(const BlastSeqSrc* seq_src, uint16_t** lookuptable, uint16_t** prefixTable, uint16_t** seedInformationTable){
    *lookuptable = seq_src->qlookupTable;
    *prefixTable = seq_src->qlookupPrefix;
    *seedInformationTable = seq_src->actualLookupTable;
}

void BlastSeqSrcPassDBWorkInfo(BlastSeqSrc* seq_src, uint64_t* indexArray, uint64_t size, uint64_t num_parts, uint64_t* startIndex, uint64_t* endIndex, uint64_t * dbOffset) {
    seq_src->indexPointer = indexArray;
    seq_src->sizeOfPointerArray = size;
    seq_src->startChunkArray = startIndex;
    seq_src->endChunkArray = endIndex;
    seq_src->sizeOfChunkArray = num_parts;
    seq_src->dbOffset = dbOffset;
}

void BlastSeqSrcGetDBWorkInfo(const BlastSeqSrc* seq_src, uint64_t** indexArray, uint64_t* size, uint64_t* num_parts, uint64_t* startIndex, uint64_t* endIndex, uint64_t * dbOffset) {
    *indexArray = seq_src->indexPointer;
    *size = seq_src->sizeOfPointerArray;
    *startIndex = seq_src->startChunkArray;
    *endIndex = seq_src->endChunkArray;
    *num_parts = seq_src->sizeOfChunkArray;
    *dbOffset = seq_src->dbOffset;
}

void BlastSeqSrcPassDBWorkInfoGPU(BlastSeqSrc* seq_src, uint64_t* indexArray, uint64_t size, uint64_t num_parts, uint64_t* startIndex, uint64_t* endIndex, uint64_t * dbOffset) {
    seq_src->indexPointerGPU = indexArray;
    seq_src->sizeOfPointerArrayGPU = size;
    seq_src->startChunkArrayGPU = startIndex;
    seq_src->endChunkArrayGPU = endIndex;
    seq_src->sizeOfChunkArrayGPU = num_parts;
    seq_src->dbOffsetGPU = dbOffset;
}

void BlastSeqSrcGetDBWorkInfoGPU(const BlastSeqSrc* seq_src, uint64_t** indexArray, uint64_t* size, uint64_t* num_parts, uint64_t** startIndex, uint64_t** endIndex, uint64_t ** dbOffset) {
    *indexArray = seq_src->indexPointerGPU;
    *size = seq_src->sizeOfPointerArrayGPU;
    *startIndex = seq_src->startChunkArrayGPU;
    *endIndex = seq_src->endChunkArrayGPU;
    *num_parts = seq_src->sizeOfChunkArrayGPU;
    *dbOffset = seq_src->dbOffsetGPU;
}

#ifdef KAPPA_PRINT_DIAGNOSTICS

static const size_t kInitialGiListSize = 10;

Blast_GiList*
Blast_GiListNew(void)
{
    return Blast_GiListNewEx(kInitialGiListSize);
}

Blast_GiList*
Blast_GiListNewEx(size_t list_size)
{
    Blast_GiList* retval = (Blast_GiList*) calloc(1, sizeof(Blast_GiList));
    if ( !retval ) {
        return NULL;
    }

    retval->data = (Int4*) calloc(list_size, sizeof(Int4));
    if ( !retval->data ) {
        return Blast_GiListFree(retval);
    }
    retval->num_allocated = list_size;

    return retval;
}

Blast_GiList*
Blast_GiListFree(Blast_GiList* gilist)
{
    if ( !gilist ) {
        return NULL;
    }
    if (gilist->data) {
        sfree(gilist->data);
    }
    sfree(gilist);
    return NULL;
}

Int2
Blast_GiList_ReallocIfNecessary(Blast_GiList* gilist)
{
    ASSERT(gilist);

    if (gilist->num_used+1 > gilist->num_allocated) {
        /* we need more room for elements */
        gilist->num_allocated *= 2;
        gilist->data = (Int4*) realloc(gilist->data, 
                                       gilist->num_allocated * sizeof(Int4));
        if ( !gilist->data ) {
            return kOutOfMemory;
        }
    }
    return 0;
}

Int2
Blast_GiList_Append(Blast_GiList* gilist, Int4 gi)
{
    Int2 retval = 0;
    ASSERT(gilist);

    if ( (retval = Blast_GiList_ReallocIfNecessary(gilist)) != 0) {
        return retval;
    }
    gilist->data[gilist->num_used++] = gi;
    return retval;
}

Blast_GiList*
BlastSeqSrcGetGis(const BlastSeqSrc* seq_src, void* oid)
{
    ASSERT(seq_src);
    ASSERT(seq_src->GetSeqLen);
    return (*seq_src->GetGis)(seq_src->DataStructure, oid);
}

#endif /* KAPPA_PRINT_DIAGNOSTICS */

/******************** BlastSeqSrcIterator API *******************************/

BlastSeqSrcIterator* BlastSeqSrcIteratorNew()
{
    return BlastSeqSrcIteratorNewEx(0);
}

const unsigned int kBlastSeqSrcDefaultChunkSize = 1024;

BlastSeqSrcIterator* BlastSeqSrcIteratorNewEx(unsigned int chunk_sz)
{
    BlastSeqSrcIterator* itr = NULL;

    if (chunk_sz == 0)
       chunk_sz = kBlastSeqSrcDefaultChunkSize;

    itr = (BlastSeqSrcIterator*) calloc(1, sizeof(BlastSeqSrcIterator));
    if (!itr) {
        return NULL;
    }

    /* Should employ lazy initialization? */
    itr->oid_list = (int*)malloc(chunk_sz * sizeof(int));
    if (!itr->oid_list) {
        sfree(itr);
        return NULL;
    }

    itr->chunk_sz = chunk_sz;
    itr->current_pos = UINT4_MAX;   /* mark iterator as uninitialized */

    return itr;
}

BlastSeqSrcIterator* BlastSeqSrcIteratorFree(BlastSeqSrcIterator* itr)
{
    if (!itr) {
        return NULL;
    }
    if (itr->oid_list) {
        sfree(itr->oid_list);
    }

    sfree(itr);
    return NULL;
}

Int4 BlastSeqSrcIteratorNext(const BlastSeqSrc* seq_src, 
                             BlastSeqSrcIterator* itr)
{
    ASSERT(seq_src);
    ASSERT(itr);
    ASSERT(seq_src->IterNext);

    return (*seq_src->IterNext)(seq_src->DataStructure, itr);
}

void
BlastSeqSrcResetChunkIterator(BlastSeqSrc* seq_src)
{
    ASSERT(seq_src);
    ASSERT(seq_src->ResetChunkIterator);
    (*seq_src->ResetChunkIterator)(seq_src->DataStructure);
}

BlastSeqSrcSetRangesArg *
BlastSeqSrcSetRangesArgNew(Int4 num_ranges)
{
    BlastSeqSrcSetRangesArg * retv = (BlastSeqSrcSetRangesArg *)
                          malloc(sizeof(BlastSeqSrcSetRangesArg));
    retv->capacity = num_ranges;
    retv->num_ranges = 0;
    retv->ranges = (Int4 *) malloc(2*num_ranges*sizeof(Int4));
    return retv;
}

BlastSeqSrcSetRangesArg *
BlastSeqSrcSetRangesArgFree(BlastSeqSrcSetRangesArg *arg)
{
	if (arg) {
		if (arg->ranges) {
			sfree(arg->ranges);
		}
		sfree(arg);
	}
    return NULL;
}

Int2
BlastSeqSrcSetRangesArgAddRange(BlastSeqSrcSetRangesArg *arg,
                                Int4 begin, Int4 end, Int4 len)
{
    ASSERT(arg);
    if ((arg->num_ranges+2) > arg->capacity) {
        Int4 new_size = arg->capacity*2;
        arg->ranges = (Int4*)realloc((void*)arg->ranges, 
                                     (sizeof(Int4)*new_size*2));
        if (!arg->ranges) {
            return 1;
        }
        arg->capacity = new_size;
    }
    begin = MAX(0, begin - BLAST_SEQSRC_OVERHANG);
    end = MIN(end + BLAST_SEQSRC_OVERHANG, len);
    arg->ranges[arg->num_ranges++] = begin;
    arg->ranges[arg->num_ranges++] = end;
    return 0;
}

static int
BeginCompareHSPs(const void* x, const void* y)
{
    Int4 *r1 = (Int4 *)x;
    Int4 *r2 = (Int4 *)y;
    return (*r1)-(*r2);
}

void
BlastSeqSrcSetRangesArgBuild(BlastSeqSrcSetRangesArg *arg)
{
    Int4 i, j;
    ASSERT(arg);
    arg->num_ranges /= 2;
    if (arg->num_ranges <= 1) return;
    qsort(arg->ranges, arg->num_ranges, 2*sizeof(Int4), BeginCompareHSPs);
    i=0;
    for (j=1; j<arg->num_ranges; ++j) {
        Int4 begin = arg->ranges[j*2];
        Int4 end = arg->ranges[j*2+1];
        ASSERT(begin >= arg->ranges[2*i]);
        if (begin > arg->ranges[2*i+1] + BLAST_SEQSRC_MINGAP) {
            /* insert as a new range */
            ++i;
            arg->ranges[i*2] = begin;
            arg->ranges[i*2+1] = end;
        } else if (end > arg->ranges[2*i+1]) {
            /* merge into the previous range */
            arg->ranges[i*2+1] = end;
        }
    }
    arg->num_ranges = i+1;
}

/*****************************************************************************/

/* The following macros implement the "member functions" of the BlastSeqSrc
 * structure so that its various fields can be accessed by the implementors of
 * the BlastSeqSrc interface */

#ifndef SKIP_DOXYGEN_PROCESSING

#define DEFINE_BLAST_SEQ_SRC_MEMBER_FUNCTIONS(member_type, member) \
DEFINE_BLAST_SEQ_SRC_ACCESSOR(member_type, member) \
DEFINE_BLAST_SEQ_SRC_MUTATOR(member_type, member)

#define DEFINE_BLAST_SEQ_SRC_ACCESSOR(member_type, member) \
member_type \
_BlastSeqSrcImpl_Get##member(const BlastSeqSrc* var) \
{ \
    if (var) \
        return var->member; \
    else \
        return (member_type) NULL; \
}

#define DEFINE_BLAST_SEQ_SRC_MUTATOR(member_type, member) \
void \
_BlastSeqSrcImpl_Set##member(BlastSeqSrc* var, member_type arg) \
{ if (var) var->member = arg; }

#endif

/* Note there's no ; after these macros! */
DEFINE_BLAST_SEQ_SRC_MEMBER_FUNCTIONS(BlastSeqSrcConstructor, NewFnPtr)
DEFINE_BLAST_SEQ_SRC_MEMBER_FUNCTIONS(BlastSeqSrcDestructor, DeleteFnPtr)
DEFINE_BLAST_SEQ_SRC_MEMBER_FUNCTIONS(BlastSeqSrcCopier, CopyFnPtr)

DEFINE_BLAST_SEQ_SRC_MEMBER_FUNCTIONS(void*, DataStructure)
DEFINE_BLAST_SEQ_SRC_MEMBER_FUNCTIONS(char*, InitErrorStr)
DEFINE_BLAST_SEQ_SRC_MEMBER_FUNCTIONS(SetInt4FnPtr, SetNumberOfThreads)
DEFINE_BLAST_SEQ_SRC_MEMBER_FUNCTIONS(GetInt4FnPtr, GetNumSeqs)
DEFINE_BLAST_SEQ_SRC_MEMBER_FUNCTIONS(GetInt4FnPtr, GetNumSeqsStats)
DEFINE_BLAST_SEQ_SRC_MEMBER_FUNCTIONS(GetInt4FnPtr, GetMaxSeqLen)
DEFINE_BLAST_SEQ_SRC_MEMBER_FUNCTIONS(GetInt4FnPtr, GetMinSeqLen)
DEFINE_BLAST_SEQ_SRC_MEMBER_FUNCTIONS(GetInt4FnPtr, GetAvgSeqLen)
DEFINE_BLAST_SEQ_SRC_MEMBER_FUNCTIONS(GetInt8FnPtr, GetTotLen)
DEFINE_BLAST_SEQ_SRC_MEMBER_FUNCTIONS(GetInt8FnPtr, GetTotLenStats)

DEFINE_BLAST_SEQ_SRC_MEMBER_FUNCTIONS(GetStrFnPtr, GetName)
DEFINE_BLAST_SEQ_SRC_MEMBER_FUNCTIONS(GetBoolFnPtr, GetIsProt)

DEFINE_BLAST_SEQ_SRC_MEMBER_FUNCTIONS(GetBoolFnPtr, GetSupportsPartialFetching)
DEFINE_BLAST_SEQ_SRC_MEMBER_FUNCTIONS(SetSeqRangeFnPtr, SetSeqRange)

DEFINE_BLAST_SEQ_SRC_MEMBER_FUNCTIONS(GetSeqBlkFnPtr, GetSequence)
DEFINE_BLAST_SEQ_SRC_MEMBER_FUNCTIONS(GetInt4FnPtr, GetSeqLen)
DEFINE_BLAST_SEQ_SRC_MEMBER_FUNCTIONS(ReleaseSeqBlkFnPtr, ReleaseSequence)

DEFINE_BLAST_SEQ_SRC_MEMBER_FUNCTIONS(AdvanceIteratorFnPtr, IterNext)
#ifdef KAPPA_PRINT_DIAGNOSTICS
DEFINE_BLAST_SEQ_SRC_MEMBER_FUNCTIONS(GetGisFnPtr, GetGis)
#endif /* KAPPA_PRINT_DIAGNOSTICS */
DEFINE_BLAST_SEQ_SRC_MEMBER_FUNCTIONS(ResetChunkIteratorFnPtr, 
                                      ResetChunkIterator)
