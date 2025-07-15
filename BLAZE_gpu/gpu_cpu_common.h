
#ifndef GPU_CPU_COMMON_H
#define GPU_CPU_COMMON_H

#include <algo/blast/core/blast_extend.h>
#include <algo/blast/core/blast_seqsrc_impl.h>
#include <algo/blast/core/blast_seqsrc.h>

#ifdef __cplusplus
extern "C"
{
#endif

    #define GPU_MEM_PERCENT 0.05

    #define BITMASK 256

    #define BYTE_SIZE 8

    #define TABLESIZE 21952

    #define CPU_OFFLOAD_SIZE_GT_EQ 1024
    #define CPU_OFFLOAD_BIN 10 // Sequences with over 1024bp will fall in bin 10 and 11

    #define API_READ_THREADS_GPU 8

    #define NUM_BUFFERS 5

    int BLAZE_GPU(BlastSeqSrc* seq_src, uint32_t threadIndex, uint32_t totalThreads, uint32_t queryLen, uint32_t * survivingSeqs, uint32_t * num_survivors, int ungapped_cutoff, int gapped_cutoff);

    int mainCall(uint8_t * dbArray_cpu, uint8_t * gpuDbArray, size_t dbArraySize, uint8_t * survivals_cpu, uint8_t *gpuSurvivals , uint32_t * sequenceLengths, uint32_t * gpuSequenceLengths, uint8_t * gpuQueryDbArray, char * gpuQuery, uint16_t querySize, size_t OIDs_processed, void* stream, uint32_t bitmask_read_en);

    int mainGappedCall(uint8_t * dbArray_cpu, uint8_t * gpuDbArray, size_t dbArraySize, uint8_t * survivals_cpu, uint8_t *gpuSurvivals , uint32_t * sequenceLengths, uint32_t * gpuSequenceLengths, uint8_t * gpuQueryDbArray, char * gpuQuery, uint16_t querySize, size_t OIDs_processed, void* stream, uint32_t bitmask_read_en);

    void createCudaStream(void** streamPtr);

    void copySurvivals(uint8_t * survivals_cpu, uint8_t *gpuSurvivals, size_t OIDs_processed, void* stream);

    void deleteCudaStream(void* stream);

    void deviceMemCpyToDeviceAsync(void *dst, void *src, size_t size, void *stream);

    void copyConstantMemory(uint8_t * gpuQuerybitmaskArray, size_t querybitmaskarray_gpu_size);

    void copyConstantMemoryLUT(uint16_t * seedLUT, size_t seedLUTSize);


    void streamSynchronize(void *stream);

    void hostAlloc(void** pointer, size_t size) ;

    void deviceAlloc(void** pointer, size_t size);

    void deviceMemCpyToDevice(void*dst, void* src, size_t size);

    void hostFree(void* pointer);

    void deviceFree(void* pointer);

    size_t getTotalGPUMemory();

    size_t getFreeGPUMemory();
#ifdef __cplusplus
}
#endif

#endif /* GPU_CPU_COMMON_H */
