/*  $Id: prelim_search_runner.hpp 514950 2016-09-27 14:52:47Z rackerst $
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
 */

/** @file prelim_search_runner.hpp
 * Defines internal auxiliary functor object to run the preliminary stage of
 * the BLAST search.
 */

#ifndef ALGO_BLAST_API___PRELIM_SEARCH_RUNNER__HPP
#define ALGO_BLAST_API___PRELIM_SEARCH_RUNNER__HPP

/** @addtogroup AlgoBlast
 *
 * @{
 */

#include <corelib/ncbithr.hpp>                  // for CThread
#include <algo/blast/api/setup_factory.hpp>
#include "blast_memento_priv.hpp"

// CORE BLAST includes
#include <algo/blast/core/blast_engine.h>

// Include the GPU header file
#include <algo/blast/core/gpu_cpu_common.h>

#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <functional>
#include <vector>
#include <atomic>

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(blast)

// class ThreadPool {
// public:
//     ThreadPool(size_t numThreads) : stop(false) {
//         for(size_t i = 0; i < numThreads; i++) {
//             workers.emplace_back([this] {
//                 while(true) {
//                     std::function<void()> task;
//                     {
//                         std::unique_lock<std::mutex> lock(queueMutex);
//                         cv.wait(lock, [this]{ 
//                             return stop || !tasks.empty(); 
//                         });
//                         if(stop && tasks.empty()) {
//                             return;
//                         }
//                         task = std::move(tasks.front());
//                         tasks.pop();
//                     }
//                     task(); // do work
//                 }
//             });
//         }
//     }

//     void enqueue(std::function<void()> f) {
//         {
//             std::unique_lock<std::mutex> lock(queueMutex);
//             tasks.push(std::move(f));
//         }
//         cv.notify_one();
//     }

//     void waitAll() {
//         // Wait until all tasks have finished
//         // This is a no-op unless you have a separate mechanism
//         // to track how many tasks are in flight, e.g. an atomic counter.
//     }

//     ~ThreadPool() {
//         {
//             std::unique_lock<std::mutex> lock(queueMutex);
//             stop = true;
//         }
//         cv.notify_all();
//         for(auto &worker : workers) {
//             worker.join();
//         }
//     }

// private:
//     std::vector<std::thread> workers;
//     std::queue<std::function<void()>> tasks;
//     std::mutex queueMutex;
//     std::condition_variable cv;
//     bool stop;
// };




/// Functor to run the preliminary stage of the BLAST search
class CPrelimSearchRunner : public CObject
{
public:
    CPrelimSearchRunner(SInternalData& internal_data,
                        const CBlastOptionsMemento* opts_memento)
        : m_InternalData(internal_data), m_OptsMemento(opts_memento)
    {}
    ~CPrelimSearchRunner() {}
    int operator()() {
        _ASSERT(m_OptsMemento);
        _ASSERT(m_InternalData.m_Queries);
        _ASSERT(m_InternalData.m_QueryInfo);
        _ASSERT(m_InternalData.m_SeqSrc);
        _ASSERT(m_InternalData.m_ScoreBlk);
        _ASSERT(m_InternalData.m_LookupTable);
        _ASSERT(m_InternalData.m_HspStream);

        SBlastProgressReset(m_InternalData.m_ProgressMonitor->Get());

        Int2 retval = 0;

        if(m_InternalData.pflag == 1) {

            // Set the is GPU flag to 1
            BlastSeqSrcSetIsGPU(m_InternalData.m_SeqSrc->GetPointer(), 1);

            // size_t projectedSurvivor = 0.1 * BlastSeqSrcGetNumSeqs(m_InternalData.m_SeqSrc->GetPointer());

            // uint32_t * survivingSequences = new uint32_t[projectedSurvivor];

            // uint32_t num_survivors = 0;

            // BLAZE_GPU(m_InternalData.m_SeqSrc->GetPointer(), m_InternalData.threadIndex, m_InternalData.threadCount, m_InternalData.m_Queries->length, survivingSequences, &num_survivors);

            // // Set the survivors pointer to the BlastSeqSrc
            // BlastSeqSrcSetSurvivors(m_InternalData.m_SeqSrc->GetPointer(), survivingSequences, num_survivors);

            // // Create a file and append the thread index to the file name and write the survivors to the file as a CSV file
            // std::string fileName = "survivors_" + std::to_string(m_InternalData.threadIndex) + ".csv";

            // std::ofstream outFile(fileName, std::ios::binary);

            // // Write the number of survivors to the file as a csv file
            // for (uint32_t i = 0; i < num_survivors; i++) {
            //     outFile << survivingSequences[i] << std::endl;
            // }

            // // Close the file
            // outFile.close();

            // Print the size of the vector to stderr
            // std::cerr << " Potential Sequences: " << num_survivors << std::endl;

            retval =            Blast_RunPreliminarySearchWithInterruptMT(
                                 m_InternalData.threadCount,
                                 m_InternalData.threadIndex,
                                 m_OptsMemento->m_ProgramType,
                                 m_InternalData.m_Queries,
                                 m_InternalData.m_QueryInfo,
                                 m_InternalData.m_SeqSrc->GetPointer(),
                                 m_OptsMemento->m_ScoringOpts,
                                 m_InternalData.m_ScoreBlk->GetPointer(),
                                 m_InternalData.m_LookupTable->GetPointer(),
                                 m_OptsMemento->m_InitWordOpts,
                                 m_OptsMemento->m_ExtnOpts,
                                 m_OptsMemento->m_HitSaveOpts,
                                 m_OptsMemento->m_EffLenOpts,
                                 m_OptsMemento->m_PSIBlastOpts,
                                 m_OptsMemento->m_DbOpts,
                                 m_InternalData.m_HspStream->GetPointer(),
                                 m_InternalData.m_Diagnostics->GetPointer(),
                                 m_InternalData.m_FnInterrupt,
                                 m_InternalData.m_ProgressMonitor->Get());

                
        }else {
            //fHello_Two(m_InternalData.threadIndex);

            // Set the is GPU flag to 0
            BlastSeqSrcSetIsGPU(m_InternalData.m_SeqSrc->GetPointer(), 0);

            retval = Blast_RunPreliminarySearchWithInterruptMT(
                                 m_InternalData.threadCount,
                                 m_InternalData.threadIndex,
                                 m_OptsMemento->m_ProgramType,
                                 m_InternalData.m_Queries,
                                 m_InternalData.m_QueryInfo,
                                 m_InternalData.m_SeqSrc->GetPointer(),
                                 m_OptsMemento->m_ScoringOpts,
                                 m_InternalData.m_ScoreBlk->GetPointer(),
                                 m_InternalData.m_LookupTable->GetPointer(),
                                 m_OptsMemento->m_InitWordOpts,
                                 m_OptsMemento->m_ExtnOpts,
                                 m_OptsMemento->m_HitSaveOpts,
                                 m_OptsMemento->m_EffLenOpts,
                                 m_OptsMemento->m_PSIBlastOpts,
                                 m_OptsMemento->m_DbOpts,
                                 m_InternalData.m_HspStream->GetPointer(),
                                 m_InternalData.m_Diagnostics->GetPointer(),
                                 m_InternalData.m_FnInterrupt,
                                 m_InternalData.m_ProgressMonitor->Get());
        }


        // Int2 retval = Blast_RunPreliminarySearchWithInterruptMT(
        //                          m_InternalData.threadCount,
        //                          m_InternalData.threadIndex,
        //                          m_OptsMemento->m_ProgramType,
        //                          m_InternalData.m_Queries,
        //                          m_InternalData.m_QueryInfo,
        //                          m_InternalData.m_SeqSrc->GetPointer(),
        //                          m_OptsMemento->m_ScoringOpts,
        //                          m_InternalData.m_ScoreBlk->GetPointer(),
        //                          m_InternalData.m_LookupTable->GetPointer(),
        //                          m_OptsMemento->m_InitWordOpts,
        //                          m_OptsMemento->m_ExtnOpts,
        //                          m_OptsMemento->m_HitSaveOpts,
        //                          m_OptsMemento->m_EffLenOpts,
        //                          m_OptsMemento->m_PSIBlastOpts,
        //                          m_OptsMemento->m_DbOpts,
        //                          m_InternalData.m_HspStream->GetPointer(),
        //                          m_InternalData.m_Diagnostics->GetPointer(),
        //                          m_InternalData.m_FnInterrupt,
        //                          m_InternalData.m_ProgressMonitor->Get());

        

        return static_cast<int>(retval);
    }

private:
    /// Data structure containing all the needed C structures for the
    /// preliminary stage of the BLAST search
    SInternalData& m_InternalData;

    /// Pointer to memento which this class doesn't own
    const CBlastOptionsMemento* m_OptsMemento;


    /// Prohibit copy constructor
    CPrelimSearchRunner(const CPrelimSearchRunner& rhs);
    /// Prohibit assignment operator
    CPrelimSearchRunner& operator=(const CPrelimSearchRunner& rhs);
};

/// Thread class to run the preliminary stage of the BLAST search
class CPrelimSearchThread : public CThread
{
public:
    CPrelimSearchThread(SInternalData& internal_data,
                        const CBlastOptionsMemento* opts_memento)
        : m_InternalData(internal_data), m_OptsMemento(opts_memento)
    {

        // // Print that you were here!
        // cerr << "Printing from the constructor!!!" << endl;

        // The following fields need to be copied to ensure MT-safety
        BlastSeqSrc* seqsrc =
            BlastSeqSrcCopy(m_InternalData.m_SeqSrc->GetPointer());
        m_InternalData.m_SeqSrc.Reset(new TBlastSeqSrc(seqsrc,
                                                       BlastSeqSrcFree));
        // The progress field must be copied to ensure MT-safety
        if (m_InternalData.m_ProgressMonitor->Get()) {
            SBlastProgress* bp =
                SBlastProgressNew(m_InternalData.m_ProgressMonitor->Get()->user_data);
            m_InternalData.m_ProgressMonitor.Reset(new CSBlastProgress(bp));
        }
        // The BlastQueryInfo field needs to be copied to silence Thread
        // Sanitizer warnings, and probably to ensure MT-safety too.
        BlastQueryInfo* queryInfo =
                BlastQueryInfoDup(m_InternalData.m_QueryInfo);
        m_InternalData.m_QueryInfo = queryInfo;
    }

protected:
    virtual ~CPrelimSearchThread(void) {
        BlastQueryInfoFree(m_InternalData.m_QueryInfo);
    }

    virtual void* Main(void) {

        return (void*)
            ((intptr_t) CPrelimSearchRunner(m_InternalData, m_OptsMemento)());
    }

private:
    SInternalData m_InternalData;
    const CBlastOptionsMemento* m_OptsMemento;
};

END_SCOPE(blast)
END_NCBI_SCOPE

/* @} */

#endif /* ALGO_BLAST_API___PRELIM_SEARCH_RUNNER__HPP */
