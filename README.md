# BLAZE: Exploiting Hybrid Parallelism and Size-customized Kernels to Accelerate BLASTP on GPUs

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

The original NCBI BLAST+ license and copyright notices are retained for all modified components.

## Requirements

- **Operating System**: Ubuntu 24.04.2 LTS

- **Dependencies**:
  - `libboost-all-dev` (v1.83.0.1)
  - `seqkit` (v2.10.0)
  - CUDA Toolkit (release 12.6)
  - CUDA Driver (version 570.124.06)
  - `ncbi-blast+` (version 2.13.0)
  - python packages
      - numpy
      - biopython
      - matplotlib

> **Note:** Refer to the [NCBI-BLAST+ toolkit documentation](https://www.ncbi.nlm.nih.gov/books/NBK279690/) for additional build dependencies.

## Parameters 
There are several preprocessor defines across the src that control various aspects of the code.

Important ones you’ll likely want to adjust are listed below. <br/>
API\_READ\_THREADS\_GPU - Defines the number of CPU threads that copy the portion of the database assigned to the GPU. (Default 8, which is half of the total 16 threads that are passed in `blaze_runner.sh`) <br/>
GPU\_MEM\_PERCENT - Defines the percentage of the GPU’s total global memory reserved per CPU thread that manages data copies to the GPU. (Default 0.05)<br/>
NUM\_BUFFERS - Defines the number of buffers that are asynchronously prefilled with GPU DB data while the GPU executes kernels. (Default 5.0)

> Note: These three are related; ensure they’re configured so the total GPU memory reserved does not exceed available VRAM.

MAX\_SEEDS - Maximum number of ungapped seeds stored in shared memory before the sequence pair is offloaded to the CPU for reevaluation. (Default 8)<br/>
MAX\_SEEDS\_GAPPED - Maximum number of gapped seeds stored in shared memory before the sequence pair is offloaded to the CPU for reevaluation. (Default 8)<br/>
bbins - Defines the DB bin boundaries. (Default 12)<br/>
TQ\_VALUES - Defines the query bitset sizes. <br/>
SUBJECT\_LEN\_BINS - Defines the subject-length bin boundaries; if modified, also update X(Q) within `util.hpp` and the unordered\_maps kernelMap and kernelMapEverything.



## Navigating the code

```bash
.
├── BLAZE_gpu # Contains the gpu src code
│   ├── gpuBLAZE.cu
│   ├── gpu_cpu_common.h
│   └── util.hpp
├── BLAZE_mk # Contains modified makefiles
│   ├── src_app_blastdb_Makefile.in
│   ├── src_app_blastdb_Makefile.makeblastdb.app
│   ├── src_app_blast_Makefile.blastp.app
│   ├── src_app_blast_Makefile.in
│   └── src_app_Makefile.in
├── BLAZE_ncbi # Contains modified ncbi src files
│   ├── blast_engine.c
│   ├── blast_engine.h
│   ├── blast_seqsrc.c
│   ├── blast_seqsrc.h
│   ├── prelim_search_runner.hpp
│   ├── prelim_stage.cpp
│   ├── seqdbimpl.cpp
│   └── setup_factory.hpp
├── BLAZE_queries # Contains the benchmark queries
│   ├── NP_085504.fa
│   ├── NP_174344.fa
│   ├── NP_214548.fa
│   ├── NP_245243.fa
│   ├── NP_308413.fa
│   ├── NP_345015.fa
│   ├── NP_349135.fa
│   ├── NP_436070.fa
│   ├── NP_566269.fa
│   ├── NP_611552.fa
│   ├── NP_633288.fa
│   ├── NP_670887.fa
│   ├── NP_706440.fa
│   ├── NP_733286.fa
│   ├── NP_781168.fa
│   ├── NP_784621.fa
│   ├── NP_818972.fa
│   ├── NP_864591.fa
│   ├── NP_872878.fa
│   ├── NP_894055.fa
│   ├── NP_931067.fa
│   ├── NP_937947.fa
│   ├── NP_952492.fa
│   ├── NP_965272.fa
│   ├── XP_237451.fa
│   ├── XP_499151.fa
│   ├── XP_610726.fa
│   ├── YP_000338.fa
│   ├── YP_054627.fa
│   ├── YP_062020.fa
│   ├── YP_069013.fa
│   ├── YP_081219.fa
│   ├── YP_121326.fa
│   ├── YP_129633.fa
│   ├── YP_137100.fa
│   ├── YP_141079.fa
│   ├── YP_182226.fa
│   ├── YP_189837.fa
│   ├── YP_209371.fa
│   └── YP_244862.fa
├── blaze_runner.sh
├── build.sh # Builds BLAZE
├── generate_graphs.sh # Generate results graph
├── graphs.py
├── LICENSE
├── perf_tests.sh # Runs st, mt and BLAZE
└── README.md
```

## Installation and Deployment

1. **Clone repository:**

   ```bash
   git clone https://github.com/181charan/BLAZE-SC-ad
   ```

2. **Build BLAZE executable:**

   ```bash
   ./build.sh
   ```

The built executable supports executing default searches in the same way as the standard NCBI BLASTP binary.

3. **Run performance tests:**

   ```bash
   ./perf_test.sh
   ```
> **Important Notes: 1:** Before use, sort the database FASTA by sequence length and run the default `makeblastdb` from ncbi-blast-2.13.0+; BLAZE expects this format.

> **2:** The nr database can be downloaded from the NCBI ftp server located at https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/

> **3:** Ensure to edit `NUM_THREADS` in blaze_runner.sh and `API_READ_THREADS_GPU` depending on the CPU.

> **4:** The initial BLAZE execution functions as a warm-up run, generating one-time indices. The corresponding runtime, approximately 6 minutes, is discarded from the reported measurements.

4. **Generate result graphs:**

   ```bash
   ./generate_graphs.sh
   ```
   The graph is saved as `main_blaze.png`

