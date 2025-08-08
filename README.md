# BLAZE: Exploiting Hybrid Parallelism and Size-customized Kernels to Accelerate BLASTP on GPUs

## Requirements

- **Operating System**: Ubuntu 24.04.2 LTS

- **Dependencies**:
  - `libboost-all-dev` (v1.83.0.1)
  - `seqkit` (v2.10.0)
  - CUDA Toolkit (release 12.6)
  - CUDA Driver (version 570.124.06)
  - `ncbi-blast+` (version 2.13.0)

> **Note:** Refer to the [NCBI-BLAST+ toolkit documentation](https://www.ncbi.nlm.nih.gov/books/NBK279690/) for additional build dependencies.

BLAZE is free for academic and non-commercial use. The original NCBI BLAST+ license and copyright notices are retained for all modified components.

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

4. **Generate result graphs:**

   ```bash
   ./generate_graphs.sh
   ```

