
# Validated hardware, operating systems and compilers

hpcscan has been successfully tested on the hardware, operating systems and compilers listed below.

Operating system | Compiler | MPI | Host | Device | Test modes 
|----------------|----------|-----|------|------------|----------
Ubuntu 24.04.4 LTS | g++ (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 | mpirun (Open MPI) 5.0.9 | AMD Ryzen 9 7940HX with Radeon Graphics | - |  Baseline, CacheBlk
Rocky Linux 9.7 (Blue Onyx) | g++ (GCC) 15.2.0 | mpirun (Open MPI) 5.0.9 | AMD EPYC 9555 64-Core Processor | - |  Baseline, CacheBlk
Ubuntu 22.04.1 LTS | g++ (Ubuntu 11.3.0-1ubuntu1~22.04) 11.3.0 | mpirun (Open MPI) 4.1.2 | Intel(R) Core(TM) i5-7200U CPU @ 2.50GHz **(Intel Kaby Lake)** | - |  Baseline, CacheBlk
Ubuntu 22.04.1 LTS | Intel icpc (ICC) 2021.7.0 20220726 | Intel MPI Version 2021.7 | Intel(R) Core(TM) i5-7200U CPU @ 2.50GHz **(Intel Kaby Lake)** | - |  Baseline, CacheBlk
Red Hat 4.8.5-39 | Intel oneAPI DPC++/C++ Compiler 2022.1.0 | Intel MPI Version 2021.6 | Intel(R) Xeon(R) Gold 6240L CPU @ 2.60GHz **(Intel Cascade Lake)** | - |  Baseline, CacheBlk
Red Hat 4.8.5-39 | <li> Intel oneAPI DPC++/C++ Compiler 2022.1.0 </li> <li> NVIDIA nvcc release 11.7 </li> | Intel MPI Version 2021.6 | Intel(R) Xeon(R) Gold 6240L CPU @ 2.60GHz **(Intel Cascade Lake)** | Tesla V100S-PCI **(NVIDIA GPU)** |  Baseline, CacheBlk, Cuda, Cuda_Opt, Cuda_Ref
Red Hat 8.5.0-10 | NEC nc++ (NCC) 4.0.0 | NEC MPI 3.1.0 | Intel(R) Xeon(R) Gold 6126 CPU @ 2.60GHz **(Intel Skylake)** | NEC SX-Aurora TSUBASA 20B-P **(NEC Vector Engine)** | Baseline, CacheBlk, NEC, NEC_SCA
Red Hat 8.5.0-10 | Intel oneAPI DPC++/C++ Compiler 2022.1.0 | Intel MPI Version 2021.6 | Intel(R) Xeon(R) Gold 6126 CPU @ 2.60GHz **(Intel Skylake)** | - | Baseline, CacheBlk
SUSE Linux Enterprise Server 15 | Intel icpc (ICC) 19.0.5.281 20190815 | - | Intel(R) Xeon(R) CPU E5-2698 v3 @ 2.30GHz **(Intel Haswell)** | - |  - |
Red Hat 4.8.5-39 | Intel icpc version 19.1.2.254 | - | Intel(R) Xeon(R) Gold 6248 CPU @ 2.50GHz **(Intel Cascade Lake)** | - |  - 
Ubuntu 20.04.1 LTS |  <li> gcc version 9.3.0 </li> <li> NVIDIA nvcc release 11.3, V11.3.109 </li> | - | Intel(R) Core(TM) i7-1065G7 CPU @ 1.30GHz **(Intel Ice Lake)** | GP108M [GeForce MX330] **(NVIDIA GPU)** | - 
CentOS Linux release 7.7.1908 |  <li> Intel icpc (ICC) 19.1.0.166 20191121 </li> <li> NVIDIA nvcc release 11.0, V11.0.167 </li> | - | Intel(R) Xeon(R) Gold 6142 CPU @ 2.60GHz **(Intel Skylake)** | GV100GL [Tesla V100 SXM2 32GB] **(NVIDIA GPU)** |  - 
Ubuntu 20.04.1 LTS | Intel(R) oneAPI DPC++ Compiler 2021.2.0 | - | Intel(R) Core(TM) i7-1065G7 CPU @ 1.30GHz **(Intel Ice Lake)** | - | - 
Ubuntu 20.04.1 LTS | <li> g++ 9.3.0 </li> <li> AMD hipcc 4.2.21155-37cb3a34 </li> | - | AMD EPYC 7742 64-Core Processor @ 2.25GHz **(AMD Rome)** | [AMD Instinct MI100] **(AMD GPU)** | - 
