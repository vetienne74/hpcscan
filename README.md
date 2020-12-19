
![hpcscan logo](misc/hpcscanLogo/hpcscanLogo.jpg)

**TABLE OF CONTENTS**

- [Welcome to hpcscan](#welcome-to-hpcscan)
- [Overview](#overview)
  * [Description](#description)
  * [Project directories](#project-directories)
  * [List of test cases](#list-of-test-cases)
  * [List of test modes](#list-of-test-modes)
- [Environment set-up](#environment-set-up)
  * [Basic requirements](#basic-requirements)
  * [Optional requirements](#optional-requirements)
  * [Environment script (mandatory)](#environment-script-mandatory)
- [Compilation](#compilation)
  * [Makefile](#makefile)
  * [Enabled test modes](#enabled-test-modes)
- [Validation](#validation)
  * [Validation tests](#validation-tests)
  * [Validated hardware, operating systems and compilers](#validated-hardware-operating-systems-and-compilers)
- [Execution](#execution)
  * [Usage](#usage)
  * [Input and output](#input-and-output)
- [Performance benchmarks](#performance-benchmarks)
- [Versions](#versions)
- [Have fun!](#have-fun)

# Welcome to hpcscan

Version 1.1

Contact: Vincent Etienne (Saudi Aramco / EXPEC ARC / GPT) / Email: vetienne@rocketmail.com

Contributors
* Vincent Etienne
* Suha Kayum
* Marcin Rogowski
* Laurent Gatineau

# Overview

## Description

hpcscan is a tool for benchmarking algorithms/kernels that are found in many scientific applications on various architectures/systems.

It features several categories of test cases aiming to measure memory, computation and communication bandwidths.

All test cases are validated with embedded reference solutions.

## Project directories

* `bin` this directory is created during compilation and contains hpcscan executable
* `build` hpcscan can be compiled from here
* `doc` documents
* `env` scripts to initialize hpcscan environment
* `mics` output samples and studies 
* `script` scripts for validation and performance benchmarks
* `src` all hpcscan source files

## List of test cases

Test case name | Description | Remark
------------ | ----------- | ------------
Comm         | **MPI communications bandwidth** <ul><li>Uni-directional (Half-duplex with MPI_Send) proc1 -> proc2</li><li>Bi-directional (Full-duplex with MPI_Sendrecv) proc1 <-> proc2</li><li>Grid halos exchange (MPI_Sendrecv) all procs <-> all procs</li></ul> | <p>This case requires at least 2 MPI processes<br>Width of halos depends on the selected FD stencil order <br> <font color="blue"> **Validation is done with filling grids with predefined values** </font> <br> **Measures GPoints/s and GBytes/s** </p>
FD_D2        | **Finite-difference (second derivatives in space) computations bandwidth** <ul><li> <img src="https://render.githubusercontent.com/render/math?math=U={\partial^2}/{\partial x_1^2} \: (V)"> (for grid dim. 1, 2 or 3) </li> <li> <img src="https://render.githubusercontent.com/render/math?math=U={\partial^2}/{\partial x_2^2} \: (V)"> (for grid dim. 2 or 3) </li>  <li> <img src="https://render.githubusercontent.com/render/math?math=U={\partial^2}/{\partial x_3^2} \: (V)"> (for grid dim. 3) </li> <li> <img src="https://render.githubusercontent.com/render/math?math=U= \Delta (V)"> (for grid dim 2 or 3) </li> </ul> | <p>Accuracy is checked against multi-dimensional sine function<br>Accuracy depends on the selected FD stencil order <br> <font color="blue"> **Computes L1 Error against analitic solution** </font> <br> **Measures GPoints/s, GBytes/s and GFlop/s** </p> 
Grid         | **Grid operations bandwidth** <ul> <li> Fill grid U with constant value </li> <li> Max. diff. between grids U and V </li> <li> L1 norm between U and V </li> <li> Sum of abs(U) </li> <li> Max. of U </li> <li> Min. of U </li> <li> Complex grid manipulation (pressure update in propagator) U = 2 x V - U + C x L </li> <li> Boundary condition (free surface) at all edges of U </li> </ul> | <p>Operation can be done in portions of the grid (for instance, excluding halos) <br> <font color="blue"> **Validation is done with filling grids with predefined values** </font> <br> **Measures GPoints/s and GBytes/s**<p>
Memory       | **Memory operations bandwidth** <ul> <li> Fill array A with constant value </li> <li> Copy array A = B </li> <li> Add 2 arrays A = B + C </li> <li> Multiply 2 arrays A = B * C </li> <li> Add 2 arrays and update array A = A + B </li> </ul>| <p>Conversely to Test Case Grid, operations are done on continuous memory arrays <br> <font color="blue"> **Validation is done with filling arrays with predefined values** </font> <br> **Measures GPoints/s and GBytes/s**<p>
Propa        | **Acoustic wave propagator bandwidth** <ul> <li> 2nd order wave equation </li> <li> <img src="https://render.githubusercontent.com/render/math?math={\partial^2}/{\partial t^2} (P)=c^2 \: \Delta (P)"> </li> <li> Free surface boundary condition is applied to all edges of the domain </li> <li> Wavefield is initialized at t=-dt and t=-2dt with a particular solution </li> </ul> | <p>Accuracy is checked against Eigen mode analytic solution of the wave equation<br>Number of modes can be parametrized differently in every dimension<br>Time step can be set arbitrarily or set to the stability condition<br>Dimension, grid size, and number of time steps can be set arbitrarily<br>Accuracy depends on the selected FD stencil order <br> <font color="blue"> **Computes L1 Error against analitic solution** </font> <br> **Measures GPoints/s, GBytes/s and GFlop/s** </p> 
Template     | Test case template | Copy this template to create a new test case
Util         | Utility tests to check internal functions | Reserved for developpers

## List of test modes

Test mode name | Description | Remark
------------ | ----------- | ------------
Baseline     | CPU standard implementation | -
CacheBlk     | CPU with cache blocking optimization techniques | -
Cuda         | GPU with CUDA without optimization | -
NEC          | NEC with compiler directives | Only available on NEC SX-Aurora TSUBASA
NEC_SCA      | NEC with Stencil Code Accelerator | Only available on NEC SX-Aurora TSUBASA
OpenAcc      | GPU with OpenAcc without optimization | -

# Environment set-up

## Basic requirements

* C++ compiler with OpenMP support
* MPI library

## Optional requirements

* C++ compiler with OpenAcc support 
* CUDA compiler
* NEC compiler

## Environment script (mandatory)

In order to compile and run hpcscan, you need to source one of the files in ./env

Example

`source ./env/setEnvMarsGccCuda.sh`

> **For a new system, you would need to create a file for your system (take example from one of the existing files)**


# Compilation

## Makefile

Go to ./build, and use the command

`make`

[Display command output](misc/fileForReadme/make.txt)

> If hpcscan environment has not been set (see [Environment script (mandatory)](#environment-script-mandatory)), compilation will abort.

By default, hpcscan is compiled in single in precision

To compile in double precision

`make precision=double`

Executable is ./bin/hpcscan

## Enabled test modes

To check the test modes that are enabled in your hpcscan binary, use the command

`./bin/hpcscan -v`

[Display command output](misc/fileForReadme/version.txt)

# Validation

## Validation tests

To check hpcscan has been correctly built and works fine, go to ./script and launch

`sh runValidationTests.sh`

[Display command output](misc/fileForReadme/runValidationTests.txt)

This script runs a set a light test cases and should complete within few minutes (even on a laptop).

You should get in the ouptput report (displayed on the terminal)

* All tests marked as PASSED (370 tests passed per test mode enabled)
* No test marked as FAILED

Check the summary at the end of report to have a quick look on this.

## Validated hardware, operating systems and compilers

hpcscan has been successfully tested on the hardware, operating systems and compilers listed below

Operating system | Compiler | Host (H) | Device (D) | Baseline | CacheBlk | Cuda | NEC | NEC_SCA | OpenAcc
|----------------|----------|----------|------------|----------|----------|------|-----|---------|--------
Ubuntu 20.04.1 LTS |  gcc version 9.3.0 / nvcc release 10.1, V10.1.243 | Intel(R) Core(TM) i7-1065G7 CPU @ 1.30GHz **(Intel Ice Lake)** | GP108M [GeForce MX330] **(NVIDIA GPU)** | OK (H) | OK (H) | ON GOING (D) | - | - | ON GOING (D)
SUSE Linux Enterprise Server 15 | icpc (ICC) 19.0.5.281 20190815 | Intel(R) Xeon(R) CPU E5-2698 v3 @ 2.30GHz **(Intel Haswell)** | - | OK (H) | OK (H) | - | - | - | -
CentOS Linux release 8.1.1911 | nc++ (NCC) 3.1.0 | Intel(R) Xeon(R) Gold 6126 CPU @ 2.60GHz **(Intel Skylake)** | NEC SX-Aurora TSUBASA **(NEC Vector Engine)** | OK (D) | OK (D) | - | OK (D) | OK (D) | -
CentOS Linux release 7.7.1908 |  icpc (ICC) 19.1.0.166 20191121 / nvcc release 11.0, V11.0.167 | Intel(R) Xeon(R) Gold 6142 CPU @ 2.60GHz **(Intel Skylake)** | V100 **(NVIDIA GPU)** | OK (H) | OK (H) | ON GOING (D) | - | - | -
CentOS Linux release 7.7.1908 |  pgc++ 20.1-0 LLVM 64-bit / nvcc release 11.0, V11.0.167 | Intel(R) Xeon(R) Gold 6142 CPU @ 2.60GHz **(Intel Skylake)** | V100 **(NVIDIA GPU)** | OK (H) | OK (H) | ON GOING (D) | - | - | ON GOING (D)

# Execution

## Usage

hpcscan can be launched from a terminal with all configuration parameters within a single line.

**To get help on the parameters**

`./bin/hpcscan -h`

[Display command output](misc/fileForReadme/commandLineParam.txt)

**Execution with a unique MPI process**

`mpirun -n 1 ./bin/hpcscan -testCase <TESTCASE> -testMode <TESTMODE>`

where
* TESTCASE is the name of the test case (see [List of test cases](#list-of-test-cases))
* TESTMODE is the name of the test mode (see [List of test modes](#list-of-test-modes))

Example

`mpirun -n 1 ./bin/hpcscan -testCase Propa -testMode CacheBlk`

[Display command output](misc/fileForReadme/runPropaTestCase.txt)

> If you omit to specify `-testMode <TESTMODE>`, the Baseline mode it assumed.

Example

`mpirun -n 1 ./bin/hpcscan -testCase Propa`

**Execution with a multiples MPI processes**

`mpirun -n <N> ./bin/hpcscan -testCase <TESTCASE> -testMode <TESTMODE> -nsub1 <NSUB1> -nsub2 <NSUB2> -nsub3 <NSUB3>`

> When several MPI processes are used, subdomain decomposition is activated. The product NSUB1 x NSUB2 x NSUB3 must be equal to N (no. of MPI processes).
> You may omit to specify the number of subdomains along an axis if that number is 1.

Example

`mpirun -n 2 ./bin/hpcscan -testCase Comm -nsub1 2`

**Configuration of the grid size and dimension**

Simply add on the command line

`-n1 <N1> -n2 <N2> -n3 <N3> -dim <DIM>`

Where `N1, N2, N3` are the number of grid points along axis 1, 2 and 3.

And `DIM` = 1,2 or 3 (1D, 2D or 3D grids). By default 3D grid is assumed.

Example

`mpirun -n 1 ../bin/hpcscan -testCase Grid -dim 2 -n1 200 -n2 300`

## Input and output

**Input**

hpcscan does not require any input file. All data are built internally.

**Output**

Information regarding results validation and performances are sent to the terminal output.

Performance log file (Description to be added)

Grid output (Description to be added)

Debug trace output (Description to be added)

# Performance benchmarks

These tests are intended to measure various bandwidths

> **They are intensive tests that require to run on HPC platforms**

Performance measurements and scripts to reproduce results can be found in ./doc/TestCases/TestCases.pdf

# Versions

Version      | Description | Release date
------------ | ----------- | ------------
v1.0          | Initial version with test modes Baseline, CacheBlk and NEC_SCA  | Nov 28, 2020
v1.1          | Added test modes NEC, OpenAcc (ON GOING) and Cuda (ON GOING) | Coming soon

# Have fun!

Please share and send your feedbacks to vetienne@rocketmail.com


