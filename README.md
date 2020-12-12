
![hpcscan logo](https://github.com/vetienne74/hpcscan/blob/KAUST-GPU2020/misc/hpcscanLogo/hpcscanLogo.jpg)

**TABLE OF CONTENTS**

- [Welcome to hpcscan](#welcome-to-hpcscan)
- [Overview](#overview)
  * [Description](#description)
  * [List of test cases](#list-of-test-cases)
  * [List of test modes](#list-of-test-modes)
- [Environment set-up](#environment-set-up)
  * [Basic requirements](#basic-requirements)
  * [Optional requirements](#optional-requirements)
  * [Environment script (mandatory)](#environment-script-mandatory)
- [Compilation](#compilation)
  * [Makefile](#makefile)
  * [Enabled test modes](#enabled-test-modes)
  * [Operating systems and compilers](#operating-systems-and-compilers)
- [Validation](#validation)
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

hpcscan is a tool for benchmarking scientific computing kernels on various platforms.

It features several categories of test cases aiming to measure memory, computation and interconnect bandwidth.

All cases are validated with embedded reference solutions.

## List of test cases

Test case name | Description | Remark
------------ | ----------- | ------------
Comm         | MPI communications bandwidth | This case requires at least 2 MPI processes
FD_D2        | Finite-difference computations bandwidth | - 
Grid         | Grid operations bandwidth | -
Memory       | Memory operations bandwidth | -
Propa        | Acoustic wave propagator bandwidth | - 
Template     | Test case template | Copy this template to create a new test case
Util         | Utility tests to check internal functions | Reserved for developpers

## List of test modes

Test mode name | Description | Remark
------------ | ----------- | ------------
Baseline     | CPU standard implementation | -
CacheBlk     | CPU with cache blocking optimization techniques | -
Cuda         | GPU with CUDA without optimization | -
NEC_SCA      | NEC Aurora with Stencil Code Accelerator | -
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

> If hpcscan environment has not been set (see [Environment script (mandatory)](#environment-script-mandatory)), compilation will abort.

By default, hpcscan is compiled in single in precision

To compile in double precision

`make precision=double`

Executable is ./bin/hpcscan

## Enabled test modes

To check the test modes that are enabled, use the command

`./bin/hpcscan -v`

## Operating systems and compilers

hpcscan has been successfully tested on the operating systems and compilers listed below

Operating system | Compiler  | Remark
---------------- | --------- | ------
Ubuntu 20.04     | g++ 9.3.0 | -

# Validation

To check hpcscan has been correctly built and works fine, go to ./script and launch

`sh runValidationTests.sh`

This script runs a set a light test cases and should complete within few minutes (even on a laptop).

You should get in the ouptput report (displayed on the terminal)

* All tests marked as PASSED (370 tests passed per test mode enabled)
* No test marked as FAILED

Check the summary at the end of report to have a quick look on this.

# Execution

## Usage

hpcscan can be launched from a terminal with all configuration parameters within a single line.

**To get help on the parameters**

`./bin/hpcscan -h`

**Execution with a unique MPI process**

`mpirun -n 1 ./bin/hpcscan -testCase <TESTCASE> -testMode <TESTMODE>`

where
* TESTCASE is the name of the test case (see [List of test cases](#list-of-test-cases))
* TESTMODE is the name of the test mode (see [List of test modes](#list-of-test-modes))

Example

`mpirun -n 1 ./bin/hpcscan -testCase Propa -testMode CacheBlk`

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

TO DO

# Performance benchmarks

These tests are intended to measure various bandwidths

> **They are intensive tests that require to run on HPC platforms**

Performance measurements and scripts to reproduce results can be found in ./doc/TestCases/TestCases.pdf

# Versions

Version      | Description | Release date
------------ | ----------- | ------------
1.0          | Initial version with test modes Baseline, CacheBlk and NEC_SCA  | Nov 28, 2020
1.1          | On going work on test modes OpenAcc and Cuda | Coming soon

# Have fun!

Please share and send your feedbacks to vetienne@rocketmail.com


