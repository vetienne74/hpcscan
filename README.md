
===========================================================
=                                                         =
=                  WELCOME TO HPCSCAN!                    =
=                                                         =
=   Version 1.0                                           =
=   Contact: Vincent Etienne                              =
=            Saudi Aramco / EXPEC ARC / GPT               =
=            vetienne@rocketmail.com                      =
=                                                         =
===========================================================

# Description

HPCSCAN is an HPC Test Cases Suite for Benchmarking Algorithms and Computing Platforms

# Set up your environment

You just need a C++ compiler with MPI & OpenMP

# Compile

Go to ./build, and use the command "make"
The Makefile adapts to the machine you are connected
Have a look in the Makefile if you need to add a new machine

By default, hpcscan is compiled in single in precision
To compile in double precision "make precision=double"

Executable is ./bin/hpcscan

# Validate

To check good the behavior of hpcscan
Go to ./script
And "sh runValidationTests.sh"
This script runs a set a light test cases within a short time (even on a laptop)

You should get in the ouptput report (displayed on the terminal)
All tests marked as PASSED
and zero test marked as FAILED
Check the summary at the end of report to have a quick look on this

# Run a specific test case

All parameters are placed on the command line

## To get help on the parameters
./bin/hpcscan -h

## Example to run the Comm test case

mpirun ./bin/hpcscan -testCase Comm

# Run performance test cases

These tests are intended to measure various bandwidths
They are intensive tests that require to run on HPC platforms 

To run all test cases
Go to ./script

TO BE UPDATED...

===========================================================
=                  H A V E    F U N !                     =
===========================================================
