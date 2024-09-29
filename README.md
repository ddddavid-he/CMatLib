# CMatLib

A C library of basic matrix operations.

Defining **basic**: Calculations are performed in the most trivial methods, nothing about acceleration is considered. 
And only `double` type is supported for contents of matrices.

## Usage

Usage and samples of each function are written in source files and header files.

Name of each function starts with a prefix `cml_` using camel naming system.

### difference between `basic.c` and `CMatLib.c`

 - `basic.c` and `basic.h` defines function for minimal usage. 
    Functions are defined without validation for input values. Incorrect values such as NULL or out-of-range index may 
    cause unexpected behaviours.
 
 - `CMatLib.c` and `CMatLib.h` defines functions for common usage.
   Functions are actually those in `basic.c` wrapped with input value validation logics.
   If incorrect values are passed to functions, an error message will pop up with debug messages.

### using CMatLib as C Library

For convenience, users can build the `CMatLib.c` or `basic.c` into a static library files using command:

 ```bash
  # library for CMatLib.c
  make libcml  
 ```
which generates file `lib/libcml.a`
 
 ```bash
  # library for basic.c
  # target in lib/libcml-basic.a
  make libcml-basic  
 ```
 generates file `lib/libcml-basic.a`.
 
After the above operations, users can use this library simply through `#include "basic.h"` or `#include "CMatLib.h"` in 
the source code and compile that with FLAGS including `-Llib -lcml-basic` or `-Llib -lcml`.

### embed the code into your own project

Do whatever you want to, tired of writing docs.

