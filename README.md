<p align="left">
<a href="https://github.com/bkosanovic/fcmclt/blob/master/LICENSE.md"><img alt="MIT License" src="https://img.shields.io/github/license/bkosanovic/fcmclt" /></a>
<a href="https://github.com/bkosanovic/fcmclt/releases/latest"><img alt="GitHub Release" src="https://img.shields.io/github/v/release/bkosanovic/fcmclt?display_name=tag" /></a>
<a href="https://github.com/bkosanovic/fcmclt/releases"><img alt="GitHub All Releases" src="https://img.shields.io/github/downloads/bkosanovic/fcmclt/total" /></a>
<a href="https://www.mathworks.com/matlabcentral/fileexchange/120333-fcmclt"><img alt="View fcmclt on MATLAB File Exchange" src="https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg" /></a>
</p>

# fcmclt
Fuzzy clustering C implementation for MATLAB (FCM, Gustafson-Kessel, clustering validity,
extrapolation with presumed cluster centers)

## Project status

The only planned changes at this time are to improve the documentation, i.e. this README file.

## Description

This repository provides efficient C implementation of two types of Fuzzy c-Means Clustering
algorithms and
several related tools that can be used in clustering evaluation and validation. Once installed,
the provided
functions are directly callable within MATLAB. To use them outside of MATLAB, one would need to
modify the existing C code.
Additional implementation details that may be useful are included
in the [Implementation section](/README.md#Implementation) below.

The following functionality is supported:
- [Fuzzy c-Means Clustering](https://en.wikipedia.org/wiki/Fuzzy_clustering)
  (based on the work of _James C. Bezdek_ and _Enrique Ruspini_.)
  - Euclidean, diagonal, and Mahalanobis distance metric
- **Gustafson-Kessel** fuzzy c-means clustering
- **Initial fuzzy partition matrix** generation
- **12 Clustering validity functionals:**
  - partition coefficient, partition entropy, nonfuzzy index,
    minimum and mean hard tendencies, minimum and maximum relative fuzziness,
    the minimum nearest maximum membership cardinality, compactness
    and separation index, fuzzy hypervolume, average partition
    density, and the partition density of a resulting fuzzy partition
- **Fuzzy partition matrix extrapolation** that uses the presumed cluster centers
- **Fuzzy scatter and covariance matrices** calculation for fuzzy clusters
- **Test code** (in MATLAB) to verify the proper installation of software
- **Detailed help text** describing how to use the provided functions

The primary motivation for this work was to estimate membership functions of the
[Temporal Fuzzy Sets](https://ieeexplore.ieee.org/abstract/document/375073)
which span the
[Fuzzy Information Space](https://www.academia.edu/39731003/SIGNAL_AND_SYSTEM_ANALYSIS_IN_FUZZY_INFORMATION_SPACE).
Nevertheless, the provided software package may be used for solving any problem that
requires fuzzy c-means clustering, e.g. tracking of speech and noise power changes,
echo path estimation, analysis, and characterization, etc. Fuzzy clustering has been applied
within a wide range of engineering and research disciplines including geology, agriculture,
social sciense, etc.

## Installation

The installation may be done either by downloading the latest release which includes generated
binaries for the Windows x64 MEX files (MATLAB 2022b used), or by cloaning or downloading
the source code. The process of generating your own binaries from the source code
is straighforward.

### Prerequisites

To use this 
software one must have [MATLAB](https://www.mathworks.com/) installed and if you plan
to generate the MEX binaries you must be able to run the
[mex command.](https://www.mathworks.com/help/matlab/matlab_external/build-an-executable-mex-file.html) 
This command is used for generating
the accelerated versions of the three main algorithms from their C source code.

### Steps to install

#### From the release

Once you verified you have MATLAB installed you can download the latest release. All the
necessary files are contained within the `fcmclt_M_m_p_w64.zip` file. When you unpack the
zip archive, the MATLAB scripts and the binaries will be contained in the
`./src` folder. You may optionally copy the contents of ./src directory to a folder where you
want MATLAB
to access it. As an example, let us assume that you copied all the files
from `./src` to `D:\MyMatlab\fcmclt` on a Windows system.

You may now proceed to testing your installation [as described below](#testing-and-finalizing-your-installation).

#### Manual from the source with generating MEX binaries

Once you verified you have MATLAB installed and you can successfully run mex command to build
MATLAB functions implemented in C language, you may proceed with the following steps:

- Clone or download this repository to the system where you have MATLAB installed
- Copy the contents of `./src` directory to a folder where you want MATLAB to access it
  - As an example, let us assume you copied all the files from `./src` to `D:\MyMatlab\fcmclt`
    on a Windows system.
- Start MATLAB and within MATLAB change directory to `D:\MyMatlab\fcmclt`
  - `>> cd D:\MyMatlab\fcmclt`
- Run the following three `mex` commands:
  ```matlab
  mex extfpm.c
  mex fcmc.c
  mex gkfcmc.c
  ```
- You should see three new files created: `extfpm.<mexext>`, `fcmc.<mexext>`,
  `gkfcmc.<mexext>`, where `<mexext>` reflects the mex file extension on your system. For
  64-bit Windows system, it should be `mexw64`.

#### Testing and finalizing your installation

- You can test your installation by running test routines that were provided in `./src` folder
  and that you copied to `D:\MyMatlab\fcmclt`,
   e.g. try running `>> testfcmc` from that folder.
- Try all other test routines. They have prefix `test` or `tst`.
  - NOTE: Some may display a warning indicating the dataset may be too small for reliable
    clustering. That is normal and expected behavior for those examples.

If you want to make the `fcmclt` package available to you no matter which folder you are in,
you may
add the `D:\MyMatlab\fcmclt` to your MATLAB path. You may consult MATLAB documentation
if you do not know how to do it.

## Usage

The best way to learn how to use the functions provided in this package is to read their
documentation. You may run `>> help <function-name>` to get the function description.
In addition to the help text for the three main algorithms (fcmc, gkfcmc, and extfpm),
you should consult the MATLAB source code for the test routines provided in `.m` files. That
will provide you with examples of how to integrate these routines into your programs
and how to initialize
or evaluate the results.

The MATLAB source code may also provide references to additional reading material. The notation
used in the source code is based primarily on this book:

- **J.C. Bezdek,** _"Pattern Recognition with Fuzzy Objective
              Function Algorithms,"_ Plenum Press, New York, 1981.

## Support

The support is not provided for this package.

## Implementation

The initial implementation of this package was done between 1992-1995. There were only a handful
of minor changes done since that time. Those were limited to accomodate the MEX API 
changes introduced by [MathWorks, Inc.](https://www.mathworks.com/)
over a long period of time (25+ years).

The code in `./src` folder contains C and MATLAB sources. The C files are:

- `extfpm.c`: Extrapolates the fuzzy partition matrix using the presumed cluster centers
- `fcmc.c`: Implements the Fuzzy c-Means Clustering (FCMC) algorithm
- `gkfcmc.c`: Implements Gustafson-Kessel (GK) variant of the FCMC algorithm

The relevant MATLAB files are:

- `Contents.m`: Contains help text for the `fcmclt` package.
- `fcmcinit.m`: Very important routine that is used to generate the initial fuzzy partition
   matrix U<sub>0</sub>.
- `cltvalid.m`: Calculates clustering validity functionals
- `fscat.m`: Calculates fuzzy scatter and covariance matrices for a fuzzy cluster
- `testextfpm.m`: Example of how to extrapolate fuzzy partition matrix from the presumed
   cluster centers
- `testfcmc.m`: Example of how to use FCMC algorithm
- `testgk.m`: Example of GK variant for two Gaussian classes
- `testgk2.m`: Example of GK for Gustafson's cross
- `tstvalid.m`: Example of how to use validity functionals for the FCMC algorithm

### C code organization

The C functions include only two header files. The `math.h` and `mex.h`. The only
two core MATLAB algorithms that are invoked within the C code are the matrix inverse and
determinant.

For brevity, only the `fcmc.c` code structure is outlined here:

- The main function called by MATLAB is `mexFunction()`
  - Checks the input arguments and applies the defaults for optional ones as needed
  - Fetches pointers to vector data
  - Allocates some memory
  - Prepares for the FCMC routine
    - Selects the distance metric (may allocate additional memory and may invoke
      MATLAB inverse function in case of Mahalanobis metric)
  - Executes the FCMC algorithm (`do_fcm()`)
  - Checks and creates output variables
- All other functions are written in plain C code

In essense, the `mexFunction()` is a MATLAB wrapper for the algorithms that are to be
accelerated.

NOTE: If you would like to use the C code outside of MATLAB, you may need to do a few things:

- Make it re-entrant. That is, you would need to instantiate the variables that are currently
  global or static. That was not a problem for the MATLAB implementation, but may create
  problems if you try to integrate this code somewhere else.
- You would need to provide compatible
  implementations for matrix inverse and determinant MATLAB functions in case you have
  to use the parts of code that depend on these.
- Finally, you would need to replace the `mexFunction()` implementation with your own
  wrapper that would handle memory management, optional input arguments, etc.
  The memory management you will need to implement will likely differ from the
  way MATLAB manages memory. Especially when it comes to allocating and deallocating
  memory.

The code in `gkfcmc.c` and `extfpm.c` files is structured in a similar way.

## Roadmap

There are no plans to expand this package with additional features.

## Contributing

Pull requests will only be considered for the following contributions:

- Bug fixes (if you can find any)
- Interesting examples that show how to use the provided functions
  - Please try to provide a single self-contained MATLAB file that does not require any
    data files. If additional MATLAB toolboxes are required, please list them all.
  - ***Document your code!*** If your code cannot be easily reviewed
    it will be swiftly rejected no matter how great you think it is.

## Authors

The software provided here has been developed by Bogdan Kosanovic in the early 1990s during his 
Ph.D. work at the University of Pittsburgh.

## Acknowledgments

The author would like to acknowledge _Dr. Brahim Hamadicharef_
who was at the time
(August, 2003) with the
University of Plymouth for submitting a fix to make the `fcmclt` package work for
MATLAB 6.5 R13. If we are to trust the LinkedIn, Dr. Hamadicharef is now (2022) with
the _Institute of High Performance Computing (IHPC)_ in Singapore.

## References

For general information on fuzzy clustering, partition coefficient (F), and partition
entropy (H), the best starting point is

    J.C. Bezdek, "Pattern Recognition with Fuzzy Objective
    Function Algorithms", Plenum Press, New York, 1981.

For more information on nonfuzzy index (NFI), please refer to

    M. Roubens, "Pattern Classification Problems and Fuzzy Sets",
    Fuzzy Sets and Systems, 1:239-253, 1978.

For more information on minimum and mean hard tendencies (MinHT and MeanHT), please refer to

    F.F. Rivera, E.L. Zapata, and J.M. Carazo, "Cluster validity
    based on the hard tendency of the fuzzy classification",
    Pattern Recognition Letters, 11:7-12, 1990.

For more information on minimum and maximum relative fuzziness (MinRF and MaxRF), please
refer to

    H.L. Gordon and R.L. Somorjai, "Fuzzy Cluster Analysis of
    Molecular Dynamics Trajectories", Proteins: Structure,
    Function, and Genetics, 14:249-264, 1992.

For more information on minimum nearest maximum membership cardinality (MinNMMcard)
and related functionals including a way to go about the selection of fuzzy
exponent $m$, please refer to

    B.R. Kosanovic, "Signal and System Analysis in Fuzzy Information Space",
    Ph.D. Dissertation, University of Pittaburgh, 1995.

For more information on compactness and separation index (S), please refer to

    X.L. Xie and G. Beni, "A Validity Measure for Fuzzy Clustering",
    IEEE Trans. PAMI, 13(8):841-847, 1991.

For more information on fuzzy hypervolume (Fhv), average partition
         density (Dpa), and partition density (Pd), please refer to

    I. Gath and A.B. Geva, "Fuzzy clustering for the estimation of
    the parameters of the components of mixtures of normal
    distributions", Pattern Recognition Letters, 9:77-86, 1989.

## License

[MIT License](/LICENSE.md)
