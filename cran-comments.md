## Test environments
* local ubuntu 20.04 install, R 4.0.5
* Debian Linux, R-devel, GCC ASAN/UBSAN
* Fedora Linux, R-devel, clang, gfortran
* win-builder (devel and release)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* check_win_devel()
* check_rhub()

## R CMD check results
There were no ERRORs, or WARNINGs. There was one NOTE, because an optional dependency on package GreedyEPL was recently removed from CRAN. I am contacting the maintainer to see if the package could be put back online, but as it is an optional dependency, it should not prevent the release of BNPdensity. 

## Downstream dependencies
There are currently no downstream dependencies for this package.