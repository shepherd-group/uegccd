# uegCCD
Code for calculations on the uniform electron gas

# authors
2024: William Z. Van Benschoten, James J. Shepherd

2023: William Z. Van Benschoten, Tina N. Mihm, James J. Shepherd

2018: Tom Henderson, James J. Shepherd, Gustavo Scuseria

# references
Please cite:

1. J. Chem. Phys. 154, 024113 (2021) 10.1063/5.0033408
2. J. Chem. Phys. 50, 191101 (2019) 10.1063/1.5091445
3. J. Chem. Phys. 140, 124102 (2014) 10.1063/1.4867783
4. Phys. Rev. Lett. 112, 133002 (2014) 10.1103/PhysRevLett.112.133002

# how to use: Mac
1. On a mac with xCode (& command line tools) installed and the GNU Fortran (GCC) 6.3.0 compiler, run the make script (./make) to compile. 
2. ZCode inside directory ZRun will execute with a file labelled Input in the same directory.
3. Output will be produced including MP2 and CCD calculations

On a different system, you need to modify the make script to include
a reference to Lapack. This code will compile with Lapack, but by
default it's just set up to compile on a mac.

# licence
MIT License. Copyright (c) 2018 Tom Henderson, James J. Shepherd, Gustavo Scuseria

MIT License. Copyright (c) 2023

MIT License. Copyright (c) 2024
