# uegCCD
Code for calculations on the uniform electron gas

## disclaimer
This is a beta release

## authors
2024: William Z. Van Benschoten, James J. Shepherd

2023: William Z. Van Benschoten, Tina N. Mihm, James J. Shepherd

2018: Tom Henderson, James J. Shepherd, Gustavo Scuseria

## references
Please cite:

1. J. Chem. Phys. 154, 024113 (2021) 10.1063/5.0033408
2. J. Chem. Phys. 50, 191101 (2019) 10.1063/1.5091445
3. J. Chem. Phys. 140, 124102 (2014) 10.1063/1.4867783
4. Phys. Rev. Lett. 112, 133002 (2014) 10.1103/PhysRevLett.112.133002

## compiling
uegCCD is built using a script in the `source` directory:

```bash
cd source
./make.sh
```

The script is set up to auto-detect if you are compiling on Mac or Linux.
If compiling on Linux, it will attempt to link against the OpenBLAS library for LAPACK support.
You may need to adjust this command if using a different LAPACK implementation.

By default, uegCCD will compile at the highest optimization.
The `make.sh` script supports a number of (combinable) flags to alter this behavior:

| Flag | Description |
|------|-------------|
| `-d` | Compile with debug symbols (`-Og -g`) |
| `-p` | Compile for profiling (`-pg`) |
| `-t` | Compile with OpenMP threading (`-fopenmp`) |

## licence
MIT License. Copyright (c) 2018 Tom Henderson, James J. Shepherd, Gustavo Scuseria

MIT License. Copyright (c) 2023

MIT License. Copyright (c) 2024
