This is the distribution of the QR^2-code package (QR^2: Quantum-theory-based Resonant Raman).
It is a driver to compute the Resonant Raman spectrum including single, double and defect-induced double resonant Raman, 
using Quantum ESPRESSO (QE) package and EPW module therein as the underlying engine.

NB: The current QR^2-code is developed and tested based on QE of version 7.3.1. 
We strongly recommend using it with QE-7.3.1 for compatibility.

The QR^2 package consists of two folders, "examples" and "source code", in the current release v1.0.0
(https://github.com/JoeyyHuang/QR2-code/archive/refs/tags/v1.0.0.tar.gz).

## PACKAGES
```md
QR2-code
├── examples
│   ├── Graphene: input files for single, double and defect resonant Raman calculations of semimetal material
│   ├── MoS2: input files for single and double resonant Raman calculations of semiconductor material with spin-orbit coupling under consideration
│   └── hBN: input files for single and double resonant Raman calculations of insulator material
└── source code
    ├── EPW-src: the core codes for Raman spectrum calculating, modified based on the original files with the same names in QE/EPW/src
    ├── PHonon-PH: phonon-related code for sorting phonon branches in terms of symmetry, modified based on the original file with the same name in QE/PHonon/PH
    ├── Raman_PP: post-processing code for Raman spectrum plot and analysis
    └── phonon_sort: script for generating the input files to perform the phonon-sorting calculation
```

examples: 
--Graphene: input files for single, double and defect resonant Raman calculations of semimetal material
--MoS2: input files for single and double resonant Raman calculations of semiconductor material with spin-orbit coupling under consideration
--hBN: input files for single and double resonant Raman calculations of insulator material
source code:
--EPW-src: the core codes for Raman spectrum calculating, modified based on the original files with the same names in QE/EPW/src
--PHonon-PH: phonon-related code for sorting phonon branches in terms of symmetry, modified based on the original file with the same name in QE/PHonon/PH
--Raman_PP: post-processing code for Raman spectrum plot and analysis
--phonon_sort: script for generating the input files to perform the phonon-sorting calculation

## USAGE
1. put the source code files in EPW-src inside QE/EPW/src directory and repalce the original files with the same names;
    put the source code files in PHonon-PH inside QE/PHonon/PH directory and repalce the original file with the same name;
    
    cd to the main QE directory and type
    ```
    ./configure [options]
    make pw ph
    make epw
    ```
   Optionally, `make -jN pw ph` and `make -jN epw` runs parallel compilation on `N` processors.

2. put the source code files in Raman_PP and phonon_sort inside QE/bin directory;

    cd to the QE/bin directory and type
    ```
    gfortran phonon_sort.f90 -o phonon_sort.x
    mpiifort raman_pp.f90 -o raman_pp.x
    ```
Now you should have the executables (pw.x, ph.x, q2r.x, matdyn.x, phonon_sort.x, epw.x and raman_pp.x) in QE/bin directory required for resonant Raman calculations
and you can run the examples.

For more information about the relevant theory and variables explanation in the calculations, please see the web site [http://qr2-code.com](http://qr2-code.com/).

# References and citing
The theory behind the QERaman code is described in our pre-print:
> J. Huang, R. Liu, Y. Zhang, N. T. Hung, H. Guo, R. Saito and T. Yang, [QR2-code: An open-source program for double resonance Raman spectra](https://arxiv.org/abs/2505.10041), *arXiv.2505.10041*

# License
GNU General Public License (v3)