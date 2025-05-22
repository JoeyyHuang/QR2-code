This QR2-code program is intended to calculation resonant Raman spectrum, including the single-resonant, double-resonant, as well as defect-induced double-resonant Raman, based on the DFT calculation. In current release, you need the 7.3.1 version of the Quantum ESPRESSO (QE) package and the EPW module therein for preliminary calculations to obtain the electron-photon and electron-phonon coupling matrices. The QR2-code can integrate the previous imformation and generate the Raman spectrascopy at any given Laser energy with any incident- and scattered light polarization configuration. In double resonant Raman case, the QR2-code can further analyze the contribution of each possible pair of phonons to the Raman peak intensity of interest and therefore give the assignments of double resonant Raman quantitatively which is a vital task in double resonant Raman studying.

NB: The current QR2-code is developed and tested based on QE of version 7.3.1. 
We strongly recommend using it with QE-7.3.1 for compatibility.

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

# Documentation
For further introduction, installing, input variables description, and detailed calculating steps, please refer the official web site: [http://qr2-code.com](http://qr2-code.com/).

If you encounter any issues while using it, please feel free to contact us using the GitHub discussions page: [https://github.com/JoeyyHuang/QR2-code/discussions](https://github.com/JoeyyHuang/QR2-code/discussions). We will get back to you as soon as possible!

# References and citing
The theory behind the QERaman code is described in our pre-print:
> J. Huang, R. Liu, Y. Zhang, N. T. Hung, H. Guo, R. Saito and T. Yang, [QR2-code: An open-source program for double resonance Raman spectra](https://arxiv.org/abs/2505.10041), *arXiv.2505.10041*

# Contact

Jianqi Huang

Liaoning Academy of Materials, Shenyang 110167, China

E-mail:&nbsp;&nbsp;jqhuang@lam.ln.cn

Nguyen Tuan Hung

Department of Physics, Tohoku University, Sendai 980-8578, Japan

E-mail:&nbsp;&nbsp;nguyen.tuan.hung.e4@tohoku.ac.jp

Teng Yang

Institute of Metal Research, Chinese Academy of Sciences, Shenyang 110016, China

E-mail:&nbsp;&nbsp;yanghaiteng@msn.com

# License
GNU General Public License (v3)