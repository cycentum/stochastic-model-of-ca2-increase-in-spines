# A stochastic model of Ca2+ Increase in Spines 

Source codes of the model in the manuscript: Koumura T, Urakubo H, Ohashi K, Fujii M, Kuroda S. "Stochasticity in Ca2+ Increase in Spines Enables Robust and Sensitive Information Coding". PLoS One. 2014. e99040. doi:10.1371/journal.pone.0099040

(C)opyright 2013, by Kuroda-lab, The University of Tokyo, Japan
Date: 30 March 2013
Also available in http://kurodalab.bs.s.u-tokyo.ac.jp/info/Ca_Increases/

## Reference
Please cite the following article if you use this program in your work.

Koumura T, Urakubo H, Ohashi K, Fujii M, Kuroda S. "Stochasticity in Ca2+ Increase in Spines Enables Robust and Sensitive Information Coding". PLoS One. 2014. e99040. doi:10.1371/journal.pone.0099040

## Files
+ README.md: this file.
+ src/*: source codes, including the makefile.

## Compiling
The model is written in C++. The make command is available. If you would like to copy the files to another directory, please do not change the relative location of the files in the src directory.

## Parameters
The parameters of the model can be manipulated in the main function in
src/Koumura/main.cpp.

+ volume: volume of the spine (0.1 um3), for stochastic model only.
+ PFCVInterval: PF-CF interval (sec).
+ ampPF: vector which contains the amplitude of each PF input, relative to the amplitude written in the manuscript. The number of the PF inputs corresponds to the size of the vector.
+ randSeed: random seed, for stochastic model only.
+ file: file name of the output file.

The runStochastic function and the runDeterministic function are for the running of the stochastic model and the deterministic model.

## Output
The output is the concentration (uM) of the all molecules in the model, sampled every 1 msec. The output is saved in the file specified by the "file" parameter in the main function.

## Libraries
This program includes dSFMT developed by Mutsuo Saito and Makoto Matsumoto, licensed under the new BSD License.

## License
MIT License.
