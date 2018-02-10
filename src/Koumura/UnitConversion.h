/*
 * (C)opyright 2013, by Kuroda-lab, The University of Tokyo, Japan.
 * http://kurodalab.bs.s.u-tokyo.ac.jp/
 * 
 * This file is a part of a program: "A stochastic model of Ca2+ Increase in Spines".
 * 
 * Licensed under the MIT License.
 * 
 * Please cite the following article if you use this program in your work:
 * Koumura T, Urakubo H, Ohashi K, Fujii M, Kuroda S. "Stochasticity in Ca2+ Increase in Spines Enables Robust and Sensitive Information Coding". PLoS One. 2014. e99040. doi:10.1371/journal.pone.0099040.
 */

#ifndef UNITCONVERSION_H
#define UNITCONVERSION_H

class UnitConversion
{
public:
	static const double AVOGADRO;

	static double um3_to_l(double value);
	static double num_to_umol(double num);
	static double umol_to_num(double umol);
	static double num_um3_to_uM(double num, double volume);
	static double uM_um3_to_num(double uM, double volume);
};

#endif
