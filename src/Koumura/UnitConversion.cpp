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

#include "UnitConversion.h"

const double UnitConversion::AVOGADRO=6e+23;

double UnitConversion::um3_to_l(double value)
{
	return value*1e-15;
}

double UnitConversion::num_to_umol(double num)
{
	return num/AVOGADRO*1e+6;
}

double UnitConversion::umol_to_num(double umol)
{
	return umol*AVOGADRO/1e+6;
}

double UnitConversion::num_um3_to_uM(double num, double volume)
{
	return num_to_umol(num)/um3_to_l(volume);
}

double UnitConversion::uM_um3_to_num(double uM, double volume)
{
	return umol_to_num(uM)*um3_to_l(volume);
}
