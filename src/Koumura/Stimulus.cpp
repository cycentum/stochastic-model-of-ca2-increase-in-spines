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

#include "Stimulus.h"
#include "Molecule.h"
#include <iostream>

Stimulus::Stimulus(Molecule* m, double tb, double d, int n) : __molecule(m), __timeBegin(tb), __num(n)
{
	__numPerTime=__num/d;
	__numFinished=0;
}

bool Stimulus::step(double time)
{
	double n=__numPerTime*(time-__timeBegin);
	if(__stochasticity>0) n=((int)(n/__stochasticity))*__stochasticity;
	if(n>__num) n=__num;
	if(n-__numFinished>0)
	{
		__molecule->updateNum(n-__numFinished);
		__numFinished=n;
		return true;
	}
	return false;
}

double Stimulus::timeBegin() {return __timeBegin;}

void Stimulus::stochasticity(double s){__stochasticity=s;}

bool Stimulus::finished(){return __numFinished>=__num;}

bool Stimulus_compare(const Stimulus* stimulus0, const Stimulus* stimulus1)
{
	return stimulus0->__timeBegin < stimulus1->__timeBegin;
}
