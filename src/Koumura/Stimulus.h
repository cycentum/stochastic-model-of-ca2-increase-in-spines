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

#ifndef STIMULUS_H
#define STIMULUS_H

class Molecule;

class Stimulus
{
friend bool Stimulus_compare(const Stimulus* stimulus0, const Stimulus* stimulus1);

private:
	Molecule* __molecule;
	double __timeBegin, __numPerTime, __stochasticity, __numFinished;
	int __num;

public:
	Stimulus(Molecule* molecule, double timeBegin, double duration, int num);
	bool step(double time);
	double timeBegin();
	void stochasticity(double s);
	bool finished();
};

bool Stimulus_compare(const Stimulus* stimulus0, const Stimulus* stimulus1);

#endif
