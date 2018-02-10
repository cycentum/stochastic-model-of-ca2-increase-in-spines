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

#ifndef MOLECULE_H
#define MOLECULE_H

#include "../General/General/String.h"

class Molecule
{
private:
	double __num, __numPrev, __volume, __k[6];
	String __name;
	bool __fix;

	double __mu, __sigma2, __error;
	
public:
	double num();
	double numPrev();
	String name();
	double volume();
	void fix(bool f);
	
	Molecule(String name, double volume, double num);
	
	void updateNum(double delta);
	void copyToNumPrev();

	double k(int index);
	void initK();
	void addK(int index, double value);
	void multiplyDtToK(int index, double dt);
	double kSum(const double weight[]);
	void updateNum(const double kWeight[]);

	void restoreNumPrev();

	//tau leap
	void initMuSigmaError();
	void updateError(double error);
	bool reactantOfNonCritical();
	double tau1(double errorCoef);
	void addMuSigma(double mu, double sigma2);
};

#endif
