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

#ifndef REACTION_H
#define REACTOIN_H

#include "../General/General/String.h"
#include "../General/Random/Random.h"
#include <list>
#include <map>
#include <set>

class Molecule;

class ReactionEntry
{
friend class Reaction;
private:
	Molecule* molecule;
	int num, delta;
	ReactionEntry(Molecule* molecule);

	//tau leap
	double (*errorScale)(double, double);
	double numCritical;
};

class Reaction
{
friend bool Reaction_compare(const Reaction* reac0, const Reaction* reac1);

private:
	list<ReactionEntry> __molecule;
	String __name;
	double __rateConstant, __flux, __stochasticity;

public:
	Reaction(String name, double rateConstant);
	
	void stochasticity(double stochasticity);
	double stochasticity();

	void addForward(Molecule* molecule, int num);
	
	void addBackward(Molecule* molecule, int num);
	
	double calcFlux();
	
	double flux(const double kWeight[]);
	double flux();
	String name();
	double rateConstant();
	void execute(double tau, Random* random);

	void execute();
	void setK(int index);

	set<Molecule*> molecule();

	//tau leap
	void setNumCritical(int n);
	bool critical();
	void setErrorScale();
	void updateError(double errorCoef);
	void addMuSigma();
};

Reaction* Reaction_new_from_conc_uM_volume_um3(String name, double rateConstant, map<Molecule*, int>& forward, map<Molecule*, int>& backward);

bool Reaction_compare(const Reaction* reac0, const Reaction* reac1);

#endif
