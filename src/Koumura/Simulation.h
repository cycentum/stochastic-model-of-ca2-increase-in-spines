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

#ifndef SIMULATION_H
#define SIMULATION_H

#include "../General/Collection/Container2.h"
#include <list>
#include <vector>
#include <map>
#include <set>
using namespace std;

class Molecule;
class Reaction;
class Stimulus;
class Random;
class Result;

class Simulation
{
private:
	list<Molecule*> __molecule;
	list<Reaction*> __reaction;
	vector<Stimulus*> __stimulus;
	double __time;
	Result* __result;
	int __timeResultIndex, __stimulusIndex;
	bool __resultUpdated;
	vector<double> __timeResult;
	//Runge Kutta Fehlberg
	double __errorCoef, __dTimeMin, __dTimeMax, __dTime;
	static const double __K_WEIGHT_MATRIX[][4], __K_WEIGHT[], __K_WEIGHT_FOR_ERROR[];
	//Gillespie
	Random* __random;
	//tau leap
	double __tau1Min;
	bool __abandonTauLeap;
	void __restoreNumPrev();

	void __init();
	void __copyToNumPrev();
	bool __stepStim(double time);
	void __stepResult(double timeNext);
	void __updateMoleculeRungeKuttaFehlberg();

public:
	static const double NUM_K;
	Simulation(list<Molecule*> molecule, list<Reaction*> reaction, vector<Stimulus*> stimulus, vector<double> timeResult, double dtMin, double dtMax, double acceptableError);
	Simulation(list<Molecule*> molecule, list<Reaction*> reaction, vector<Stimulus*> stimulus, vector<double> timeResult, Random* random);
	Simulation(list<Molecule*> molecule, list<Reaction*> reaction, vector<Stimulus*> stimulus, vector<double> timeResult, Random* random, double acceptableError, double tau1Min, int numCritical);
	~Simulation();
	void dt(double dtMin, double dtMax, double acceptableError);
	void add(Molecule *molecule);
	void add(Reaction *reaction);
	void add(Stimulus *stimulus);
	void addTimeResult(double timeResult);
	
	void clearMolecule();
	void clearReaction();
	void clearStimulus();
	void clearTimeResult();
	
	void stepGillespie();
	void stepRungeKutta();
	void stepTauLeap();
	bool abandonTauLeap();

	bool finished();
	bool resultUpdated();
	Result* result();
	double time();
	void random(Random* random);
};

#endif
