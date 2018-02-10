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

#ifndef MODEL_H
#define MODEL_H

#include "../General/Random/Random.h"
#include "../General/General/String.h"
#include <vector>
#include <map>
#include <list>

class Simulation;
class Stimulus;
class Result;
class Molecule;
class Reaction;

class Model
{
private:
	map<String, Molecule*> __mole;
	map<String, Reaction*> __reac;
	vector<Stimulus*> __stim;
	list<Molecule*> __moleList;
	Simulation* __simulation;
//	double __clockS;
	bool __tauLeap;

	void __loadMole();
	void __loadReac();

public:
	Model();
	~Model();

	void setStim(double timing, double ampPFCa, double ampPFGlu, int numPF);
	void setStim(double timing, vector<double> amp);
	void setStimGlu(double numPerMS, double durMS);
	void setStimTest();
	void setStocahsticity(double s);

	void run(Random* random);
	Result* result();
//	void saveConc(int trial, String path);
	void tauLeap(bool tl);
};

#endif
