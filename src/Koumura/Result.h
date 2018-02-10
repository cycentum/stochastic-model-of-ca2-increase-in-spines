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

#ifndef RESULT_H
#define RESULT_H

//#include "../General/IO/BinaryWriter.h"
#include <map>
#include <list>
#include <vector>
#include <fstream>
using namespace std;

class Molecule;

class Result
{
private:
	list<list<double> > __num;
	list<double> __time;
	list<Molecule*> __molecule;

public:
	Result(list<Molecule*> molecule);
//	void add(double timePrev, double timeNext, double timeResult, list<Molecule*> molecule);
	void add(double time, bool numPrev);
//	void saveConc(BinaryWriter& writer);
//	void saveNum(BinaryWriter& writer);
	void saveConc(ofstream& ofs);
	void saveNum(ofstream& ofs);
	double timeLast();
	int size();
};

#endif
