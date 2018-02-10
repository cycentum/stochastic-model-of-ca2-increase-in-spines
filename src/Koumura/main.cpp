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

#include "Model.h"
#include "Result.h"
#include "../General/Random/Random.h"
#include <fstream>

void runStochastic(double volume, double PFCVInterval, vector<double> ampPF, int randSeed, const char* file)
{
	Model model;
	model.setStim(PFCVInterval, ampPF);
	model.setStocahsticity(1.0/volume);
	model.tauLeap(true);
	Random random(randSeed);
	model.run(&random);
	ofstream ofs(file);
	model.result()->saveConc(ofs);
}

void runDeterministic(double PFCVInterval, vector<double> ampPF, const char* file)
{
	Model model;
	model.setStim(PFCVInterval, ampPF);
	model.setStocahsticity(0);
	model.tauLeap(true);
	model.run(NULL);
	ofstream ofs(file);
	model.result()->saveConc(ofs);
}

int main(int argc, char* argv[])
{
	double volume=1;
	double PFCVInterval=0.16;
	vector<double> ampPF;
	ampPF.push_back(1);
	ampPF.push_back(1);
	ampPF.push_back(1);
	ampPF.push_back(1);
	ampPF.push_back(1);
	int randSeed=0;
	char file[]="output.txt";

	runStochastic(volume, PFCVInterval, ampPF, randSeed, file);
	//runDeterministic(PFCVInterval, ampPF, file);
}
