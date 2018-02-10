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

#include "Molecule.h"
#include "UnitConversion.h"
#include "Simulation.h"
#include <iostream>
#include <cmath>

double errorScaleA(double num);
double errorScaleAA(double num);
double errorScaleAB(double num);
double errorScaleAAA(double num);
double errorScaleAAB(double num);
double errorScaleABC(double num);

double Molecule::num() {return __num;}
double Molecule::numPrev() {return __numPrev;}
String Molecule::name() {return __name;}
double Molecule::volume() {return __volume;}

void Molecule::fix(bool f){__fix=f;}

Molecule::Molecule(String na, double v, double nu) : __name(na), __volume(v), __num(nu), __fix(false){}

void Molecule::updateNum(double delta)
{
	if(!__fix) __num+=delta;
/*	if(__num<0)
	{
//		__num=0;
		cerr<<"Negative number: "<<name()<<endl;
		exit(1);
	}*/
}

void Molecule::copyToNumPrev() {__numPrev=__num;}

double Molecule::k(int index){return __k[index];}
void Molecule::initK()
{
	for(int k=0; k<Simulation::NUM_K; ++k) __k[k]=0;
}
void Molecule::addK(int index, double value){__k[index]+=value;}
void Molecule::multiplyDtToK(int index, double dt){__k[index]*=dt;}

double Molecule::kSum(const double weight[])
{
	double sum=0;
	for(int kIndex=0; kIndex<Simulation::NUM_K; ++kIndex) if(weight[kIndex]!=0) sum+=weight[kIndex]*k(kIndex);
	return sum;
}

void Molecule::updateNum(const double kWeight[]){updateNum(kSum(kWeight));}

void Molecule::restoreNumPrev(){__num=__numPrev;}

void Molecule::initMuSigmaError()
{
	__mu=0;
	__sigma2=0;
	__error=-1;
}

void Molecule::updateError(double error)
{
	if(__error<0||error<__error) __error=error;
}

bool Molecule::reactantOfNonCritical()
{
	return __error>0;
}

double Molecule::tau1(double errorCoef)
{
	double m=__error/abs(__mu);
	double v=__error*__error/__sigma2;
	return m<v?m:v;
}

void Molecule::addMuSigma(double mu, double sigma2)
{
	__mu+=mu;
	__sigma2+=sigma2;
}
