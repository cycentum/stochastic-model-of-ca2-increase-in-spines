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

#include "Reaction.h"
#include "Molecule.h"
#include "Simulation.h"
#include "UnitConversion.h"
#include <iostream>

double errorScaleA(double num, double stochasticity);
double errorScaleAA(double num, double stochasticity);
double errorScaleAB(double num, double stochasticity);
double errorScaleAAA(double num, double stochasticity);
double errorScaleAAB(double num, double stochasticity);
double errorScaleABC(double num, double stochasticity);

Reaction::Reaction(String n, double rc) : __name(n), __rateConstant(rc), __stochasticity(1) {}

void Reaction::stochasticity(double s)
{
	__stochasticity=s;
}
double Reaction::stochasticity() {return __stochasticity;}

void Reaction::addForward(Molecule* molecule, int num)
{
	bool found=false;
	for(list<ReactionEntry>::iterator m=__molecule.begin(); m!=__molecule.end(); ++m)
	{
		if(m->molecule==molecule)
		{
			m->num+=num;
			m->delta-=num;
			found=true;
			break;
		}
	}
	if(!found)
	{
		__molecule.push_back(ReactionEntry(molecule));
		__molecule.back().num=num;
		__molecule.back().delta=-num;
	}	
}

void Reaction::addBackward(Molecule* molecule, int num)
{
	bool found=false;
	for(list<ReactionEntry>::iterator m=__molecule.begin(); m!=__molecule.end(); ++m)
	{
		if(m->molecule==molecule)
		{
			m->delta+=num;
			found=true;
			break;
		}
	}
	if(!found)
	{
		__molecule.push_back(ReactionEntry(molecule));
		__molecule.back().num=0;
		__molecule.back().delta=num;
	}
}

double Reaction::calcFlux()
{
	if(__stochasticity==0 || __rateConstant==0)
	{
		__flux=0;
		return 0;
	}
	for(list<ReactionEntry>::iterator it=__molecule.begin(); it!=__molecule.end(); ++it)
	{
		if(it->molecule->num()<it->num*__stochasticity)
		{
			__flux=0;
			return 0;
		}
	}
	__flux=__rateConstant/__stochasticity;
	for(list<ReactionEntry>::iterator it=__molecule.begin(); it!=__molecule.end(); ++it)
	{
		Molecule* m=it->molecule;
		int num=it->num;
		for(int n=0; n<num; ++n) __flux*=m->num()-n*__stochasticity;
	}
	return __flux;
}
double Reaction::flux(){return __flux;}

double Reaction::flux(const double kWeight[])
{
	__flux=__rateConstant;
	if(__rateConstant==0) return 0;
	for(list<ReactionEntry>::iterator it=__molecule.begin(); it!=__molecule.end(); ++it)
	{
		Molecule* m=it->molecule;
		int num=it->num;
		for(int n=0; n<num; ++n) __flux*=m->num()+m->kSum(kWeight);
	}
	return __flux;
}

String Reaction::name() {return __name;}
double Reaction::rateConstant(){return __rateConstant;}

void Reaction::execute()
{
	for(list<ReactionEntry>::iterator it=__molecule.begin(); it!=__molecule.end(); ++it) it->molecule->updateNum(it->delta*__stochasticity);
}

void Reaction::setK(int index)
{
	for(list<ReactionEntry>::iterator it=__molecule.begin(); it!=__molecule.end(); ++it) it->molecule->addK(index, __flux*it->delta);
}

set<Molecule*> Reaction::molecule()
{
	set<Molecule*> mole;
	for(list<ReactionEntry>::iterator e=__molecule.begin(); e!=__molecule.end(); ++e) mole.insert(e->molecule);
	return mole;
}

void Reaction::setNumCritical(int n)
{
	for(list<ReactionEntry>::iterator it=__molecule.begin(); it!=__molecule.end(); ++it) if(it->delta<0)
	{
		it->numCritical=-it->delta*n*__stochasticity;
	}
}

bool Reaction::critical()
{
	for(list<ReactionEntry>::iterator it=__molecule.begin(); it!=__molecule.end(); ++it) if(it->delta<0)
	{
		if(it->molecule->num()<=it->numCritical) return true;
	}
	return false;
}

void Reaction::setErrorScale()
{
	int order=0;
	for(list<ReactionEntry>::iterator it=__molecule.begin(); it!=__molecule.end(); ++it) if(it->delta<0) order-=it->delta;
	for(list<ReactionEntry>::iterator it=__molecule.begin(); it!=__molecule.end(); ++it) if(it->delta<0)
	{
		switch(order)
		{
		case 1:
			it->errorScale=errorScaleA;
			break;
		case 2:
			if(it->delta==-1) it->errorScale=errorScaleAB;
			else it->errorScale=errorScaleAA;
			break;
		case 3:
			if(it->delta==-1) it->errorScale=errorScaleABC;
			else if(it->delta==-2) it->errorScale=errorScaleAAB;
			else it->errorScale=errorScaleAAA;
			break;
		}
	}
}

void Reaction::updateError(double errorCoef)
{
	for(list<ReactionEntry>::iterator it=__molecule.begin(); it!=__molecule.end(); ++it) if(it->delta<0)
	{
		double error=errorCoef*it->molecule->num()/it->errorScale(it->molecule->num(), __stochasticity);
		if(error<__stochasticity) error=__stochasticity;
		it->molecule->updateError(error);
	}
}

void Reaction::addMuSigma()
{
	for(list<ReactionEntry>::iterator it=__molecule.begin(); it!=__molecule.end(); ++it) if(it->delta<0)
	{
		double d=it->delta*__stochasticity;
		it->molecule->addMuSigma(d*flux(), d*d*flux());
	}
}

void Reaction::execute(double tau, Random* random)
{
	int num=random->getPoisson(flux()*tau);
	for(list<ReactionEntry>::iterator it=__molecule.begin(); it!=__molecule.end(); ++it) it->molecule->updateNum(it->delta*__stochasticity*num);
}

Reaction* Reaction_new_from_conc_uM_volume_um3(String name, double rateConstant, map<Molecule*, int>& forward, map<Molecule*, int>& backward)
{
	int degree=0;
	double volumeMin=0;
	for(map<Molecule*, int>::iterator it=forward.begin(); it!=forward.end(); ++it)
	{
		double v=it->first->volume();
		int d=it->second;
		for(int i=0; i<d; ++i) rateConstant/=v;
		if(volumeMin==0 || v<volumeMin) volumeMin=v;
		degree+=d;
	}
	for(int d=0; d<degree-1; ++d) rateConstant=UnitConversion::num_to_umol(rateConstant);
	
	Reaction* reaction=new Reaction(name, rateConstant*volumeMin);
	for(map<Molecule*, int>::iterator it=forward.begin(); it!=forward.end(); ++it) reaction->addForward(it->first, it->second);
	for(map<Molecule*, int>::iterator it=backward.begin(); it!=backward.end(); ++it) reaction->addBackward(it->first, it->second);
	return reaction;
}


ReactionEntry::ReactionEntry(Molecule* m) : molecule(m) {}


bool Reaction_compare(const Reaction* reac0, const Reaction* reac1)
{
	return reac0->__flux < reac1->__flux;
}


double errorScaleA(double num, double stochasticity){return 1;}
double errorScaleAA(double num, double stochasticity){return 2+stochasticity/(num-stochasticity);}
double errorScaleAB(double num, double stochasticity){return 2;}
double errorScaleAAA(double num, double stochasticity){return 3+stochasticity/(num-stochasticity)+2*stochasticity/(num-2*stochasticity);}
double errorScaleAAB(double num, double stochasticity){return 3+3*stochasticity/2/(num-stochasticity);}
double errorScaleABC(double num, double stochasticity){return 3;}
