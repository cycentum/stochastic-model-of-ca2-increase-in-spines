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

#include "Simulation.h"
#include "Molecule.h"
#include "Reaction.h"
#include "Stimulus.h"
#include "Result.h"
#include "UnitConversion.h"
#include "../General/Random/Random.h"
#include <cmath>
#include <algorithm>
#include <iostream>

/*const double Simulation::NUM_K=6;
const double Simulation::__K_WEIGHT_MATRIX[][5]=
	{{0, 0, 0, 0, 0},
	{1.0/4, 0, 0, 0, 0},
	{3.0/32, 9.0/32, 0, 0, 0},
	{1932.0/2197, -7200.0/2197, 7296.0/2197, 0, 0},
	{439.0/216, 8, 3680.0/513, 845.0/4104, 0},
	{-8.0/27, 2, -3544.0/2565, 1859.0/4104, -11.0/40}};
const double Simulation::__K_WEIGHT[]={25.0/216, 0, 1408.0/2565, 2197.0/4104, -1.0/5, 0};
const double Simulation::__K_WEIGHT_FOR_ERROR[]={1.0/360, 0, -128.0/4275, -2197.0/75240, 1.0/50, 2.0/55};*/

const double Simulation::NUM_K=4;
const double Simulation::__K_WEIGHT_MATRIX[][4]=
	{{0, 0, 0},
	{1.0/2, 0, 0},
	{0, 3.0/4, 0},
	{2.0/9, 1.0/3, 4.0/9}};
const double Simulation::__K_WEIGHT[]={7.0/24, 1.0/4, 1.0/3, 1.0/8};
const double Simulation::__K_WEIGHT_FOR_ERROR[]={5.0/72, -1.0/12, -1.0/9, 1.0/8};

Simulation::Simulation(list<Molecule*> m, list<Reaction*> r, vector<Stimulus*> s, vector<double> tr, double dtmi, double dtma, double acceptableError) : __molecule(m), __reaction(r), __stimulus(s), __dTimeMin(dtmi), __dTimeMax(dtma), __timeResult(tr)
{
	__errorCoef=pow(acceptableError/2, 1.0/4);
	__dTime=__dTimeMin;
	__init();
}
Simulation::Simulation(list<Molecule*> m, list<Reaction*> r, vector<Stimulus*> s, vector<double> tr, Random* ra) : __molecule(m), __reaction(r), __stimulus(s), __random(ra), __timeResult(tr)
{
	__init();
}
Simulation::Simulation(list<Molecule*> m, list<Reaction*> r, vector<Stimulus*> s, vector<double> tr, Random* ra, double acceptableError, double t1Min, int numCritical) : __molecule(m), __reaction(r), __stimulus(s), __random(ra), __timeResult(tr), __errorCoef(acceptableError), __tau1Min(t1Min)
{
	for(list<Reaction*>::iterator re=__reaction.begin(); re!=__reaction.end(); ++re)
	{
		(*re)->setNumCritical(numCritical);
		(*re)->setErrorScale();
	}
	__init();
}

Simulation::~Simulation()
{
	delete __result;
}

void Simulation::__init()
{
	__result=new Result(__molecule);
	__timeResultIndex=0;
	__stimulusIndex=0;
	__time=0;
	sort(__stimulus.begin(), __stimulus.end(), Stimulus_compare);
	__resultUpdated=false;
}

void Simulation::stepGillespie()
{
	__copyToNumPrev();

	double fluxAll=0;
	for(list<Reaction*>::iterator r=__reaction.begin(); r!=__reaction.end(); ++r) fluxAll+=(*r)->calcFlux();

	double timeNext;
	if(fluxAll==0)
	{
		//nothing will change if continued
		if(__stimulusIndex==__stimulus.size()) timeNext=__result->timeLast();
		//go to next stimulus
		else timeNext=__stimulus[__stimulusIndex]->timeBegin();
	}
	else
	{
		__reaction.sort(Reaction_compare);
		Reaction* reactionNext;
		double rand=__random->getDouble()*fluxAll;
		double cumFlux=0;
		for(list<Reaction*>::iterator r=__reaction.begin(); r!=__reaction.end(); ++r)
		{
			if((*r)->flux()==0) continue;
			cumFlux+=(*r)->flux();
			if(cumFlux>=rand)
			{
				reactionNext=*r;
				break;
			}
		}
		if(reactionNext==NULL) reactionNext=*(__reaction.rbegin());
		reactionNext->execute();
		timeNext=__time-1.0/fluxAll*log(__random->getDouble());
	}

	__stepStim(timeNext);
	__stepResult(timeNext);

	__time=timeNext;
}

void Simulation::stepTauLeap()
{
	__copyToNumPrev();
	__abandonTauLeap=false;

	list<Reaction*> critical, nonCritical;
	double timeNext;
	for(list<Reaction*>::iterator r=__reaction.begin(); r!=__reaction.end(); ++r)
	{
		(*r)->calcFlux();
		if((*r)->flux()==0) continue;
		if((*r)->critical()) critical.push_back(*r);
		else nonCritical.push_back(*r);
	}
	if(critical.size()==0&&nonCritical.size()==0)
	{
		//nothing will change if continued
		if(__stimulusIndex==__stimulus.size()) timeNext=__result->timeLast();
		//go to next stimulus
		else timeNext=__stimulus[__stimulusIndex]->timeBegin();
	}
	else
	{
		double tau1=-1;
		if(nonCritical.size()>0)
		{
			for(list<Molecule*>::iterator m=__molecule.begin(); m!=__molecule.end(); ++m) (*m)->initMuSigmaError();
			for(list<Reaction*>::iterator r=nonCritical.begin(); r!=nonCritical.end(); ++r)
			{
				(*r)->updateError(__errorCoef);
				(*r)->addMuSigma();
			}
			for(list<Molecule*>::iterator m=__molecule.begin(); m!=__molecule.end(); ++m)
			{
				if((*m)->reactantOfNonCritical())
				{
					double t=(*m)->tau1(__errorCoef);
					if(tau1<0||t<tau1) tau1=t;
				}
			}
			if(tau1!=-1&&tau1<__tau1Min)
			{
				__abandonTauLeap=true;
				return;
				tau1=-1;
			}
		}
		double tau2=-1, flux2=0;
		if(critical.size()>0)
		{
			for(list<Reaction*>::iterator r=critical.begin(); r!=critical.end(); ++r) flux2+=(*r)->flux();
			tau2=-1.0/flux2*log(__random->getDouble());
		}
		if(tau2<0||tau1<tau2)
		{
			for(list<Reaction*>::iterator r=nonCritical.begin(); r!=nonCritical.end(); ++r) (*r)->execute(tau1, __random);
			timeNext=__time+tau1;
		}
		else
		{
			critical.sort(Reaction_compare);
			Reaction* reactionNext;
			double rand=__random->getDouble()*flux2;
			double cumFlux=0;
			for(list<Reaction*>::iterator r=critical.begin(); r!=critical.end(); ++r)
			{
				if((*r)->flux()==0) continue;
				cumFlux+=(*r)->flux();
				if(cumFlux>=rand)
				{
					reactionNext=*r;
					break;
				}
			}
			if(reactionNext==NULL) reactionNext=*(critical.rbegin());
			reactionNext->execute();
		
			for(list<Reaction*>::iterator r=nonCritical.begin(); r!=nonCritical.end(); ++r) (*r)->execute(tau2, __random);
			timeNext=__time+tau2;
		}
		for(list<Molecule*>::iterator m=__molecule.begin(); m!=__molecule.end(); ++m)
		{
			if((*m)->num()<0)
			{
				__abandonTauLeap=true;
				__restoreNumPrev();
				return;
			}
		}
	}

	__stepStim(timeNext);
	__stepResult(timeNext);

	__time=timeNext;
}

void Simulation::stepRungeKutta()
{
	__copyToNumPrev();

	for(list<Molecule*>::iterator m=__molecule.begin(); m!=__molecule.end(); ++m) (*m)->initK();
	for(int k=0; k<NUM_K; ++k)
	{
		for(list<Reaction*>::iterator r=__reaction.begin(); r!=__reaction.end(); ++r) (*r)->flux(__K_WEIGHT_MATRIX[k]);
		for(list<Reaction*>::iterator r=__reaction.begin(); r!=__reaction.end(); ++r) (*r)->setK(k);
		for(list<Molecule*>::iterator m=__molecule.begin(); m!=__molecule.end(); ++m) (*m)->multiplyDtToK(k, __dTime);
	}
	for(list<Molecule*>::iterator m=__molecule.begin(); m!=__molecule.end(); ++m) (*m)->updateNum(__K_WEIGHT);

	for(list<Molecule*>::iterator m=__molecule.begin(); m!=__molecule.end(); ++m)
	{
		if((*m)->num()<0)
		{
			cerr<<"Negative number: "<<(*m)->name()<<endl;
			exit(1);
		}
	}

	double timeNext=__time+__dTime;
	if(__stepStim(timeNext)) __dTime=__dTimeMin;
	else
	{
		double error=0;
		for(list<Molecule*>::iterator m=__molecule.begin(); m!=__molecule.end(); ++m)
		{
			double e=(*m)->kSum(__K_WEIGHT_FOR_ERROR);
			error+=e*e;
		}
		double dtPrev=__dTime;
		__dTime=__errorCoef*pow(__dTime, 5.0/4)/pow(error, 1.0/8);
		if(__dTime<__dTimeMin) __dTime=__dTimeMin;
		else if(__dTime>__dTimeMax) __dTime=__dTimeMax;
	}

	__stepResult(timeNext);
	__time=timeNext;
}

void Simulation::__copyToNumPrev()
{
	for(list<Molecule*>::iterator m=__molecule.begin(); m!=__molecule.end(); ++m) (*m)->copyToNumPrev();
}

bool Simulation::__stepStim(double time)
{
	bool updated=false;
	for(int s=__stimulusIndex; s<__stimulus.size(); ++s)
	{
		Stimulus* stim=__stimulus[s];
		if(stim->timeBegin()>=time) break;
		if(stim->step(time)) updated=true;
	}
	while(__stimulusIndex<__stimulus.size()&&__stimulus[__stimulusIndex]->finished()) ++__stimulusIndex;
	return updated;
}

void Simulation::__stepResult(double timeNext)
{
	int riFirst, num=1;
	if(timeNext>=__timeResult[__timeResultIndex])
	{
		riFirst=__timeResultIndex;
		while(riFirst+num<__timeResult.size()&&timeNext>=__timeResult[riFirst+num]) ++num;
		if(num==1)
		{
			if(timeNext-__timeResult[__timeResultIndex]<=__timeResult[__timeResultIndex]-__time)
				__result->add(timeNext, false);
			else if(__result->size()==0||__result->timeLast()<__time) 
				__result->add(__time, true);
		}
		else 
		{
			if(__result->size()==0||__result->timeLast()<__time)
				__result->add(__time, true);
			__result->add(timeNext, false);
		}
		__resultUpdated=true;
		__timeResultIndex+=num;
	}
	else __resultUpdated=false;
}

bool Simulation::finished() {return __timeResultIndex==__timeResult.size();}
bool Simulation::resultUpdated() {return __resultUpdated;}

Result* Simulation::result() {return __result;}
double Simulation::time() {return __time;}

bool Simulation::abandonTauLeap()
{
	return __abandonTauLeap;
}

void Simulation::__restoreNumPrev()
{
	for(list<Molecule*>::iterator m=__molecule.begin(); m!=__molecule.end(); ++m) (*m)->restoreNumPrev();
}
