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
#include "Simulation.h"
#include "Molecule.h"
#include "Reaction.h"
#include "Stimulus.h"
#include "UnitConversion.h"
#include "Param.h"
//#include "../General/IO/BinaryWriter.h"
#include <iostream>
#include <ctime>

Model::Model()
{
	__loadMole();
	__loadReac();
}

Model::~Model()
{
	for(list<Molecule*>::iterator m=__moleList.begin(); m!=__moleList.end(); ++m) delete *m;
	for(map<String, Reaction*>::iterator r=__reac.begin(); r!=__reac.end(); ++r) delete r->second;
	for(vector<Stimulus*>::iterator s=__stim.begin(); s!=__stim.end(); ++s) delete *s;
	delete __simulation;
}

void Model::__loadMole()
{
	list<Mole> model=getMole();
	for(list<Mole>::iterator m=model.begin(); m!=model.end(); ++m)
	{
		Molecule *mo=new Molecule(m->name, UnitConversion::um3_to_l(m->volume), m->num);
		if(mo->name().size()>=4 && mo->name().compare(mo->name().size()-3, 3, "Fix")==0) mo->fix(true);
		__mole.insert(pair<String, Molecule*>(m->name, mo));
		__moleList.push_back(mo);
	}
}

void Model::__loadReac()
{
	list<Reac> model=getReac();
	for(list<Reac>::iterator r=model.begin(); r!=model.end(); ++r)
	{
		String nameF=r->name;
		double rateConstF=r->rate;
		map<Molecule*, int> moleF;
		for(list<Container2<String, int> >::iterator n=r->num.begin(); n!=r->num.end(); ++n)
			moleF.insert(pair<Molecule*, int>(__mole[n->value0], n->value1));
		
		++r;
		String nameB=r->name;
		double rateConstB=r->rate;
		map<Molecule*, int> moleB;
		for(list<Container2<String, int> >::iterator n=r->num.begin(); n!=r->num.end(); ++n)
			moleB.insert(pair<Molecule*, int>(__mole[n->value0], n->value1));

		{
			Reaction* re=Reaction_new_from_conc_uM_volume_um3(nameF, rateConstF, moleF, moleB);
			__reac.insert(pair<String, Reaction*>(nameF, re));
		}
		{
			Reaction* re=Reaction_new_from_conc_uM_volume_um3(nameB, rateConstB, moleB, moleF);
			__reac.insert(pair<String, Reaction*>(nameB, re));
		}
	}
}

void Model::setStim(double timing, double ampPFCa, double ampPFGlu, int numPF)
{
	map<String, Stim> model=getStim();

	for(int n=0; n<numPF; ++n)
	{
		if(ampPFCa>0)
		{
			Stim s=model["PFCa"];
			__stim.push_back(new Stimulus(__mole[s.mole], s.time+n*0.01, s.duration, s.num*ampPFCa));
		}
		if(ampPFGlu>0)
		{
			Stim s=model["PFGlu"];
			__stim.push_back(new Stimulus(__mole[s.mole], s.time+n*0.01, s.duration, s.num*ampPFGlu));
		}
	}
	if(timing<1.5)
	{
		Stim s=model["CF"];
		__stim.push_back(new Stimulus(__mole[s.mole], s.time+timing, s.duration, s.num));
	}
}
void Model::setStim(double timing, vector<double> amp)
{
	map<String, Stim> model=getStim();

	for(int a=0; a<amp.size(); ++a)
	{
		if(amp[a]<=0) continue;
		{
			Stim s=model["PFCa"];
			__stim.push_back(new Stimulus(__mole[s.mole], s.time+a*0.01, s.duration, s.num*amp[a]));
		}
		{
			Stim s=model["PFGlu"];
			__stim.push_back(new Stimulus(__mole[s.mole], s.time+a*0.01, s.duration, s.num*amp[a]));
		}
	}
	if(timing<1.5)
	{
		Stim s=model["CF"];
		__stim.push_back(new Stimulus(__mole[s.mole], s.time+timing, s.duration, s.num));
	}
}

void Model::setStimGlu(double numPerMS, double durMS)
{
	if(numPerMS>0&&durMS>0)
	{
		map<String, Stim> model=getStim();
		Stim s=model["PFGlu"];
		__stim.push_back(new Stimulus(__mole[s.mole], s.time, durMS/1000.0, numPerMS*durMS));
	}
}

void Model::setStimTest()
{
	map<String, Stim> model=getStim();
	Stim s=model["A"];
	__stim.push_back(new Stimulus(__mole[s.mole], s.time, s.duration, s.num));
}

void Model::setStocahsticity(double stc)
{
	for(map<String, Reaction*>::iterator r=__reac.begin(); r!=__reac.end(); ++r) r->second->stochasticity(stc);
	for(vector<Stimulus*>::iterator s=__stim.begin(); s!=__stim.end(); ++s) (*s)->stochasticity(stc);
}

void Model::tauLeap(bool tl){__tauLeap=tl;}

void Model::run(Random* random)
{
	list<Reaction*> reac;
	for(map<String, Reaction*>::iterator r=__reac.begin(); r!=__reac.end(); ++r) reac.push_back(r->second);
	vector<double> timeResult;
	for(int t=0; t<2001; ++t) timeResult.push_back(t*0.001);
	if(random==NULL) __simulation=new Simulation(__moleList, reac, __stim, timeResult, 1e-7, 1e-5, 100);
	else if(!__tauLeap) __simulation=new Simulation(__moleList, reac, __stim, timeResult, random);	
	else __simulation=new Simulation(__moleList, reac, __stim, timeResult, random, 0.03, 1e-9, 10);

//	clock_t clockBegin=clock();
	int tauLeapAbandoned=0;
	while(!__simulation->finished())
	{
		if(random==NULL) __simulation->stepRungeKutta();
		else if(!__tauLeap||tauLeapAbandoned>0)
		{
			__simulation->stepGillespie();
			if(tauLeapAbandoned>0) --tauLeapAbandoned;
		}
		else
		{
			__simulation->stepTauLeap();
			if(__simulation->abandonTauLeap()) tauLeapAbandoned=100;
		}
		if(__simulation->resultUpdated()) cout<<"Model::run: time="<<__simulation->time()<<endl;
	}
//	clock_t clockEnd=clock();
//	__clockS=(double)(clockEnd-clockBegin)/CLOCKS_PER_SEC;
}

Result* Model::result(){return __simulation->result();}
/*
void Model::saveConc(int trial, String path)
{
	BinaryWriter writer(path, 1024*16);
	writer.write(1);
	writer.write(trial);
	writer.write(trial);
	result()->saveConc(writer);
	writer.close();
}
*/
