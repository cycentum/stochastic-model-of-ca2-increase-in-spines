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

#include "Result.h"
#include "Molecule.h"
#include "UnitConversion.h"
#include <iostream>

Result::Result(list<Molecule*> mole) : __molecule(mole){}

/*
void Result::add(double timePrev, double timeNext, double timeResult, list<Molecule*> molecule)
{
	double dt=timeNext-timePrev;
	for(list<Molecule*>::iterator m=molecule.begin(); m!=molecule.end(); ++m)
	{
		double n=((*m)->numPrev()*timeNext-(*m)->num()*timePrev+((*m)->num()-(*m)->numPrev())*timeResult)/dt;
		__num[*m].push_back(n);
	}
}*/
void Result::add(double time, bool numPrev)
{
	__time.push_back(time);
	list<double> n;
	for(list<Molecule*>::iterator m=__molecule.begin(); m!=__molecule.end(); ++m)
	{
		if(numPrev) n.push_back((*m)->numPrev());
		else n.push_back((*m)->num());
	}
	__num.push_back(n);
}
/*
void Result::saveConc(BinaryWriter& writer)
{
	writer.write((int)__time.size());
	list<double>::iterator t=__time.begin();
	for(list<list<double> >::iterator n=__num.begin(); n!=__num.end(); ++n)
	{
		writer.write(*t);
		++t;
		list<Molecule*>::iterator m=__molecule.begin();
		for(list<double>::iterator v=n->begin(); v!=n->end(); ++v)
		{
			writer.write(UnitConversion::num_to_umol(*v/(*m)->volume()));
			++m;
		}
	}
}

void Result::saveNum(BinaryWriter& writer)
{
	writer.write((int)__time.size());
	list<double>::iterator t=__time.begin();
	for(list<list<double> >::iterator n=__num.begin(); n!=__num.end(); ++n)
	{
		writer.write(*t);
		++t;
		for(list<double>::iterator v=n->begin(); v!=n->end(); ++v)
		{
			writer.write(*v);
		}
	}
}*/

double Result::timeLast()
{
	return *__time.rbegin();
}

int Result::size(){return __time.size();}

void Result::saveConc(ofstream& ofs)
{
	ofs<<"\t";
	for(list<Molecule*>::iterator m=__molecule.begin(); m!=__molecule.end(); ++m)
	{
		ofs<<(*m)->name()<<"\t";
	}
	ofs<<endl;
	list<double>::iterator t=__time.begin();
	for(list<list<double> >::iterator n=__num.begin(); n!=__num.end(); ++n)
	{
		ofs<<(*t)<<"\t";
		++t;
		list<Molecule*>::iterator m=__molecule.begin();
		for(list<double>::iterator v=n->begin(); v!=n->end(); ++v)
		{
			ofs<<UnitConversion::num_to_umol(*v/(*m)->volume())<<"\t";
			++m;
		}
		ofs<<endl;
	}
}

void Result::saveNum(ofstream& ofs)
{
	ofs<<"\t";
	for(list<Molecule*>::iterator m=__molecule.begin(); m!=__molecule.end(); ++m)
	{
		ofs<<(*m)->name()<<"\t";
	}
	ofs<<endl;
	list<double>::iterator t=__time.begin();
	for(list<list<double> >::iterator n=__num.begin(); n!=__num.end(); ++n)
	{
		ofs<<(*t)<<"\t";
		++t;
		list<Molecule*>::iterator m=__molecule.begin();
		for(list<double>::iterator v=n->begin(); v!=n->end(); ++v)
		{
			ofs<<(*v)<<"\t";
			++m;
		}
		ofs<<endl;
	}
}
