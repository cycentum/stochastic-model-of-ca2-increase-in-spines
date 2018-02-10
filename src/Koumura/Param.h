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

#ifndef PARAM_H
#define PARAM_H

#include "../General/General/String.h"
#include "../General/Collection/Container2.h"
#include <list>
#include <map>

class Mole
{
public:
	String name;
	int num;
	double volume;
	Mole(String n, int nu, double vol) : name(n), num(nu), volume(vol){}
	Mole(){}
};
list<Mole> getMole();

class Reac
{
public:
	String name;
	double rate;
	list<Container2<String, int> > num;
	Reac(String n, double r) : name(n), rate(r){}
	Reac(){}
	void add(String n, int nu)
	{
		num.push_back(Container2<String, int>(n, nu));
	}
};
list<Reac> getReac();

class Stim
{
public:
	String name;
	String mole;
	double time, duration;
	int num;
	Stim(String na, String m, double t, double d, int nu) : name(na), mole(m), time(t), duration(d), num(nu){}
	Stim(){}
};
map<String, Stim> getStim();

#endif
