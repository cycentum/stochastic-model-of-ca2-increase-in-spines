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

#ifndef CONTAINER2_H
#define CONTAINER2_H

template <class T0, class T1> class Container2
{
public:
	T0 value0;
	T1 value1;

	Container2(T0 v0, T1 v1) : value0(v0), value1(v1){}
	Container2();
};

#endif
