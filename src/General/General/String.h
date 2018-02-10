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

#ifndef STRING_H
#define STRING_H

#define _CRT_SECURE_NO_WARNINGS

#include <string>
#include <vector>
using namespace std;

class String : public string
{
public:
	String();
	String(string value);
	String(const char *value);
	String(int size, const char *format, ...);
	vector<String> split(char deliminator);
	String sub(int begin, int size);
	int toInt();
	double toDouble();
};

#endif
