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

#include "String.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

String::String()
{}
String::String(string value) : string(value)
{}
String::String(const char *value) : string(value)
{}
String::String(int size, const char *format, ...)
{
	char *tmp = new char[size];
	va_list list;
    va_start(list, format);
    vsprintf(tmp, format, list);
    va_end(list);
	*this=tmp;
	delete []tmp;
}

vector<String> String::split(char deliminator)
{
	vector<String> token;
	int tokenBegin=0;
	for(int pos=0; pos<size(); pos++)
	{
		if((*this)[pos]==deliminator)
		{
			token.push_back(String(substr(tokenBegin, pos-tokenBegin)));
			tokenBegin=pos+1;
		}
	}
	token.push_back(String(substr(tokenBegin)));
	return token;
}

String String::sub(int begin, int size)
{
	return substr(begin, size);
}

int String::toInt()
{
	return atoi(c_str());
}
double String::toDouble()
{
	return atof(c_str());
}
