#ifndef RANDOM_H
#define RANDOM_H

#include "dSFMT.h"
//#include "../IO/BinaryWriter.h"
//#include "../IO/BinaryReader.h"

/**
data structure of binary
{
 double (status[2])[DSFMT_N+1]
 int index;
}
*/
class Random
{
public:
	dsfmt_t __dsfmt;
	unsigned long __finishedNum;
	int __seed;
	Random(dsfmt_t dsfmt);
	double __getGammaLargeA(double a);

public:
	Random(int seed);
	double getDouble();
	double getDoubleOpenOpen();
	int getPoisson(double mean);
	double getGamma(int a);
	double getBinomial(double p, double n);
	unsigned long finishedNum();
	int seed();
//	void save(BinaryWriter& writer);
//	static Random load(BinaryReader& reader);
};

#endif
