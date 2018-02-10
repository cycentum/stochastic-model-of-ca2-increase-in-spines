#include "Random.h"
#include <math.h>
#include <stdlib.h>

#define M_PI 3.14159265358979323846264338328

double stirling(double);
double pow(double x, unsigned int n);
#define SMALL_MEAN 14           /* If n*p < SMALL_MEAN then use BINV
                                   algorithm. The ranlib
                                   implementation used cutoff=30; but
                                   on my computer 14 works better */

#define BINV_CUTOFF 110         /* In BINV, do not permit ix too large */

#define FAR_FROM_MEAN 20        /* If ix-n*p is larger than this, then
                                   use the "squeeze" algorithm.
                                   Ranlib used 20, and this seems to
                                   be the best choice on my machine as
                                   well */

Random::Random(int seed) : __seed(seed), __finishedNum(0)
{
	dsfmt_init_gen_rand(&__dsfmt, __seed);
}
Random::Random(dsfmt_t d) : __dsfmt(d){}

double Random::getDouble()
{
	++__finishedNum;
	return dsfmt_genrand_close_open(&__dsfmt);
}

double Random::getDoubleOpenOpen()
{
	++__finishedNum;
	return dsfmt_genrand_open_open(&__dsfmt);
}

int Random::getPoisson(double mean)
{
	double emu;
	double prod = 1.0;
	unsigned int k = 0;

	while (mean > 10)
	{
		unsigned int m = mean * (7.0 / 8.0);

		double X = getGamma(m);

		if (X >= mean)
		{
			return k + getBinomial(mean / X, m - 1);
		}
		else
		{
			k += m;
			mean -= X; 
		}
	}

	/* This following method works well when mu is small */

	emu = exp (-mean);

	do
	{
		prod *= getDouble();
		k++;
	}
	while (prod > emu);

	return k - 1;
}

double Random::getGamma(int a)
{
	if (a < 12)
	{
		unsigned int i;
		double prod = 1;

		for (i = 0; i < a; i++)
		{
			prod *= getDoubleOpenOpen();
		}

		/* Note: for 12 iterations we are safe against underflow, since
		the smallest positive random number is O(2^-32). This means
		the smallest possible product is 2^(-12*32) = 10^-116 which
		is within the range of double precision. */

		return -log (prod);
	}
	else
	{
		return __getGammaLargeA((double) a);
	}
}

double Random::__getGammaLargeA(double a)
{
	/* Works only if a > 1, and is most efficient if a is large

	This algorithm, reported in Knuth, is attributed to Ahrens.  A
	faster one, we are told, can be found in: J. H. Ahrens and
	U. Dieter, Computing 12 (1974) 223-246.  */

	double sqa, x, y, v;
	sqa = sqrt (2 * a - 1);
	do
	{
		do
		{
			y = tan (M_PI * getDouble());
			x = sqa * y + a - 1;
		}
		while (x <= 0);
		v = getDouble();
	}
	while (v > (1 + y * y) * exp ((a - 1) * log (x / (a - 1)) - sqa * y));

	return x;
}

double Random::getBinomial(double p, double n)
{
	int ix;                       /* return value */
	int flipped = 0;
	double q, s, np;

	if (n == 0)
		return 0;

	if (p > 0.5)
	{
		p = 1.0 - p;              /* work with small p */
		flipped = 1;
	}

	q = 1 - p;
	s = p / q;
	np = n * p;

	/* Inverse cdf logic for small mean (BINV in K+S) */

	if (np < SMALL_MEAN)
	{
		double f0 = pow (q, n);   /* f(x), starting with x=0 */

		while (1)
		{
			/* This while(1) loop will almost certainly only loop once; but
			* if u=1 to within a few epsilons of machine precision, then it
			* is possible for roundoff to prevent the main loop over ix to
			* achieve its proper value.  following the ranlib implementation,
			* we introduce a check for that situation, and when it occurs,
			* we just try again.
			*/

			double f = f0;
			double u = getDouble();

			for (ix = 0; ix <= BINV_CUTOFF; ++ix)
			{
				if (u < f)
					goto Finish;
				u -= f;
				/* Use recursion f(x+1) = f(x)*[(n-x)/(x+1)]*[p/(1-p)] */
				f *= s * (n - ix) / (ix + 1);
			}

			/* It should be the case that the 'goto Finish' was encountered
			* before this point was ever reached.  But if we have reached
			* this point, then roundoff has prevented u from decreasing
			* all the way to zero.  This can happen only if the initial u
			* was very nearly equal to 1, which is a rare situation.  In
			* that rare situation, we just try again.
			*
			* Note, following the ranlib implementation, we loop ix only to
			* a hardcoded value of SMALL_MEAN_LARGE_N=110; we could have
			* looped to n, and 99.99...% of the time it won't matter.  This
			* choice, I think is a little more robust against the rare
			* roundoff error.  If n>LARGE_N, then it is technically
			* possible for ix>LARGE_N, but it is astronomically rare, and
			* if ix is that large, it is more likely due to roundoff than
			* probability, so better to nip it at LARGE_N than to take a
			* chance that roundoff will somehow conspire to produce an even
			* larger (and more improbable) ix.  If n<LARGE_N, then once
			* ix=n, f=0, and the loop will continue until ix=LARGE_N.
			*/
		}
	}
	else
	{
		/* For n >= SMALL_MEAN, we invoke the BTPE algorithm */

		int k;

		double ffm = np + p;      /* ffm = n*p+p             */
		int m = (int) ffm;        /* m = int floor[n*p+p]    */
		double fm = m;            /* fm = double m;          */
		double xm = fm + 0.5;     /* xm = half integer mean (tip of triangle)  */
		double npq = np * q;      /* npq = n*p*q            */

		/* Compute cumulative area of tri, para, exp tails */

		/* p1: radius of triangle region; since height=1, also: area of region */
		/* p2: p1 + area of parallelogram region */
		/* p3: p2 + area of left tail */
		/* p4: p3 + area of right tail */
		/* pi/p4: probability of i'th area (i=1,2,3,4) */

		/* Note: magic numbers 2.195, 4.6, 0.134, 20.5, 15.3 */
		/* These magic numbers are not adjustable...at least not easily! */

		double p1 = floor (2.195 * sqrt (npq) - 4.6 * q) + 0.5;

		/* xl, xr: left and right edges of triangle */
		double xl = xm - p1;
		double xr = xm + p1;

		/* Parameter of exponential tails */
		/* Left tail:  t(x) = c*exp(-lambda_l*[xl - (x+0.5)]) */
		/* Right tail: t(x) = c*exp(-lambda_r*[(x+0.5) - xr]) */

		double c = 0.134 + 20.5 / (15.3 + fm);
		double p2 = p1 * (1.0 + c + c);

		double al = (ffm - xl) / (ffm - xl * p);
		double lambda_l = al * (1.0 + 0.5 * al);
		double ar = (xr - ffm) / (xr * q);
		double lambda_r = ar * (1.0 + 0.5 * ar);
		double p3 = p2 + c / lambda_l;
		double p4 = p3 + c / lambda_r;

		double var, accept;
		double u, v;              /* random variates */

TryAgain:

		/* generate random variates, u specifies which region: Tri, Par, Tail */
		u = getDouble() * p4;
		v = getDouble();

		if (u <= p1)
		{
			/* Triangular region */
			ix = (int) (xm - p1 * v + u);
			goto Finish;
		}
		else if (u <= p2)
		{
			/* Parallelogram region */
			double x = xl + (u - p1) / c;
			v = v * c + 1.0 - fabs (x - xm) / p1;
			if (v > 1.0 || v <= 0.0)
				goto TryAgain;
			ix = (int) x;
		}
		else if (u <= p3)
		{
			/* Left tail */
			ix = (int) (xl + log (v) / lambda_l);
			if (ix < 0)
				goto TryAgain;
			v *= ((u - p2) * lambda_l);
		}
		else
		{
			/* Right tail */
			ix = (int) (xr - log (v) / lambda_r);
			if (ix > (double) n)
				goto TryAgain;
			v *= ((u - p3) * lambda_r);
		}

		/* At this point, the goal is to test whether v <= f(x)/f(m) 
		*
		*  v <= f(x)/f(m) = (m!(n-m)! / (x!(n-x)!)) * (p/q)^{x-m}
		*
		*/

		/* Here is a direct test using logarithms.  It is a little
		* slower than the various "squeezing" computations below, but
		* if things are working, it should give exactly the same answer
		* (given the same random number seed).  */

#ifdef DIRECT
		var = log (v);

		accept =
			LNFACT (m) + LNFACT (n - m) - LNFACT (ix) - LNFACT (n - ix)
			+ (ix - m) * log (p / q);

#else /* SQUEEZE METHOD */

		/* More efficient determination of whether v < f(x)/f(M) */

		k = abs (ix - m);

		if (k <= FAR_FROM_MEAN)
		{
			/* 
			* If ix near m (ie, |ix-m|<FAR_FROM_MEAN), then do
			* explicit evaluation using recursion relation for f(x)
			*/
			double g = (n + 1) * s;
			double f = 1.0;

			var = v;

			if (m < ix)
			{
				int i;
				for (i = m + 1; i <= ix; i++)
				{
					f *= (g / i - s);
				}
			}
			else if (m > ix)
			{
				int i;
				for (i = ix + 1; i <= m; i++)
				{
					f /= (g / i - s);
				}
			}

			accept = f;
		}
		else
		{
			/* If ix is far from the mean m: k=ABS(ix-m) large */

			var = log (v);

			if (k < npq / 2 - 1)
			{
				/* "Squeeze" using upper and lower bounds on
				* log(f(x)) The squeeze condition was derived
				* under the condition k < npq/2-1 */
				double amaxp =
					k / npq * ((k * (k / 3.0 + 0.625) + (1.0 / 6.0)) / npq + 0.5);
				double ynorm = -(k * k / (2.0 * npq));
				if (var < ynorm - amaxp)
					goto Finish;
				if (var > ynorm + amaxp)
					goto TryAgain;
			}

			/* Now, again: do the test log(v) vs. log f(x)/f(M) */

#if USE_EXACT
			/* This is equivalent to the above, but is a little (~20%) slower */
			/* There are five log's vs three above, maybe that's it? */

			accept = LNFACT (m) + LNFACT (n - m)
				- LNFACT (ix) - LNFACT (n - ix) + (ix - m) * log (p / q);

#else /* USE STIRLING */
			/* The "#define Stirling" above corresponds to the first five
			* terms in asymptoic formula for
			* log Gamma (y) - (y-0.5)log(y) + y - 0.5 log(2*pi);
			* See Abramowitz and Stegun, eq 6.1.40
			*/

			/* Note below: two Stirling's are added, and two are
			* subtracted.  In both K+S, and in the ranlib
			* implementation, all four are added.  I (jt) believe that
			* is a mistake -- this has been confirmed by personal
			* correspondence w/ Dr. Kachitvichyanukul.  Note, however,
			* the corrections are so small, that I couldn't find an
			* example where it made a difference that could be
			* observed, let alone tested.  In fact, define'ing Stirling
			* to be zero gave identical results!!  In practice, alv is
			* O(1), ranging 0 to -10 or so, while the Stirling
			* correction is typically O(10^{-5}) ...setting the
			* correction to zero gives about a 2% performance boost;
			* might as well keep it just to be pendantic.  */

			{
				double x1 = ix + 1.0;
				double w1 = n - ix + 1.0;
				double f1 = fm + 1.0;
				double z1 = n + 1.0 - fm;

				accept = xm * log (f1 / x1) + (n - m + 0.5) * log (z1 / w1)
					+ (ix - m) * log (w1 * p / (x1 * q))
					+ stirling (f1) + stirling (z1) - stirling (x1) - stirling (w1);
			}
#endif
#endif
		}


		if (var <= accept)
		{
			goto Finish;
		}
		else
		{
			goto TryAgain;
		}
	}

Finish:

	return (flipped) ? (n - ix) : (unsigned int)ix;
}

unsigned long Random::finishedNum(){return __finishedNum;}

int Random::seed(){return __seed;}
/*
void Random::save(BinaryWriter& writer)
{
	for(int i=0; i<DSFMT_N+1; ++i)
		for(int j=0; j<2; ++j)
			writer.write(__dsfmt.status[i].d[j]);
	writer.write(__dsfmt.idx);
}

Random Random::load(BinaryReader& reader)
{
	dsfmt_t d;
	for(int i=0; i<DSFMT_N+1; ++i)
		for(int j=0; j<2; ++j)
			d.status[i].d[j]=reader.readDouble();
	d.idx=reader.readInt();
	return Random(d);
}*/

double stirling (double y1)
{
	double y2 = y1 * y1;
	double s =
		(13860.0 -
		(462.0 - (132.0 - (99.0 - 140.0 / y2) / y2) / y2) / y2) / y1 / 166320.0;
	return s;
}

double gsl_pow_int(double x, int n)
{
	double value = 1.0;

	/* repeated squaring method 
	* returns 0.0^0 = 1.0, so continuous in x
	*/
	do {
		if(n & 1) value *= x;  /* for n odd */
		n >>= 1;
		x *= x;
	} while (n);

	return value;
}
