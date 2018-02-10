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

#include "Param.h"

//#define DEBUG
//#define FIX_d1
//#define FIX_d2
//#define CUT_e1

#ifndef DEBUG
list<Mole> getMole()
{
	list<Mole> mole;
	mole.push_back(Mole("Glu", 0, 0.002)); 
	mole.push_back(Mole("mGluR", 10, 0.002)); 
	mole.push_back(Mole("Glu-mGluR", 0, 0.002)); 
	mole.push_back(Mole("Gq-GDP", 52, 0.002)); 
	mole.push_back(Mole("mGluR-Gq", 8, 0.002)); 
	mole.push_back(Mole("Glu-mGluR-Gq", 0, 0.002)); 
	mole.push_back(Mole("Ga-GTP", 0, 0.002)); 
	mole.push_back(Mole("Gbc", 0, 0.002)); 
	mole.push_back(Mole("Ga-GDP", 0, 0.002)); 
	mole.push_back(Mole("PIP2", 5000, 0.002)); 
	mole.push_back(Mole("PLC-PIP2", 42, 0.002)); 
	mole.push_back(Mole("PLC-Ca-PIP2", 8, 0.002)); 
	mole.push_back(Mole("PLC-Gq-PIP2", 0, 0.002)); 
	mole.push_back(Mole("PLC-Ca-Gq-PIP2", 0, 0.002)); 
	mole.push_back(Mole("PLC-Ca", 1, 0.002)); 
	mole.push_back(Mole("PLC-Ca-Gq", 0, 0.002)); 
	mole.push_back(Mole("DAG", 0, 0.002)); 
	mole.push_back(Mole("IP3_PSD", 0, 0.002)); 
	mole.push_back(Mole("IP3", 6, 0.1)); 
	mole.push_back(Mole("IP3_3-kinase", 52, 0.1)); 
	mole.push_back(Mole("IP3K-2Ca", 2, 0.1)); 
	mole.push_back(Mole("IP3K-2Ca-IP3", 0, 0.1)); 
	mole.push_back(Mole("IP3_5-phosphatase", 59, 0.1)); 
	mole.push_back(Mole("IP5P-IP3", 2, 0.1)); 
	mole.push_back(Mole("IP3Rec", 14, 0.1)); 
	mole.push_back(Mole("IP3R-IP3", 0, 0.1)); 
	mole.push_back(Mole("IP3R-open", 0, 0.1)); 
	mole.push_back(Mole("IP3R-Ca", 1, 0.1)); 
	mole.push_back(Mole("IP3R-2Ca", 0, 0.1)); 
	mole.push_back(Mole("IP3R-3Ca", 0, 0.1)); 
	mole.push_back(Mole("IP3R-4Ca", 0, 0.1)); 
	mole.push_back(Mole("Ca_PSD", 0, 0.002)); 
	mole.push_back(Mole("Ca", 4, 0.1)); 
	mole.push_back(Mole("SERCA", 148, 0.1)); 
	mole.push_back(Mole("SERCA-2Ca", 7, 0.1)); 
	mole.push_back(Mole("PMCA", 68, 0.1)); 
	mole.push_back(Mole("PMCA-Ca", 40, 0.1)); 
	mole.push_back(Mole("NCX", 32, 0.1)); 
	mole.push_back(Mole("NCX-2Ca", 0, 0.1)); 
	mole.push_back(Mole("CaStore", 1800, 0.02)); 
	mole.push_back(Mole("Calreticulin", 960000, 0.02)); 
	mole.push_back(Mole("Calreticulin-Ca", 72000, 0.02)); 
	mole.push_back(Mole("Ca_ext", 12000000, 10.0)); 
	mole.push_back(Mole("MgGreen", 14940, 0.1)); 
	mole.push_back(Mole("MgGreenAst", 60, 0.1)); 
	mole.push_back(Mole("parvalbumin", 1380, 0.1)); 
	mole.push_back(Mole("PV-Ca", 1620, 0.1)); 
	mole.push_back(Mole("Calbindin-D28k", 5850, 0.1)); 
	mole.push_back(Mole("CB-Ca", 150, 0.1)); 
	mole.push_back(Mole("LowAffBuf", 5997, 0.1)); 
	mole.push_back(Mole("LAB-Ca", 3, 0.1)); 
	mole.push_back(Mole("LowAffBuf2", 6000, 0.1)); 
	mole.push_back(Mole("LAB2-Ca", 0, 0.1)); 
	mole.push_back(Mole("IP4", 0, 0.1)); 
	mole.push_back(Mole("IP2", 0, 0.1)); 
	mole.push_back(Mole("Glu-Decay", 0, 0.002)); 

#ifdef FIX_d1
	mole.push_back(Mole("IP3_Fix", 6, 0.1)); 
#endif
#ifdef FIX_d2
	mole.push_back(Mole("Ca_Fix", 4, 0.1)); 
#endif
	return mole;
}

list<Reac> getReac()
{
	list<Reac> reac;
	{ 
		Reac r("a1f", 11.1); 
		r.add("mGluR", 1); 
		r.add("Glu", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("a1b", 100.0); 
		r.add("Glu-mGluR", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("a2f", 11.1); 
		r.add("mGluR-Gq", 1); 
		r.add("Glu", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("a2b", 100.0); 
		r.add("Glu-mGluR-Gq", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("a3f", 2.0); 
		r.add("mGluR", 1); 
		r.add("Gq-GDP", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("a3b", 100.0); 
		r.add("mGluR-Gq", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("a4f", 2.0); 
		r.add("Glu-mGluR", 1); 
		r.add("Gq-GDP", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("a4b", 100.0); 
		r.add("Glu-mGluR-Gq", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("a5f", 116.0); 
		r.add("Glu-mGluR-Gq", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("a5b", 0.0); 
		r.add("Glu-mGluR", 1); 
		r.add("Gbc", 1); 
		r.add("Ga-GTP", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("a6f", 1.0E-4); 
		r.add("Gq-GDP", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("a6b", 0.0); 
		r.add("Gbc", 1); 
		r.add("Ga-GTP", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("a7f", 0.02); 
		r.add("Ga-GTP", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("a7b", 0.0); 
		r.add("Ga-GDP", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("a8f", 6.0); 
		r.add("Gbc", 1); 
		r.add("Ga-GDP", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("a8b", 0.0); 
		r.add("Gq-GDP", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("b1f", 300.0); 
		r.add("PLC-PIP2", 1); 
		r.add("Ca_PSD", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("b1b", 100.0); 
		r.add("PLC-Ca-PIP2", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("b2f", 900.0); 
		r.add("PLC-Gq-PIP2", 1); 
		r.add("Ca_PSD", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("b2b", 30.0); 
		r.add("PLC-Ca-Gq-PIP2", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("b3f", 800.0); 
		r.add("PLC-PIP2", 1); 
		r.add("Ga-GTP", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("b3b", 40.0); 
		r.add("PLC-Gq-PIP2", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("b4f", 1200.0); 
		r.add("PLC-Ca-PIP2", 1); 
		r.add("Ga-GTP", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("b4b", 6.0); 
		r.add("PLC-Ca-Gq-PIP2", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("b5f", 1200.0); 
		r.add("PLC-Ca", 1); 
		r.add("Ga-GTP", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("b5b", 6.0); 
		r.add("PLC-Ca-Gq", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("b6f", 2.0); 
		r.add("PLC-Ca-PIP2", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("b6b", 0.0); 
		r.add("PLC-Ca", 1); 
		r.add("DAG", 1); 
		r.add("IP3_PSD", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("b7f", 160.0); 
		r.add("PLC-Ca-Gq-PIP2", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("b7b", 0.0); 
		r.add("PLC-Ca-Gq", 1); 
		r.add("DAG", 1); 
		r.add("IP3_PSD", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("b8f", 1.0); 
		r.add("PLC-Ca", 1); 
		r.add("PIP2", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("b8b", 170.0); 
		r.add("PLC-Ca-PIP2", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("b9f", 1.0); 
		r.add("PLC-Ca-Gq", 1); 
		r.add("PIP2", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("b9b", 170.0); 
		r.add("PLC-Ca-Gq-PIP2", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("b10f", 8.0); 
		r.add("PLC-Gq-PIP2", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("b10b", 0.0); 
		r.add("PLC-PIP2", 1); 
		r.add("Ga-GDP", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("b11f", 8.0); 
		r.add("PLC-Ca-Gq-PIP2", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("b11b", 0.0); 
		r.add("PLC-Ca-PIP2", 1); 
		r.add("Ga-GDP", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("b12f", 8.0); 
		r.add("PLC-Ca-Gq", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("b12b", 0.0); 
		r.add("PLC-Ca", 1); 
		r.add("Ga-GDP", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("c1f", 1111.1); 
		r.add("IP3_3-kinase", 1); 
		r.add("Ca", 2); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("c1b", 100.0); 
		r.add("IP3K-2Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("c2f", 100.0); 
		r.add("IP3K-2Ca", 1); 
		r.add("IP3", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("c2b", 80.0); 
		r.add("IP3K-2Ca-IP3", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("c3f", 20.0); 
		r.add("IP3K-2Ca-IP3", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("c3b", 0.0); 
		r.add("IP4", 1); 
		r.add("IP3K-2Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("c4f", 9.0); 
		r.add("IP3_5-phosphatase", 1); 
		r.add("IP3", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("c4b", 72.0); 
		r.add("IP5P-IP3", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("c5f", 18.0); 
		r.add("IP5P-IP3", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("c5b", 0.0); 
		r.add("IP2", 1); 
		r.add("IP3_5-phosphatase", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("d1f", 1000.0); 
		r.add("IP3Rec", 1); 
#ifndef FIX_d1
		r.add("IP3", 1); 
#else
		r.add("IP3_Fix", 1); 
#endif
		reac.push_back(r); 
	} 
	{ 
		Reac r("d1b", 25800.0); 
		r.add("IP3R-IP3", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("d2f", 8000.0); 
		r.add("IP3R-IP3", 1); 
#ifndef FIX_d2
		r.add("Ca", 1); 
#else
		r.add("Ca_Fix", 1); 
#endif
		reac.push_back(r); 
	} 
	{ 
		Reac r("d2b", 2000.0); 
		r.add("IP3R-open", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("d3f", 8.889); 
		r.add("IP3Rec", 1); 
		r.add("Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("d3b", 5.0); 
		r.add("IP3R-Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("d4f", 20.0); 
		r.add("IP3R-Ca", 1); 
		r.add("Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("d4b", 10.0); 
		r.add("IP3R-2Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("d5f", 40.0); 
		r.add("IP3R-2Ca", 1); 
		r.add("Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("d5b", 15.0); 
		r.add("IP3R-3Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("d6f", 60.0); 
		r.add("IP3R-3Ca", 1); 
		r.add("Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("d6b", 20.0); 
		r.add("IP3R-4Ca", 1); 
		reac.push_back(r); 
	} 
#ifndef CUT_e1
	{ 
		Reac r("e1f", 2250.0); 
		r.add("CaStore", 1); 
		r.add("IP3R-open", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("e1b", 450.0); 
		r.add("Ca", 1); 
		r.add("IP3R-open", 1); 
		reac.push_back(r); 
	} 
#endif
	{ 
		Reac r("e2f", 17147.0); 
		r.add("SERCA", 1); 
		r.add("Ca", 2); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("e2b", 1000.0); 
		r.add("SERCA-2Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("e3f", 250.0); 
		r.add("SERCA-2Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("e3b", 0.0); 
		r.add("SERCA", 1); 
		r.add("CaStore", 2); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("e4f", 1.25); 
		r.add("CaStore", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("e4b", 0.25); 
		r.add("Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("e5f", 25000.0); 
		r.add("PMCA", 1); 
		r.add("Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("e5b", 2000.0); 
		r.add("PMCA-Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("e6f", 500.0); 
		r.add("PMCA-Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("e6b", 0.0); 
		r.add("PMCA", 1); 
		r.add("Ca_ext", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("e7f", 93.827); 
		r.add("NCX", 1); 
		r.add("Ca", 2); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("e7b", 4000.0); 
		r.add("NCX-2Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("e8f", 1000.0); 
		r.add("NCX-2Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("e8b", 0.0); 
		r.add("NCX", 1); 
		r.add("Ca_ext", 2); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("e9f", 0.00166667); 
		r.add("Ca_ext", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("e9b", 0.166667); 
		r.add("Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("e10f", 0.1); 
		r.add("Calreticulin", 1); 
		r.add("CaStore", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("e10b", 200.0); 
		r.add("Calreticulin-Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("f1f", 1000.0); 
		r.add("MgGreen", 1); 
		r.add("Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("f1b", 19000.0); 
		r.add("MgGreenAst", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("f2f", 18.5); 
		r.add("parvalbumin", 1); 
		r.add("Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("f2b", 0.95); 
		r.add("PV-Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("f3f", 87.0); 
		r.add("Calbindin-D28k", 1); 
		r.add("Ca", 2); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("f3b", 11.275); 
		r.add("CB-Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("f4f", 10.0); 
		r.add("LowAffBuf", 1); 
		r.add("Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("f4b", 1000.0); 
		r.add("LAB-Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("f5f", 10.0); 
		r.add("LowAffBuf2", 1); 
		r.add("Ca", 2); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("f5b", 4000.0); 
		r.add("LAB2-Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("GluDecayf", 250.0); 
		r.add("Glu", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("GluDecayb", 0.0); 
		r.add("Glu-Decay", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("IP3Diffusionf", 980.39); 
		r.add("IP3_PSD", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("IP3Diffusionb", 19.608); 
		r.add("IP3", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("CaDiffusionf", 19.608); 
		r.add("Ca", 1); 
		reac.push_back(r); 
	} 
	{ 
		Reac r("CaDiffusionb", 980.39); 
		r.add("Ca_PSD", 1); 
		reac.push_back(r); 
	} 
	return reac;
}

map<String, Stim> getStim()
{
	map<String, Stim> stim;
	stim.insert(pair<String, Stim>("PFCa", Stim("PFCa", "Ca", 0.5, 0.001, 1500))); 
	stim.insert(pair<String, Stim>("PFGlu", Stim("PFGlu", "Glu", 0.5, 0.001, 300))); 
	stim.insert(pair<String, Stim>("CF", Stim("CF", "Ca", 0.5, 0.002, 5000))); 
	return stim;
}

#endif
#ifdef DEBUG

list<Mole> getMole()
{
	list<Mole> mole;
	mole.push_back(Mole("A", 1000, 0.1));
	mole.push_back(Mole("B", 3000, 0.1));
	mole.push_back(Mole("AB", 0, 0.1));
	mole.push_back(Mole("C", 2000, 0.1));
	mole.push_back(Mole("ACC", 0, 0.1));
	return mole;
}

list<Reac> getReac()
{
	list<Reac> reac;
	{
		Reac r("af", 1);
		r.add("A", 1);
		r.add("B", 1);
		reac.push_back(r);
	}
	{
		Reac r("ab", 0.2);
		r.add("AB", 1);
		reac.push_back(r);
	}
	{
		Reac r("bf", 0.2);
		r.add("A", 1);
		r.add("C", 2);
		reac.push_back(r);
	}
	{
		Reac r("bb", 10);
		r.add("ACC", 1);
		reac.push_back(r);
	}
	return reac;
}

map<String, Stim> getStim()
{
	map<String, Stim> stim;
	stim.insert(pair<String, Stim>("A", Stim("A", "A", 0.5, 0.05, 1000)));
	return stim;
}

#endif
