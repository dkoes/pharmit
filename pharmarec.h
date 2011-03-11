/*
Pharmer: Efficient and Exact 3D Pharmacophore Search
Copyright (C) 2011  David Ryan Koes and the University of Pittsburgh

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

/*
 * pharmarec.h
 *
 *  Created on: July 30, 2010
 *      Author: dkoes
 *
 *      Routines for pharmacophore recognition.
 */

#ifndef PHARMAREC_H_
#define PHARMAREC_H_

#include <iostream>
#include <boost/unordered_map.hpp>
#include "boost/tuple/tuple.hpp"
#include <float.h>
#include <json/json.h>
#include <openbabel/mol.h>
#include <openbabel/math/vector3.h>
#include <openbabel/obconversion.h>
#include "CommandLine2/CommandLine.h"

using namespace std;
using namespace boost;
using namespace OpenBabel;

extern cl::opt<bool> DKoesTest;

//somewhat arbitrary numerical threshold for close enoughness
#define THRESHOLD (.0001)

//generate a vector corresponding to pharma point (hbonds, aromatic?)
struct PharmaPoint;
typedef void (*genPointVectorFn)(const vector<int>& atoms_indexes,
		const OBMol& mol, PharmaPoint& pnt);

//description of interaction features
struct PharmaInteract
{
	unsigned complement; //index of interacting point
	float maxDist; //max distance that a feature can be to still be considered interacting
	unsigned minMatch; //minimum number of feature that must be close

	PharmaInteract(unsigned c, float md, unsigned mm): complement(c), maxDist(md), minMatch(mm) {}
	PharmaInteract(): complement(0), maxDist(0), minMatch(0) {}
};

//a description of a pharmacophore class
struct Pharma {
	string name;
	int atomic_number_label;
	unsigned index; //position in pharmas array
	vector<OBSmartsPattern> smarts;
	float defaultSearchRadius;
	genPointVectorFn getVectors;
	float clusterLimit;

	Pharma(): atomic_number_label(0), index(-1), defaultSearchRadius(0), getVectors(NULL), clusterLimit(0)
	{

	}

	Pharma(int indx, const char *n, const char** sm, int atomic, float r = .5, float cl = 0, unsigned nb =1, float trfr = 0):
		name(n), atomic_number_label(atomic),
				index(indx), defaultSearchRadius(r), getVectors(NULL), clusterLimit(cl){
		//const char * smiles for easier initialization
		if(sm != NULL)
		{
			while(*sm != NULL)
			{
				smarts.push_back(OBSmartsPattern());
				smarts.back().Init(*sm);
				sm++;
			}
		}
		setVectorFn();
	}

	Pharma(int indx, const string& n, const vector<string>& sm, int atomic, float r = .5, float cl = 0, unsigned nb =1, float trfr = 0):
		name(n), atomic_number_label(atomic),
				index(indx), defaultSearchRadius(r), getVectors(NULL), clusterLimit(cl) {
		unsigned ns = sm.size();
		smarts.resize(ns);
		for(unsigned i = 0; i < ns; i++)
		{
			smarts[i].Init(sm[i]);
		}
		setVectorFn();
	}

	bool operator==(const Pharma& rhs) const;
	bool operator!=(const Pharma& rhs) const { return !(*this == rhs); }
	//function for determining vectors of a point; with no arguments
	//uses built-in defaults
	void setVectorFn(genPointVectorFn fn = NULL);
};

//a collection of pharmacophore descriptors
//PharmaPoints reference back to these descriptors so the creating Pharmas
//should not be allowed to go out of scope before the pharmapoints
class Pharmas
{
private:
	//supporting copying operations is askign for trouble - there
	//really should only be one instance of pharmas that is either programmatically
	//initialized with a vector of Pharmas or by reading in pharmas from a file
	Pharmas(const Pharmas&);
	Pharmas& operator=(const Pharmas&);

protected:
	Pharma *pharmas;
	unsigned numPharmas;
	unordered_map<string, unsigned> nameLookup;
	unordered_map<unsigned, unsigned> anumLookup;
	void initialize(const vector<Pharma>& ps);

public:

	Pharmas(): pharmas(NULL), numPharmas(0) {}

	//note that ps gets copied: pointers don't stay the same
	Pharmas(const vector<Pharma>& ps): pharmas(NULL), numPharmas(0)
	{
		initialize(ps);
	}

	virtual ~Pharmas()
	{
		delete [] pharmas;
	}


	const Pharma* operator[](unsigned i) const { return &pharmas[i]; }

	bool operator==(const Pharmas& rhs) const;
	bool operator!=(const Pharmas& rhs) const { return !(*this == rhs); }

	unsigned size() const { return numPharmas; }

	const Pharma* pharmaFromAtomicNumber(int anum) const;
	const Pharma* pharmaFromName(const string& name) const;

	bool read(istream& in);
	void write(ostream& out) const;
};

extern const vector<Pharma> defaultPharmaVec;

//a single pharmacophore point
struct PharmaPoint {
	double x;
	double y;
	double z;
	vector<vector3> vecs; //directionality, can have multiple
	unsigned size;
	const Pharma *pharma; //pharmacophore descriptor

	//query attributes
	double radius;
	float vecpivot;
	enum PointRequirements {Required, Optional, NotPresent};
	PointRequirements requirements;
	unsigned minSize;
	unsigned maxSize;

	PharmaPoint() :
		x(0), y(0), z(0), size(0), pharma(NULL), radius(0), vecpivot(0), requirements(Required), minSize(0), maxSize(0) {
	}

	PharmaPoint(const Pharma* p) :
		x(0), y(0), z(0), size(0), pharma(p), radius(p->defaultSearchRadius), vecpivot(0), requirements(Required), minSize(0), maxSize(0)  {
	}

	PharmaPoint(double xc, double yc, double zc, unsigned sz, const Pharma* p) :
		x(xc), y(yc), z(zc), size(),pharma(p), radius(p->defaultSearchRadius), vecpivot(0), requirements(Required), minSize(0), maxSize(0)  {
	}

	double radiusWeight() const
	{
		return 1.0/(radius*radius);
	}
	//essentially in the same spot
	bool sameLocation(const PharmaPoint& rhs) const {
		if (fabs(x - rhs.x) > THRESHOLD)
			return false;
		if (fabs(y - rhs.y) > THRESHOLD)
			return false;
		if (fabs(z - rhs.z) > THRESHOLD)
			return false;
		return true;
	}

	const char* requirementStr() const
	{
		switch(requirements)
		{
		case Required:
			return "required";
		case Optional:
			return "optional";
		case NotPresent:
			return "notpresent";
		}
		return NULL;
	}

#define SQR(x) ((x)*(x))

	//return distance between a and b
	static double pharmaDist(const PharmaPoint& a, const PharmaPoint& b) {
		return sqrt(SQR(a.x-b.x) + SQR(a.y-b.y) + SQR(a.z-b.z));
	}
	friend ostream& operator<<(ostream &stream, const PharmaPoint& obj);

	bool read(const Pharmas& pharmas, istream &in);
};


//read from a json formatted stream
extern bool readPharmaPointsJSON(const Pharmas& pharma, Json::Value& data, vector<PharmaPoint>& points);

//identify all pharma points in mol
extern void getPharmaPoints(const Pharmas& pharmas, OBMol& mol, vector<PharmaPoint>& points);
//identify all pharma points in a multi-conformer molecule
extern void getPharmaPointsMC(const Pharmas& pharmas, OBMol& mol, vector< vector<PharmaPoint> >& points);

//accelerated pharma point recognition for proteins
extern void getProteinPharmaPoints(const Pharmas& pharmas, OBMol& protein, vector<PharmaPoint>& points);

extern void getInteractionPoints(const Pharmas& pharmas, OBMol& receptor, OBMol& ligand,
		vector<PharmaPoint>& points, vector<PharmaPoint>& screenedout);

//translate a point vector into json
extern bool convertPharmaJson(Json::Value& root, const vector<PharmaPoint>& points);


//add artificial pharma "atoms" to mol
extern void addPharmaPoints(OBMol& mol, vector<PharmaPoint>& points);

//extract pharmacophore points into json from moldata in format
extern bool jsonPharmaQuery(const Pharmas& pharmas, Json::Value& root,
		const string& moldata, OBFormat *format, const string& recdata, OBFormat *rformat);


#endif /* PHARMAREC_H_ */
