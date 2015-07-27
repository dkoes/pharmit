/*
 * OBMoleculeAnalytic.h
 *
 *  Created on: Feb 21, 2012
 *      Author: dkoes
 *
 *  This represents a molecule object.  Is supports intersection
 *  testing, file writing, and iteration over an input.
 *
 *  This version does not use grids.
 */

#ifndef OBMOLECULEA_H_
#define OBMOLECULEA_H_

#include <iostream>
#include <string>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/atom.h>
#include <cmath>
#include "MolSphere.h"
#include "Cube.h"
#include "MGrid.h"
#include "PMol.h"
#include "AnalyticCheckers.h"

using namespace OpenBabel;
using namespace std;

class OBAMolIterator;
class OBAMolOutput;


//interface between molecular shape and the gss tree
//defines iterator type, intersection, and file create/write
class OBAMolecule
{
	//a little icky, make mutable so intersect can bump intersecting spheres
	//to the front of the list; this ends up being more efficient than filtering
	//all the intersecting spheres
	mutable vector<SphereChecker> checkers;

	OBMol mol;
	float adjustAmount;
	float probe;

	OBMol& getMol()
	{
		return mol;
	}

	void set(const OBMol& m, float dimension, float resolution, float prb = 1.4,
			float adj = 0)
	{
		mol = m;
		adjustAmount = adj;
		probe = prb;
		vector<MolSphere> spheres;
		spheres.reserve(mol.NumAtoms());
		float mind = -dimension/2;
		float maxd = dimension/2;
		for (OBAtomIterator aitr = mol.BeginAtoms(); aitr != mol.EndAtoms();
				++aitr)
		{
			OBAtom* atom = *aitr;
			float r = etab.GetVdwRad(atom->GetAtomicNum())
											- adjustAmount;
			float x = atom->x();
			float y = atom->y();
			float z = atom->z();
			//omit sphere's outside of grid range; this may introduce
			//some very small artifacts at the edge of the grid
			if(x > mind-r && x < maxd+r &&
					y > mind - r && y < maxd+r &&
					z > mind - r && z < maxd+r)
			{
				spheres.push_back(MolSphere(x,y,z,r));
			}
		}

		checkers.clear();
		checkers.reserve(spheres.size());
		for (unsigned i = 0, n = spheres.size(); i < n; i++)
		{
			checkers.push_back(SphereChecker(spheres[i]));
			checkers.back().addNeighbors(spheres, i + 1, probe);
		}

		reverse(checkers.begin(), checkers.end());
	}

	friend class OBAMolOutput; //needs to call getMol
	friend class OBAMolIterator; //calls set

public:
	typedef OBAMolIterator iterator;

	OBAMolecule() :
			adjustAmount(0), probe(0)
	{
	}
	~OBAMolecule()
	{
	}

	bool intersects(const Cube& cube) const
	{
		for (unsigned i = 0, n = checkers.size(); i < n; i++)
		{
			if (checkers[i].intersectsCube(cube))
			{
				if (i > 0)
					swap(checkers[i], checkers[0]); //locality optimization

				return true;
			}
		}
		return false;
	}

	//return true if point is within object
	bool containsPoint(float x, float y, float z) const
	{
		for (unsigned i = 0, n = checkers.size(); i < n; i++)
		{
			if (checkers[i].containsPoint(x, y, z))
			{
				if (i > 0)
					swap(checkers[i], checkers[0]); //locality optimization
				return true;
			}
		}

		return false;
	}

	//computes a set of solitary grid points that represent the interaction between
	//this ligand and the provided receptor in some way
	void computeInteractionGridPoints(OBAMolecule& receptor, MGrid& grid,
			double interactionDist, double maxClusterDist,
			unsigned minClusterPoints, double interactionPointRadius);

	void write(ostream& out) const
	{
		OBMol m(mol);
		PMolCreator pmol(m, true);
		pmol.writeBinary(out);
	}

};

// use openbabel to read in from a file
class OBAMolIterator
{
	OBConversion inconv;
	OBMol mol;
	OBAMolecule currmolecule;
	bool valid;
	bool keepHydrogens;
	float dimension;
	float resolution;
	float probe;

	void readOne()
	{
		valid = inconv.Read(&mol);
		if (!keepHydrogens)
			mol.DeleteHydrogens();
		currmolecule.set(mol, dimension, resolution, probe);
	}
public:
	OBAMolIterator(const string& fname, float dim, float res, bool keepH,
			float prb) :
			valid(true), keepHydrogens(keepH), dimension(dim), resolution(res), probe(
					prb)
	{
		if (fname.size() == 0)
		{
			valid = false;
			return;
		}
		inconv.SetInFormat(inconv.FormatFromExt(fname));
		valid = inconv.ReadFile(&mol, fname);
		if (!keepHydrogens)
			mol.DeleteHydrogens();
		mol.SetHybridizationPerceived(); //otherwise segfault trying to gethyb
		currmolecule.set(mol, dimension, resolution, probe);
	}

	//validity check
	operator bool() const
	{
		return valid;
	}

	//current mol
	const OBAMolecule& operator*() const
	{
		return currmolecule;
	}

	void operator++()
	{
		readOne();
	}

};


#endif /* OBMOLECULEA_H_ */
