/*
 * RDMoleculeAnalytic.h
 *
 *  Created on: Jun 2, 2014
 *      Author: dkoes
 *
 *  Support for RDKit molecules
 */

#ifndef RDMOLECULEANALYTIC_H_
#define RDMOLECULEANALYTIC_H_



#include <iostream>
#include <string>

#include <cmath>
#include "MolSphere.h"
#include "Cube.h"
#include "MGrid.h"
#include "PMol.h"
#include "AnalyticCheckers.h"

#include <GraphMol/ROMol.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
using namespace RDKit;
using namespace std;

class RDMolIterator;
class RDMolOutput;


//interface between molecular shape and the gss tree
//defines iterator type, intersection, and file create/write
//this is based off of OBMolecule
class RDMolecule
{
	//a little icky, make mutable so intersect can bump intersecting spheres
	//to the front of the list; this ends up being more efficient than filtering
	//all the intersecting spheres
	mutable vector<SphereChecker> checkers;

	ROMol *mol;
	float adjustAmount;
	float probe;

	ROMol* getMol() //RDMol will deallocate returned memory on next set
	{
		return mol;
	}

	//takes ownership of mol ptr
	void set(ROMol* m, float dimension, float resolution, float prb = 1.4,
			float adj = 0)
	{
		const PeriodicTable *tbl = PeriodicTable::getTable();
		if(mol != NULL) delete mol;
		mol = m;
		const Conformer& conf = mol->getConformer();
		adjustAmount = adj;
		probe = prb;
		vector<MolSphere> spheres;
		spheres.reserve(mol->getNumAtoms());
		float mind = -dimension/2;
		float maxd = dimension/2;
		for (ROMol::AtomIterator aitr = mol->beginAtoms(); aitr != mol->endAtoms();
				++aitr)
		{
			Atom* atom = *aitr;
			float r = tbl->getRvdw(atom->getAtomicNum())
											- adjustAmount;
			RDGeom::Point3D pt = conf.getAtomPos(atom->getIdx());
			float x = pt.x;
			float y = pt.y;
			float z = pt.z;
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

	friend class RDMolOutput; //needs to call getMol
	friend class RDMolIterator; //calls set

public:
	typedef RDMolIterator iterator;

	RDMolecule() :
			adjustAmount(0), probe(0), mol(NULL)
	{
	}
	~RDMolecule()
	{
		if(mol) delete mol;
		mol = NULL;
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


	void write(ostream& out) const
	{
		PMolCreator pmol(*mol, true);
		pmol.writeBinary(out);
	}

};

// use rdkit to read in from a file
class RDMolIterator
{
	ForwardSDMolSupplier *sdsup;
	std::vector<std::ios*> streams;

	ROMol mol;
	RDMolecule currmolecule;
	bool valid;
	bool keepHydrogens;
	float dimension;
	float resolution;
	float probe;

	void readOne()
	{
		valid = !sdsup->atEnd();
		if(valid)
		{
			ROMol *mol = sdsup->next();
			if(!mol)
				valid = false;
			else
				currmolecule.set(mol, dimension, resolution, probe);
		}
	}
public:
	RDMolIterator(const string& fname, float dim, float res, bool keepH,
			float prb) : sdsup(NULL),
			valid(true), keepHydrogens(keepH), dimension(dim), resolution(res), probe(
					prb)
	{
		if (fname.size() == 0)
		{
			valid = false;
			return;
		}

		//open streams, checking for gzip
		std::ifstream *uncompressed_inmol = new std::ifstream(fname.c_str());
		streams.push_back(uncompressed_inmol);
		boost::iostreams::filtering_stream<boost::iostreams::input> *inmol = new boost::iostreams::filtering_stream<boost::iostreams::input>();
		streams.push_back((std::istream*) inmol);

		std::string::size_type pos = fname.rfind(".gz");
		if (pos != std::string::npos)
		{
			inmol->push(boost::iostreams::gzip_decompressor());
		}
		inmol->push(*uncompressed_inmol);


		sdsup = new ForwardSDMolSupplier(inmol,false,false,!keepH);
		if(sdsup->atEnd())
		{
			valid = false;
			return;
		}
		ROMol *mol = sdsup->next();
		currmolecule.set(mol, dimension, resolution, probe);
	}

	~RDMolIterator()
	{
		if(sdsup) delete sdsup;
		for (unsigned i = 0, n = streams.size(); i < n; i++)
		{
			delete streams[i];
		}
		streams.clear();
	}
	//validity check
	operator bool() const
	{
		return valid;
	}

	//current mol
	const RDMolecule& operator*() const
	{
		return currmolecule;
	}

	void operator++()
	{
		readOne();
	}

};


#endif /* RDMOLECULEANALYTIC_H_ */
