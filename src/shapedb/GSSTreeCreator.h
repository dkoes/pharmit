/*
 ShapeDB
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
 * GSSTreeCreator.h
 *
 *  Created on: Oct 13, 2011
 *      Author: dkoes
 *
 *      This class creates a GSS tree on disk.  It assumes the input
 *      does not fit in memory and behaves accordingly.
 *
 *      It takes an iterator over the input data, which must support
 *      an intersection (with a cube) method and a write to file method.
 *      This data is then converted to oct tree representations.  The
 *      data gets written to an indexed file while the oct trees get written
 *      to another file (storing indices to the object data).
 *
 *      The tree file is then clustered to create leaves which are
 *      appended to a leaf file as they are created.
 *
 *      The leaf nodes are similarily clustered to create a level of
 *      nodes, which are written to their own file.  These nodes are
 *      clustered into the next level's file and so on.
 *
 *      Clustering involves a top-down O(n) partitioning that splits the
 *      data until the a set small enough for an O(n^2) bottom-up packing
 *      to be performed.  Clusters are always packed to contain at least 2
 *      entries (this may be relaxed for leaves).
 *
 *      Once the final level is created, the internal nodes are all laid out
 *      in a separate file in depth-first order.  The leaves are ordered
 *      sequentially in a separate file.
 */

#ifndef GSSTREECREATOR_H_
#define GSSTREECREATOR_H_

#include <iostream>
#include <fstream>
#include <vector>

#include "GSSTypes.h"
#include "GSSTreeStructures.h"
#include "TopDownPartitioner.h"
#include "packers/Packer.h"

using namespace std;

#include "WorkFile.h"
#include "Timer.h"

//class for creating levels, follows the CM-tree bulk loading algorithm,
//but can be overridden to implement any arbitrary algorithm
class GSSLevelCreator
{
protected:
	const TopDownPartitioner *partitioner;
	const Packer *packer;

	//configuration settings
	unsigned nodePack;
	unsigned leafPack;

	//class vars used by nextlevelR
	unsigned packingSize;
	ostream *outNodes;
	ostream *outTrees;
	vector<file_index> *nodeIndices;
	vector<file_index> *treeIndices;
	virtual void createNextLevelR(TopDownPartitioner *P);

public:

	GSSLevelCreator() :
			partitioner(NULL), packer(NULL), nodePack(0), leafPack(0),
			packingSize(0), outNodes(NULL), outTrees(NULL)
	{
	}
	GSSLevelCreator(const TopDownPartitioner * part, const Packer *pack,
			unsigned np, unsigned lp) :
			partitioner(part), packer(pack), nodePack(np), leafPack(lp),
			packingSize(0), outNodes(NULL), outTrees(NULL)
	{
	}

	void initialize(const TopDownPartitioner * part, const Packer *pack,
			unsigned np = 32768, unsigned lp = 32768)
	{
		partitioner = part;
		packer = pack;
		nodePack = np;
		leafPack = lp;
	}

	virtual ~GSSLevelCreator()
	{
	}

	virtual void createNextLevel(DataViewer& data, ostream* nodefile,
			vector<file_index>& nodeindices, ostream* treefile,
			vector<file_index>& treeindices);

	unsigned getPack() const
	{
		return packer->getPack();
	}
};

class GSSTreeCreator
{
	WorkFile objects;
	WorkFile currenttrees;
	vector<file_index> treeindices;
	vector<file_index> objindices;

	vector<WorkFile> nodes;

	boost::filesystem::path dbpath;

	GSSLevelCreator *leveler;

	float dimension;
	float resolution;
	unsigned superNodeDepth;
	//some bookkeeping for analysis purposes
	unsigned numNodes;
	unsigned numLeaves;
	vector<unsigned> nodeContentDistribution;
	vector<unsigned> leafContentDistribution;
	file_index optimizeLevelsR(ostream& outnodes, ostream& outleaves,
			const GSSNodeCommon *n, unsigned level, file_index& lstart,
			file_index& lend);
	void optimizeLevels();

	void getNodesForSuperNode(const GSSInternalNode* root,
			vector<GSSInternalNode*>& newroots, unsigned curlevel,
			unsigned stoplevel);

public:
	GSSTreeCreator(GSSLevelCreator *l, unsigned sdepth = 3) :
			leveler(l), dimension(0), resolution(0), superNodeDepth(sdepth), numNodes(
					0), numLeaves(0)
	{
	}

	GSSTreeCreator() :
			leveler(NULL), dimension(0), resolution(0), superNodeDepth(3), numNodes(
					0), numLeaves(0)
	{

	}

	~GSSTreeCreator()
	{
		//workfiles must be explicitly cleared
		objects.clear();
		for (unsigned i = 0, n = nodes.size(); i < n; i++)
		{
			nodes[i].clear();
		}
	}

	float getDimension() const
	{
		return dimension;
	}

	float getResolution() const
	{
		return resolution;
	}

	bool create(boost::filesystem::path dir, boost::filesystem::path treedir,
			float dim,
			float res);

	//setup directories
	bool initialize(boost::filesystem::path dir, float dim, float res,
			GSSLevelCreator* l = NULL)
	{
		using namespace boost;
		dimension = dim;
		resolution = res;

		if (l)
			leveler = l;

		//create directory
		if (filesystem::exists(dir))
		{
			cerr << dir << " already exists.  Exiting\n";
			return false;
		}
		if (!filesystem::create_directory(dir))
		{
			cerr << "Unable to create database directory ";
			return false;
		}
		dbpath = dir;

		filesystem::path objfile = dbpath / "objs";
		string curtreesfile = filesystem::path(dbpath / "trees").string();

		//write out objects and trees
		objects.set(objfile.string().c_str());
		currenttrees.set(curtreesfile.c_str());
		treeindices.clear();
		objindices.clear();

		return true;
	}

	template<class Object> void addObject(const Object& obj)
	{
		objindices.push_back((file_index) objects.file->tellp());
		obj.write(*objects.file);

		//leaf object
		treeindices.push_back((file_index) currenttrees.file->tellp());
		MappableOctTree *tree = MappableOctTree::create(dimension, resolution,
				obj);
		tree->write(*currenttrees.file);
		delete tree;
	}

	bool createIndex();

	//return true if successful
	template<class Object, class ObjectIterator>
	bool create(boost::filesystem::path dir, ObjectIterator& itr,
			float dim, float res)
	{
		using namespace boost;
		initialize(dir, dim, res);
		Timer t;
		for (; itr; ++itr)
		{
			const Object& obj = *itr;
			addObject(obj);
		}

		cout << "Create/write trees\t" << t.elapsed() << "\n";
		t.restart();

		return createIndex();
	}

	//write out the object trees to the specified directory, with the object file
	//and also indices for reading back in later to save having to regenerate trees
	template<class Object, class ObjectIterator>
	bool createTreesOnly(boost::filesystem::path dir, ObjectIterator& itr,
			float dim,
			float res)
	{
		using namespace boost;
		dimension = dim;
		resolution = res;
		WorkFile currenttrees;
		//create directory
		if (filesystem::exists(dir))
		{
			cerr << dir << " already exists.  Exiting\n";
			return false;
		}
		if (!filesystem::create_directory(dir))
		{
			cerr << "Unable to create database directory ";
			return false;
		}
		dbpath = dir;

		filesystem::path objfile = dbpath / "objs";
		string curtreesfile = filesystem::path(dbpath / "trees").string();
		string tipath = filesystem::path(dbpath / "treeindices").string();
		string oipath = filesystem::path(dbpath / "objindices").string();

		Timer t;
		//write out objects and trees
		objects.set(objfile.string().c_str());
		currenttrees.set(curtreesfile.c_str());

		ofstream treeindices(tipath.c_str());
		ofstream objindices(oipath.c_str());

		if (!treeindices || !objindices)
			return false;

		unsigned cnt = 0;
		for (; itr; ++itr)
		{
			const Object& obj = *itr;
			file_index objindex = (file_index) objects.file->tellp();
			objindices.write((char*) &objindex, sizeof(file_index));
			obj.write(*objects.file);

			//leaf object
			file_index treeindex = (file_index) currenttrees.file->tellp();
			treeindices.write((char*) &treeindex, sizeof(file_index));
			MappableOctTree *tree = MappableOctTree::create(dim, res, obj);
			tree->write(*currenttrees.file);
			delete tree;
			cnt++;
		}
		currenttrees.clear();
		cout << "Create/write trees\t" << t.elapsed() << "\n";
		return true;
	}

	void printStats(ostream& out) const;
};

#endif /* GSSTREECREATOR_H_ */
