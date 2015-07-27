/*
 * GSSTreeSearcher.h
 *
 *  Created on: Oct 18, 2011
 *      Author: dkoes
 *
 *      Class for loading and searching a GSSTree.
 */

#ifndef GSSTREESEARCHER_H_
#define GSSTREESEARCHER_H_

#include "GSSTreeStructures.h"
#include "molecules/Molecule.h"
#include "MemMapped.h"
#include "MappableOctTree.h"
#include "molecules/ResultMolecules.h"

using namespace std;

class GSSTreeSearcher
{
	MemMapped objects; //memory mapped objects
	MemMapped internalNodes;
	MemMapped leaves;

	bool verbose; //for debugging
	unsigned total;
	float dimension;
	float resolution;

	void findTweeners(const GSSInternalNode* node, const MappableOctTree* min,
			const MappableOctTree* max, const MappableOctTree* orig,
			vector<result_info>& res,
			unsigned level, bool computeDist);
	void findTweeners(const GSSLeaf* node, const MappableOctTree* min,
			const MappableOctTree* max, const MappableOctTree* orig,
			vector<result_info>& res,
			bool computeDist);

	struct ObjDist
	{
		file_index objpos;
		double dist;

		bool operator<(const ObjDist& rhs) const
				{
			return dist < rhs.dist;
		}
	};

	//retain best k objects that are better than threshold
	class TopObj
	{
		vector<ObjDist> objs;
		unsigned k;
		double thresh;
	public:
		TopObj(unsigned _k, double t = HUGE_VAL): k(_k), thresh(t)
		{
			if(k == 0) //then no limit
				k = UINT_MAX;
			else
				objs.reserve(k+1);
		}

		void add(file_index pos, double dist);

		double worst() const
		{
			if (objs.size() < k)
				return thresh;
			else
				return objs.back().dist;
		}

		unsigned size() const
		{
			return objs.size();
		}

		const ObjDist& operator[](unsigned i) const
				{
			return objs[i];
		}

	};

	void findNearest(const GSSInternalNode* node, const MappableOctTree* obj,
			TopObj& res, unsigned level);
	void findNearest(const GSSLeaf* node, const MappableOctTree* obj,
			TopObj& res);

	void findNearest(const GSSInternalNode* node, const MappableOctTree* minobj,
			const MappableOctTree* maxobj,
			TopObj& res, unsigned level);
	void findNearest(const GSSLeaf* node, const MappableOctTree* minobj,
			const MappableOctTree* maxobj,
			TopObj& res);

	unsigned fitsCheck;
	unsigned nodesVisited;
	vector<unsigned> levelCnts;
	vector<unsigned> usefulLevelCnts;
	vector<unsigned> maxlevelCnts;
	unsigned leavesVisited;
	unsigned fullLeaves;
	bool fitsInbetween(const MappableOctTree *MIV, const MappableOctTree *MSV,
			const MappableOctTree *min, const MappableOctTree *max);

public:
	typedef shared_ptr<const MappableOctTree> ObjectTree;

	GSSTreeSearcher(bool v = false) :
			verbose(v), total(0)
	{
	}

	bool load(const filesystem::path& dbpath);
	void clear();
	~GSSTreeSearcher();

	unsigned size() const
	{
		return total;
	}

	//return everything with a shape between smallTree and bigTree
	void dc_search(ObjectTree smallTree, ObjectTree bigTree, ObjectTree refTree,
			bool loadObjs,
			Results& res);

	//linear scan
	void dc_scan_search(ObjectTree smallTree, ObjectTree bigTree,
			ObjectTree refTree, bool loadObjs,
			Results& res);

	//return k objects closest to obj and better than threshold
	void nn_search(ObjectTree objTree, unsigned k, double thresh, bool loadObjs,
			Results& res);

	//compute scores for all molecules in database
	void nn_scan(ObjectTree objTree, bool loadObjs,
			Results& res);

	//return k objects closes to small/big obj using shapeDistance
	void nn_search(ObjectTree smallTree, ObjectTree bigTree, unsigned k,
			bool loadObjs, Results& res);

	//same as above, but evaluate entire database
	void nn_scan(ObjectTree smallTree, ObjectTree bigTree, bool loadObjs,
			Results& res);

	//given an object, return a tree
	template<class Object>
	ObjectTree createTreeFromObject(const Object& obj,
			float shrink = 0, bool invert = false)
	{
		MappableOctTree *objTree = MappableOctTree::create(dimension,
				resolution, obj);

		if (shrink > 0)
		{
			//reduce the object
			MGrid grid;
			objTree->makeGrid(grid, resolution);
			grid.shrink(shrink);
			free(objTree);
			objTree = MappableOctTree::createFromGrid(grid);
		}
		else if (shrink < 0) //grow
		{
			//expand the object
			MGrid grid;
			objTree->makeGrid(grid, resolution);
			grid.grow(-shrink);
			free(objTree);
			objTree = MappableOctTree::createFromGrid(grid);
		}

		if (invert) //treat as excluded vol
			objTree->invert();

		return shared_ptr<const MappableOctTree>(objTree, free);
	}

	float getDimension() const
	{
		return dimension;
	}
	float getResolution() const
	{
		return resolution;
	}
};

#endif /* GSSTREESEARCHER_H_ */
