/*
 * GSSTreeSearcher.cpp
 *
 *  Created on: Oct 18, 2011
 *      Author: dkoes
 */

#include "GSSTreeSearcher.h"
#include <sys/types.h>
#include <sys/mman.h>
#include "MappableOctTree.h"
#include "Timer.h"
#include "ShapeDistance.h"

//load a gsstree database by mmapping files, return true if successfull
bool GSSTreeSearcher::load(const filesystem::path& dbpath)
{
	clear();
	//read in info
	filesystem::path infofile = dbpath / "info";
	ifstream info(infofile.string().c_str());
	if (!info)
		return false;
	unsigned levels;

	info >> dimension >> resolution >> levels >> total;

	//memory map files form db directory
	filesystem::path nodepath = dbpath / "nodes";
	internalNodes.map(nodepath.string(), true, true);

	filesystem::path leavespath = dbpath / "leaves";
	leaves.map(leavespath.string(), true, true);

	filesystem::path objfile = dbpath / "objs";
	if (!objects.map(objfile.string(), true, false))
	{
		clear();
		return false;
	}

	return true;
}

void GSSTreeSearcher::clear()
{
	objects.clear();
	internalNodes.clear();
	leaves.clear();
}

GSSTreeSearcher::~GSSTreeSearcher()
{
	clear();
}

//find all the shapes in the database that lie bewtween smallObj and bigObj
//if invertBig is set, than treat as an excluded volume
//if refobjtree is provided, than compute volume overlap with refobjtree as score
void GSSTreeSearcher::dc_search(ObjectTree smallobjTree,
		ObjectTree bigobjTree, ObjectTree refobjTree, bool loadObjs, Results& res)
{
	const MappableOctTree* smallTree = smallobjTree.get();
	const MappableOctTree* bigTree = bigobjTree.get();
	const MappableOctTree* origTree = refobjTree.get();
	res.clear();

	Timer t;
	vector<result_info> respos;
	fitsCheck = 0;
	fullLeaves = 0;
	nodesVisited = 0;
	leavesVisited = 0;
	levelCnts.clear();
	if (internalNodes.size() > 0)
	{
		const GSSInternalNode* root = (GSSInternalNode*) internalNodes.begin();
		findTweeners(root, smallTree, bigTree, origTree, respos, 0, loadObjs);
	}
	else
	{
		//very small tree with just a leaf
		const GSSLeaf* leaf = (GSSLeaf*) leaves.begin();
		findTweeners(leaf, smallTree, bigTree, origTree, respos, loadObjs);
	}

	Timer objload;
	if (loadObjs)
	{
		//extract objects in sequential order
		sort(respos.begin(), respos.end());
		for (unsigned i = 0, n = respos.size(); i < n; i++)
		{
			const char * addr = objects.begin() + respos[i].pos;
			res.add(addr, respos[i].val);
		}
	}

	if (verbose)
	{
		cout << "Found " << respos.size() << " objects out of " << total
				<< " in " << t.elapsed() << " s (" << objload.elapsed()
				<< " objload) with " << fitsCheck << " checks " << nodesVisited
				<< " nodes " << leavesVisited << " leaves " << fullLeaves
				<< " full leaves\n";
		for (unsigned i = 0, n = levelCnts.size(); i < n; i++)
		{
			cout << " level " << i << ": " << levelCnts[i] << " "
					<< maxlevelCnts[i] << " " << usefulLevelCnts[i] << "\n";
		}
	}
}

void GSSTreeSearcher::TopObj::add(file_index pos, double dist)
{
	ObjDist x =
			{ pos, dist };
	if(dist < thresh)
	{
		objs.insert(lower_bound(objs.begin(), objs.end(), x), x);

		if (objs.size() > k)
			objs.resize(k);
	}
}

struct ScoredChild
{
	const GSSInternalNode::Child * child;
	double score;
	double min;

	ScoredChild() :
			child(NULL), score(HUGE_VAL), min(HUGE_VAL)
	{
	}
	ScoredChild(const GSSInternalNode::Child *ch, double s, double m) :
			child(ch), score(s), min(m)
	{
	}

	bool operator<(const ScoredChild& rhs) const
			{
		return score < rhs.score;
	}
};

//explore children to find closest value
void GSSTreeSearcher::findNearest(const GSSInternalNode* node,
		const MappableOctTree* obj, TopObj& res, unsigned level)
{
	nodesVisited++;
	if (levelCnts.size() <= level)
	{
		levelCnts.resize(level + 1, 0);
		usefulLevelCnts.resize(level + 1, 0);
		maxlevelCnts.resize(level + 2, 0);
	}
	levelCnts[level]++;

	unsigned n = node->size();
	vector<ScoredChild> children;
	children.reserve(n);
	for (unsigned i = 0; i < n; i++)
	{
		const GSSInternalNode::Child *child = node->getChild(i);
		float min = 0, max = 0;
		float score = searchVolumeDist(obj, child->getMIV(), child->getMSV(),
				min, max);
		fitsCheck++;
		if (min < res.worst())
		{
			//there's hope of finding something
			children.push_back(ScoredChild(child, score, min));

		}
	}

	maxlevelCnts[level + 1] += n;

	double oldworst = res.worst();

	sort(children.begin(), children.end());
	for (unsigned i = 0, nc = children.size(); i < nc; i++)
	{
		if (children[i].min < res.worst())
		{
			const GSSInternalNode::Child *child = children[i].child;
			if (child->isLeafPosition())
			{
				const GSSLeaf* next = (const GSSLeaf*) (leaves.begin()
						+ child->position());
				findNearest(next, obj, res);
			}
			else
			{
				const GSSInternalNode* next =
						(const GSSInternalNode*) (internalNodes.begin()
								+ child->position());
				findNearest(next, obj, res, level + 1);
			}
		}

	}
	if (res.worst() < oldworst)
		usefulLevelCnts[level]++;
}

//add nearest neighbors to res if appropriate
void GSSTreeSearcher::findNearest(const GSSLeaf* node,
		const MappableOctTree* obj, TopObj& res)
{
	leavesVisited++;
	unsigned cnt = 0;
	for (unsigned i = 0, n = node->size(); i < n; i++)
	{
		const GSSLeaf::Child *child = node->getChild(i);
		double dist = volumeDist(obj, &child->tree);
		fitsCheck++;
		if (dist < res.worst())
		{
			res.add(child->object_pos, dist);
			cnt++;
		}
	}
	if (cnt == node->size())
		fullLeaves++;
}

//add nearest neighbors to res if appropriate
void GSSTreeSearcher::findNearest(const GSSLeaf* node,
		const MappableOctTree* smallTree, const MappableOctTree* bigTree,
		TopObj& res)
{
	leavesVisited++;
	unsigned cnt = 0;
	for (unsigned i = 0, n = node->size(); i < n; i++)
	{
		const GSSLeaf::Child *child = node->getChild(i);
		double dist = shapeDistance(&child->tree, &child->tree, smallTree,
				bigTree);
		fitsCheck++;
		if (dist < res.worst())
		{
			res.add(child->object_pos, dist);
			cnt++;
		}
	}
	if (cnt == node->size())
		fullLeaves++;
}

//return k objects closes to small/big obj using shapeDistance
//in order to do this efficiently we need to be able to calculate bounds
//on the best possible score when we're high up in the tree - since shapeDistance
//is generic and doesn't provide bounds, this is just implemented as a linear scan
//if it is determined that a specific shapedistance function combined with nn_search
//is highly desireable, then we should be smarter
void GSSTreeSearcher::nn_search(ObjectTree smallobjTree,
		ObjectTree bigobjTree, unsigned k,
		bool loadObjs, Results& res)
{
	const MappableOctTree* smallTree = smallobjTree.get();
	const MappableOctTree* bigTree = bigobjTree.get();
	res.clear();
	TopObj ret(k);
	Timer t;

	const GSSLeaf* leaf = (GSSLeaf*) leaves.begin();
	const GSSLeaf* end = (GSSLeaf*) leaves.end();

	for (; leaf != end; leaf = (const GSSLeaf*) ((char*) leaf + leaf->bytes()))
	{
		findNearest(leaf, smallTree, bigTree, ret);
	}

	Timer objload;

	if (loadObjs)
	{
		//extract objects in distance order
		for (unsigned i = 0, n = ret.size(); i < n; i++)
		{
			if (verbose)
			{
				cout << i << " " << ret[i].dist << "\n";
			}
			const char * addr = objects.begin() + ret[i].objpos;
			res.add(addr, ret[i].dist);
		}
	}

	if (verbose)
	{
		cout << "Found " << ret.size() << " objects out of " << total << " in "
				<< t.elapsed() << " s (" << objload.elapsed()
				<< " objload) with " << fitsCheck << " checks " << nodesVisited
				<< " nodes " << leavesVisited << " leaves " << fullLeaves
				<< " full leaves\n";
		for (unsigned i = 0, n = levelCnts.size(); i < n; i++)
		{
			cout << " level " << i << ": " << levelCnts[i] << " "
					<< maxlevelCnts[i] << " " << usefulLevelCnts[i] << "\n";
		}
	}
}

void GSSTreeSearcher::nn_scan(ObjectTree smallobjTree, ObjectTree bigobjTree,
		bool loadObjs, Results& res)
{
	//unpack pointers for convenience
	const MappableOctTree* smallTree = smallobjTree.get();
	const MappableOctTree* bigTree = bigobjTree.get();

	res.clear();
	const GSSLeaf* leaf = (GSSLeaf*) leaves.begin();
	const GSSLeaf* end = (GSSLeaf*) leaves.end();

	Timer t;
	fitsCheck = 0;
	fullLeaves = 0;
	nodesVisited = 0;
	leavesVisited = 0;
	levelCnts.clear();
	vector<result_info> respos;
	respos.reserve(size());
	for (; leaf != end; leaf = (const GSSLeaf*) ((char*) leaf + leaf->bytes()))
	{
		leavesVisited++;
		for (unsigned i = 0, n = leaf->size(); i < n; i++)
		{
			const GSSLeaf::Child *child = leaf->getChild(i);
			double dist = shapeDistance(&child->tree, &child->tree, smallTree,
					bigTree);
			fitsCheck++;
			respos.push_back(result_info(child->object_pos, dist));
		}
	}

	Timer objload;

	if (loadObjs)
	{
		//extract objects in sequential order
		sort(respos.begin(), respos.end());
		res.reserve(respos.size());
		for (unsigned i = 0, n = respos.size(); i < n; i++)
		{
			const char * addr = objects.begin() + respos[i].pos;
			res.add(addr, respos[i].val);
		}
	}

	if (verbose)
	{
		cout << "Scanned " << respos.size() << " objects out of " << total
				<< " in "
				<< t.elapsed() << " s (" << objload.elapsed()
				<< " objload) with " << fitsCheck << " checks " << nodesVisited
				<< " nodes " << leavesVisited << " leaves " << fullLeaves
				<< " full leaves\n";
		for (unsigned i = 0, n = levelCnts.size(); i < n; i++)
		{
			cout << " level " << i << ": " << levelCnts[i] << " "
					<< maxlevelCnts[i] << " " << usefulLevelCnts[i] << "\n";
		}
	}
}

void GSSTreeSearcher::nn_search(ObjectTree objectTree,	unsigned k, double thresh, bool loadObjs,
		Results& res)
{
	const MappableOctTree* objTree = objectTree.get();
	TopObj ret(k,thresh);
	Timer t;
	fitsCheck = 0;
	fullLeaves = 0;
	nodesVisited = 0;
	leavesVisited = 0;
	levelCnts.clear();
	if (internalNodes.size() > 0)
	{
		const GSSInternalNode* root = (GSSInternalNode*) internalNodes.begin();
		findNearest(root, objTree, ret, 0);
	}
	else
	{
		//very small tree with just a leaf
		const GSSLeaf* leaf = (GSSLeaf*) leaves.begin();
		findNearest(leaf, objTree, ret);
	}

	Timer objload;

	if (loadObjs)
	{
		//extract objects in distance order
		for (unsigned i = 0, n = ret.size(); i < n; i++)
		{
			if (verbose)
			{
				cout << i << " " << ret[i].dist << "\n";
			}
			const char * addr = objects.begin() + ret[i].objpos;
			res.add(addr, ret[i].dist);
		}
	}

	if (verbose)
	{
		cout << "Found " << ret.size() << " objects out of " << total << " in "
				<< t.elapsed() << " s (" << objload.elapsed()
				<< " objload) with " << fitsCheck << " checks " << nodesVisited
				<< " nodes " << leavesVisited << " leaves " << fullLeaves
				<< " full leaves\n";
		for (unsigned i = 0, n = levelCnts.size(); i < n; i++)
		{
			cout << " level " << i << ": " << levelCnts[i] << " "
					<< maxlevelCnts[i] << " " << usefulLevelCnts[i] << "\n";
		}
	}
}

void GSSTreeSearcher::nn_scan(ObjectTree objectTree, bool loadObjs, Results& res)
{
	const MappableOctTree* objTree = objectTree.get();

	res.clear();

	const GSSLeaf* leaf = (GSSLeaf*) leaves.begin();
	const GSSLeaf* end = (GSSLeaf*) leaves.end();

	Timer t;
	fitsCheck = 0;
	fullLeaves = 0;
	nodesVisited = 0;
	leavesVisited = 0;
	levelCnts.clear();
	vector<result_info> respos;
	respos.reserve(size());
	for (; leaf != end; leaf = (const GSSLeaf*) ((char*) leaf + leaf->bytes()))
	{
		leavesVisited++;
		for (unsigned i = 0, n = leaf->size(); i < n; i++)
		{
			const GSSLeaf::Child *child = leaf->getChild(i);
			double dist = volumeDist(objTree, &child->tree);
			fitsCheck++;
			respos.push_back(result_info(child->object_pos, dist));
		}
	}

	Timer objload;

	if (loadObjs)
	{
		//extract objects in sequential order
		sort(respos.begin(), respos.end());
		res.reserve(respos.size());
		for (unsigned i = 0, n = respos.size(); i < n; i++)
		{
			const char * addr = objects.begin() + respos[i].pos;
			res.add(addr, respos[i].val);
		}
	}

	if (verbose)
	{
		cout << "Scanned " << respos.size() << " objects out of " << total
				<< " in "
				<< t.elapsed() << " s (" << objload.elapsed()
				<< " objload) with " << fitsCheck << " checks " << nodesVisited
				<< " nodes " << leavesVisited << " leaves " << fullLeaves
				<< " full leaves\n";
		for (unsigned i = 0, n = levelCnts.size(); i < n; i++)
		{
			cout << " level " << i << ": " << levelCnts[i] << " "
					<< maxlevelCnts[i] << " " << usefulLevelCnts[i] << "\n";
		}
	}
}

//return true if the object(s) represented by MIV/MSV might fit in between min and max
bool GSSTreeSearcher::fitsInbetween(const MappableOctTree *MIV,
		const MappableOctTree *MSV, const MappableOctTree *min,
		const MappableOctTree *max)
{
	fitsCheck++;
	//the MSV must completely enclose min
	if (!min->containedIn(MSV))
		return false;
	//MIV must be completely enclosed by max
	if (!MIV->containedIn(max))
		return false;

	return true;
}

//find everyting between small and big using linear scan
void GSSTreeSearcher::dc_scan_search(ObjectTree smallobjTree,
		ObjectTree bigobjTree, ObjectTree origobjTree, bool loadObjs, Results& res)
{
	const MappableOctTree* smallTree = smallobjTree.get();
	const MappableOctTree* bigTree = bigobjTree.get();
	const MappableOctTree* origTree = origobjTree.get();
	res.clear();

	const GSSLeaf* leaf = (GSSLeaf*) leaves.begin();
	const GSSLeaf* end = (GSSLeaf*) leaves.end();

	Timer t;
	vector<result_info> respos;
	for (; leaf != end; leaf = (const GSSLeaf*) ((char*) leaf + leaf->bytes()))
	{
		findTweeners(leaf, smallTree, bigTree, origTree, respos, loadObjs);
	}

	Timer objload;

	if (loadObjs)
	{
		//extract objects in sequential order
		sort(respos.begin(), respos.end());
		res.reserve(respos.size());
		for (unsigned i = 0, n = respos.size(); i < n; i++)
		{
			const char * addr = objects.begin() + respos[i].pos;
			res.add(addr, respos[i].val);
		}
	}

	if (verbose)
	{
		cout << "Scanned " << respos.size() << " objects out of " << total
				<< " in " << t.elapsed() << " s (" << objload.elapsed()
				<< " objload)\n";
	}

}

void GSSTreeSearcher::findTweeners(const GSSInternalNode* node,
		const MappableOctTree* min, const MappableOctTree* max, const MappableOctTree* orig,
		vector<result_info>& respos, unsigned level, bool computeDist)
{
	nodesVisited++;
	if (levelCnts.size() <= level)
	{
		levelCnts.resize(level + 1, 0);
		usefulLevelCnts.resize(level + 1, 0);
		maxlevelCnts.resize(level + 2, 0);
	}
	levelCnts[level]++;

	vector<const GSSInternalNode::Child *> goodchildren;
	unsigned n = node->size();
	for (unsigned i = 0; i < n; i++)
	{
		const GSSInternalNode::Child *child = node->getChild(i);

		if (fitsInbetween(child->getMIV(), child->getMSV(), min, max))
		{
			goodchildren.push_back(child);
		}
	}

	maxlevelCnts[level + 1] += n;
	unsigned oldres = respos.size();
	for (unsigned i = 0, nc = goodchildren.size(); i < nc; i++)
	{
		const GSSInternalNode::Child *child = goodchildren[i];
		if (child->isLeafPosition())
		{
			const GSSLeaf* next = (const GSSLeaf*) (leaves.begin()
					+ child->position());
			findTweeners(next, min, max, orig, respos, computeDist);
		}
		else
		{
			const GSSInternalNode* next =
					(const GSSInternalNode*) (internalNodes.begin()
							+ child->position());
			findTweeners(next, min, max, orig, respos, level + 1, computeDist);
		}
	}
	if (respos.size() > oldres)
		usefulLevelCnts[level]++;

}

//identify and trees in this leaf that fit
//if computeDist is true, compute a goodness of fit for reach result
void GSSTreeSearcher::findTweeners(const GSSLeaf* node,
		const MappableOctTree* min, const MappableOctTree* max, const MappableOctTree* orig,
		vector<result_info>& respos, bool computeDist)
{
	leavesVisited++;
	unsigned cnt = 0;
	for (unsigned i = 0, n = node->size(); i < n; i++)
	{
		const GSSLeaf::Child *child = node->getChild(i);

		if (fitsInbetween(&child->tree, &child->tree, min, max))
		{
			double goodness = 0;
			if (computeDist)
			{
				if(orig)
					goodness = volumeDist(orig,&child->tree);
				else
					goodness = shapeDistance(&child->tree, &child->tree, min, max);
			}
			respos.push_back(result_info(child->object_pos, goodness));
			cnt++;
		}
	}
	if (cnt == node->size())
		fullLeaves++;
}

