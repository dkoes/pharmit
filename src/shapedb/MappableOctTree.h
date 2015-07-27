/*
 * MappableOctTree.h
 *
 *  Created on: Oct 13, 2011
 *      Author: dkoes
 *
 *  Based off of the array oct tree (several variants were tried) this
 *  is a mappable version (no virtual functions, no pointers)
 */

#ifndef MAPPABLEOCTTREE_H_
#define MAPPABLEOCTTREE_H_

#include <vector>
#include <iostream>
#include <fstream>

#include "Cube.h"
#include "MGrid.h"
#include <boost/tuple/tuple.hpp>
using namespace std;

struct Vertex
{
	float x, y, z;

	Vertex() :
			x(0), y(0), z(0)
	{
	}
	Vertex(float X, float Y, float Z) :
			x(X), y(Y), z(Z)
	{
	}

	bool operator<(const Vertex& rhs) const
	{
		if (x != rhs.x)
			return x < rhs.x;
		if (y != rhs.y)
			return y < rhs.y;
		return z < rhs.z;
	}

	bool operator==(const Vertex& rhs) const
	{
		return x == rhs.x && y == rhs.y && z == rhs.z;
	}

	float distanceSq(const Vertex& rhs) const
	{
		float X = x - rhs.x;
		float Y = y - rhs.y;
		float Z = z - rhs.z;

		return X * X + Y * Y + Z * Z;
	}
};

//this can be 15 or 31; 31 is needs for storing large, detailed objects,
//15 is sufficient for molecular shape matching
#define MINDEX_BITS 15
struct MOctNode;
struct MChildNode
{
	union
	{
		//you would think there would be a more portably way to get this bit-packed alignment,
		//but I couldn't figure one out - isLeafs should overlap
		bool isLeaf :1;
		struct
		{
			bool isLeaf :1;
			unsigned pattern :8;
			unsigned numbits :4;
		}__attribute__((__packed__)) leaf;

		struct
		{
			bool isLeaf :1;
			unsigned index :MINDEX_BITS;
		}__attribute__((__packed__)) node;
	}__attribute__((__packed__));

	MChildNode() :
			isLeaf(true)
	{
		node.index = 0;
	}

	//use to create node
	MChildNode(unsigned i) :
			isLeaf(false)
	{
		node.index = i;
	}

	//use to create leaf
	MChildNode(unsigned pat, unsigned i) :
			isLeaf(true)
	{
		leaf.pattern = pat;
		leaf.numbits = i;
	}
	void intersectUnionVolume(const MOctNode* tree, const MChildNode& rhs,
			const MOctNode* rtree, float cubeVol, float& intersectval,
			float& unionval) const;

	float
	intersectVolume(const MOctNode* tree, const MChildNode& rhs,
			const MOctNode* rtree) const;
	float
	unionVolume(const MOctNode* tree, const MChildNode& rhs,
			const MOctNode* rtree) const;

	void invert(MOctNode* tree, float maxvol);

	bool containedIn(const MOctNode* tree, const MChildNode& rhs,
			const MOctNode* rtree) const;

	float volume(const MOctNode* tree, float dim3) const; //pass volume of full node
	float volume(const vector<MOctNode>& tree, float dim3) const;

	bool equals(const MOctNode* tree, const MChildNode& rhs,
			const MOctNode* rtree) const;

	void collectVertices(vector<Vertex>& vertices, Cube box,
			const MOctNode* tree) const;

	bool checkCoord(const MOctNode* tree, unsigned i, unsigned j, unsigned k,
			unsigned max) const;
	void countLeavesAtDepths(const MOctNode* tree, unsigned depth,
			vector<unsigned>& counts) const;
}__attribute__((__packed__));

struct MOctNode
{
	float vol; //cache the volume of this subtree to speed things up
	MChildNode children[8];

	MOctNode() :
			vol(0)
	{
	}

}__attribute__((__packed__));

class MappableOctTree
{
	MChildNode root;
	float dimension;
	unsigned N; //number of octnodes
	MOctNode tree[]; //must memalloc beyond

	MappableOctTree()
	{
		//creation must be performed externally to properly allocate memory
	}

	//assume memory is allocate to tree
	MappableOctTree(float dim, const MChildNode& r,
			const vector<MOctNode>& nodes) :
			root(r), dimension(dim), N(nodes.size())
	{
		for (unsigned i = 0, n = nodes.size(); i < n; i++)
		{
			tree[i] = nodes[i];
		}
	}

	MappableOctTree(const MappableOctTree& rhs);

	static MChildNode createFrom_r(unsigned N, MChildNode* nodes,
			const MappableOctTree** trees, vector<MOctNode>& newtree,
			bool isUnion, float cubeVol);
	MChildNode createTruncated_r(const MChildNode& node, float res, bool fillIn,
			vector<MOctNode>& newtree, float cubeVol) const;
	MChildNode createRounded_r(const MChildNode& node, float res, bool fillIn,
			vector<MOctNode>& newtree, float cubeVol) const;

	static bool createRoundedSet_r(unsigned N, MChildNode* nodes,
			const MappableOctTree** trees, bool roundUp,
			vector<vector<MOctNode> >& newtrees, vector<MChildNode>& newroots,
			float cubeVol);

	static MappableOctTree* newFromVector(const vector<MOctNode>& newtree,
			const MChildNode& newroot, float dim);
public:

	MappableOctTree* clone() const;

	void invert();

	//return a tree that is rounded to resolution res, if fillIn true,
	//round up, otherwise down
	MappableOctTree* createTruncated(float res, bool fillIn) const;
	//return a tree where nodes are rounded up/down to leaves if they contain/exclude
	// <= the specified volume
	MappableOctTree* createRounded(float vol, bool fillIn) const;

	//return intersect of the n trees found in arr
	static MappableOctTree* createFromIntersection(unsigned N,
			const MappableOctTree** in);
	//union
	static MappableOctTree* createFromUnion(unsigned N,
			const MappableOctTree** in);

	//round trees in in as much as possible while maintaining local distiguishability
	//and put the result in out
	static bool createRoundedSet(unsigned N, const MappableOctTree**in,
			bool roundUp, MappableOctTree** out);

	//volume calculations that don't require creating a tmp tree
	float intersectVolume(const MappableOctTree * rhs) const;
	float unionVolume(const MappableOctTree *rhs) const;
	void intersectUnionVolume(const MappableOctTree *rhs, float& ival,
			float& uval) const;

	bool containedIn(const MappableOctTree *rhs) const;

	//return total volume contained in octtree
	float volume() const;

	//return number of leaves
	unsigned leaves() const;

	void write(ostream& out) const;

	float hausdorffDistance(const MappableOctTree* B) const;
	float relativeVolumeDistance(const MappableOctTree * B) const;
	float absoluteVolumeDistance(const MappableOctTree * B) const;

	unsigned bytes() const
	{
		return sizeof(MappableOctTree) + N * sizeof(MOctNode);
	}

	unsigned nodes() const
	{
		return N;
	}

	bool equals(const MappableOctTree* rhs) const;

	//create a grid representing this octtree
	void makeGrid(MGrid& grid, float res) const;

	void dumpGrid(ostream& out, float res) const;
	void dumpRawGrid(ostream& out, float res) const;
	void dumpMiraGrid(ostream& out, float res) const;
	void dumpAD4Grid(ostream& out, float res) const;
	void dumpSproxelGrid(ostream& out, float res, const string& voxelValue =
			"#FFFFFFFF") const;
	void countLeavesAtDepths(vector<unsigned>& counts) const;
private:
	template<class Object>
	static MChildNode create_r(float res, const Cube& cube, const Object& obj,
			vector<MOctNode>& tree)
	{
		MChildNode ret;
		if (cube.getDimension() <= res) //must be a leaf
		{
			ret.isLeaf = true;
			//ultimately, follow standard practice and only evaluate the
			//intersection with the center of the cube
			float x = 0, y = 0, z = 0;
			cube.getCenter(x, y, z);
			if (obj.containsPoint(x, y, z))
			{
				ret.leaf.pattern = 0xff;
				ret.leaf.numbits = 8;
			}
			else
			{
				ret.leaf.pattern = 0;
				ret.leaf.numbits = 0;
			}

			return ret;
		}

		//does the object overlap with this cube?
		bool intersects = obj.intersects(cube);
		if (!intersects)
		{
			//no overlap, all done
			ret.isLeaf = true;
			ret.leaf.pattern = 0;
			ret.leaf.numbits = 0;
		}
		else //subdivide into children
		{
			//assume this is an interior node
			ret.isLeaf = false;
			unsigned filledcnt = 0;
			unsigned pat = 0;
			unsigned bitcnt = 0;
			//allocate space
			unsigned pos = tree.size();
			ret.node.index = tree.size();

			if (tree.size() >= (1U << MINDEX_BITS))
			{
				cerr
						<< "Too many nodes for MINDEX_BITS. Must recompile to support larger octtress.\n";
				abort();
			}
			tree.push_back(MOctNode());
			tree[pos].vol = 0;
			for (unsigned i = 0; i < 8; i++)
			{
				Cube newc = cube.getOctant(i);
				MChildNode child = create_r(res, newc, obj, tree);
				tree[pos].children[i] = child;
				if (child.isLeaf)
				{
					if (child.leaf.pattern == 0)
						filledcnt++;
					else if (child.leaf.pattern == 0xff)
					{
						filledcnt++;
						pat |= (1 << i);
						bitcnt++;
					}
				}
				tree[pos].vol += child.volume(tree, newc.volume());
			}
			//are all the children full or empty? then can represent as a pattern
			if (filledcnt == 8)
			{
				tree.pop_back();
				ret.isLeaf = true;
				ret.leaf.pattern = pat;
				ret.leaf.numbits = bitcnt;
			}
		}
		return ret;
	}

public:

	//return a mappable oct tree, the pointer should be deallocated using free
	template<class Object>
	static MappableOctTree* create(float dim, float res, const Object& obj)
	{
		vector<MOctNode> nodes;
		MChildNode root = create_r(res, Cube(dim, -dim / 2, -dim / 2, -dim / 2),
				obj, nodes);
		unsigned N = nodes.size();
		void *mem = malloc(sizeof(MappableOctTree) + N * sizeof(MOctNode));

		MappableOctTree *ret = new (mem) MappableOctTree(dim, root, nodes);

		return ret;
	}

	static MappableOctTree* createFromGrid(const MGrid& grid)
	{
		return create(grid.getDimension(), grid.getResolution(),grid);
	}
}__attribute__((__packed__));



#endif /* MAPPABLEOCTTREE_H_ */
