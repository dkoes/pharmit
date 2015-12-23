/*
Pharmit
Copyright (c) David Ryan Koes, University of Pittsburgh and contributors.
All rights reserved.

Pharmit is licensed under both the BSD 3-clause license and the GNU
Public License version 2. Any use of the code that retains its reliance
on the GPL-licensed OpenBabel library is subject to the terms of the GPL2.

Use of the Pharmit code independently of OpenBabel (or any other
GPL2 licensed software) may choose between the BSD or GPL licenses.

See the LICENSE file provided with the distribution for more information.

*/
/*
 * MappableOctTree.cpp
 *
 *  Created on: Oct 13, 2011
 *      Author: dkoes
 */

#include "MappableOctTree.h"
#include "MiraObject.h"
#include <cstring>
#include <cassert>
#include <boost/multi_array.hpp>

using namespace boost;

MappableOctTree* MappableOctTree::clone() const
{
	unsigned sz = bytes();
	void* mem = malloc(sz);
	memcpy(mem, this, sz);
	return (MappableOctTree*) mem;
}

//return volume of this node given that the full volume of the cube is dim3
//to save space while keeping most performance, volume is cached in octnodes
float MChildNode::volume(const MOctNode* tree, float dim3) const
		{
	if (isLeaf)
	{
		return leaf.numbits * dim3 / 8;
	}
	else
	{
		return tree[node.index].vol;
	}
}

//volume version for tree under construction
float MChildNode::volume(const vector<MOctNode>& tree, float dim3) const
		{
	if (isLeaf)
	{
		return leaf.numbits * dim3 / 8;
	}
	else
	{
		return tree[node.index].vol;
	}
}

MChildNode MappableOctTree::createFrom_r(unsigned N, MChildNode* nodes,
		const MappableOctTree** trees, vector<MOctNode>& newtree, bool isUnion,
		float cubeVol)
{
	float bitvol = cubeVol / 8;
	unsigned numFilled = 0;
	unsigned andpat = 0xff, orpat = 0;
	for (unsigned i = 0; i < N; i++)
	{
		if (nodes[i].isLeaf)
		{
			numFilled++;
			andpat &= nodes[i].leaf.pattern;
			orpat |= nodes[i].leaf.pattern;
		}
	}

	if (isUnion && (numFilled == N || orpat == 0xff)) //just return a patterned node
	{
		unsigned nb = __builtin_popcount(orpat);
		MChildNode ret(orpat, nb);
		return ret;
	}
	else if (!isUnion && (numFilled == N || andpat == 0))
	{
		unsigned nb = __builtin_popcount(andpat);
		MChildNode ret(andpat, nb);
		return ret;
	}
	else //need to look at non-leaf children
	{
		unsigned pos = newtree.size();
		MChildNode ret(pos);
		newtree.push_back(MOctNode());

		unsigned numFullEmpty = 0;
		unsigned newPattern = 0;
		MChildNode nextnodes[N + 1];
		const MappableOctTree *nexttrees[N + 1];

		//setup nonleaf trees
		unsigned n = 0;
		for (unsigned i = 0; i < N; i++)
		{
			if (!nodes[i].isLeaf)
			{
				nexttrees[n] = trees[i];
				n++;
			}
		}
		if (numFilled > 0)
			nexttrees[n++] = NULL; //dummy for merge of leaves
		for (unsigned i = 0; i < 8; i++)
		{
			//create array of sub-octants from non-leafs
			unsigned cnt = 0;
			for (unsigned j = 0; j < N; j++)
			{
				if (!nodes[j].isLeaf)
				{
					nextnodes[cnt] =
							trees[j]->tree[nodes[j].node.index].children[i];
					cnt++;
				}
			}
			//summerize filled nodes
			if (numFilled > 0)
			{
				unsigned pat = isUnion ? orpat : andpat;
				if (pat & (1 << i))
					pat = 0xff;
				else
					pat = 0;
				nextnodes[cnt] = MChildNode(pat, pat ? 8 : 0);
				cnt++;
			}
			assert(cnt == n);

			MChildNode nchild = createFrom_r(n, nextnodes, nexttrees, newtree,
					isUnion, bitvol);
			newtree[pos].children[i] = nchild;

			if (nchild.isLeaf
					&& (nchild.leaf.pattern == 0 || nchild.leaf.pattern == 0xff))
			{
				numFullEmpty++;
				if (nchild.leaf.pattern == 0xff)
					newPattern |= 1 << i;
			}

			newtree[pos].vol += nchild.volume(newtree, bitvol);
		}

		if (numFullEmpty == 8)
		{
			newtree.pop_back();
			ret.isLeaf = true;
			ret.leaf.pattern = newPattern;
			ret.leaf.numbits = __builtin_popcount(newPattern);
			//volume is correct
		}
		return ret;
	}
}

//invert tree, this can be done inplace
void MappableOctTree::invert()
{
	float cubeVol = dimension * dimension * dimension;
	float expectedV = cubeVol - root.volume(tree, cubeVol);
	root.invert(tree, cubeVol);
	assert(root.volume(tree, cubeVol) == expectedV);
}

MappableOctTree* MappableOctTree::newFromVector(const vector<MOctNode>& newtree,
		const MChildNode& newroot, float dim)
{
	unsigned sz = sizeof(MappableOctTree) + newtree.size() * sizeof(MOctNode);
	void *mem = malloc(sz);
	return new (mem) MappableOctTree(dim, newroot, newtree);
}

//union
MappableOctTree* MappableOctTree::createFromUnion(unsigned N,
		const MappableOctTree** in)
{
	MChildNode roots[N];
	float dim = 0;
	for (unsigned i = 0; i < N; i++)
	{
		roots[i] = in[i]->root;
		if (i == 0)
			dim = in[i]->dimension;
		else
			assert(dim == in[i]->dimension);
	}
	vector<MOctNode> newtree;
	MChildNode newroot = createFrom_r(N, roots, in, newtree, true,
			dim * dim * dim);

	return newFromVector(newtree, newroot, dim);
}

//mutually recurse the N trees, round up/down nodes iff doing so will not
//change the local equality relationship between the N trees
bool MappableOctTree::createRoundedSet_r(unsigned N, MChildNode* nodes,
		const MappableOctTree** trees, bool roundUp,
		vector<vector<MOctNode> >& newtrees, vector<MChildNode>& newroots,
		float cubeVol)
{
	unsigned numFilled = 0;
	unsigned numFull = 0;
	unsigned numEmpty = 0;
	int lastPat = -1;
	bool differentPats = false;
	float bitvol = cubeVol / 8;
	for (unsigned i = 0; i < N; i++)
	{
		if (nodes[i].isLeaf)
		{
			numFilled++;

			if (nodes[i].leaf.pattern == 0xff)
				numFull++;
			else if (nodes[i].leaf.pattern == 0)
				numEmpty++;
			else if (lastPat == -1)
			{
				lastPat = nodes[i].leaf.pattern;
			}
			else if (nodes[i].leaf.pattern != lastPat)
				differentPats = true;
		}
	}

	if (numFilled == N) //all leaves
	{
		if (roundUp)
		{
			//can only round up equivalent leaves if remaining leaves are empty
			if (!differentPats && lastPat > 0 && numFull == 0)
			{
				for (unsigned i = 0; i < N; i++)
				{
					if (nodes[i].leaf.pattern == 0)
						newroots[i] = nodes[i];
					else
						newroots[i] = MChildNode(0xff, 8);
				}
				return true;
			}
			else //just copy
			{
				for (unsigned i = 0; i < N; i++)
				{
					newroots[i] = nodes[i];
				}
				return false;
			}
		}
		else //round down
		{
			//leaves must be equivalent or full
			if (!differentPats && lastPat > 0 && numEmpty == 0)
			{
				for (unsigned i = 0; i < N; i++)
				{
					if (nodes[i].leaf.pattern == 0)
						newroots[i] = nodes[i];
					else
						newroots[i] = MChildNode(0, 0);
				}
				return true;
			}
			else //just duplicate
			{
				for (unsigned i = 0; i < N; i++)
				{
					newroots[i] = nodes[i];
				}
				return false;
			}
		}
	}
	else //must descend non-leaf tree parts
	{
		//create new octnodes for all non-leaf trees
		unsigned positions[N];
		for (unsigned i = 0; i < N; i++)
		{
			if (nodes[i].isLeaf)
			{
				//stays a leaf
				newroots[i] = nodes[i];
			}
			else
			{
				unsigned pos = newtrees[i].size();
				positions[i] = pos;
				newroots[i] = MChildNode(pos);
				newtrees[i].push_back(MOctNode());
			}
		}

		MChildNode subnodes[N];
		bool changed = false;
		unsigned leafCnts[N];
		unsigned leafPatterns[N];
		memset(leafCnts, 0, sizeof(leafCnts));
		memset(leafPatterns, 0, sizeof(leafPatterns));

		for (unsigned o = 0; o < 8; o++)
		{
			//setup subnodes, divide leaf nodes
			for (unsigned i = 0; i < N; i++)
			{
				if (nodes[i].isLeaf)
				{
					if (nodes[i].leaf.pattern & (1 << o))
						subnodes[i] = MChildNode(0xff, 8);
					else
						subnodes[i] = MChildNode(0, 0);
				}
				else
				{
					subnodes[i] =
							trees[i]->tree[nodes[i].node.index].children[o];
				}
			}

			vector<MChildNode> newnodes(N);
			changed |= createRoundedSet_r(N, subnodes, trees, roundUp, newtrees,
					newnodes, bitvol);

			for (unsigned i = 0; i < N; i++)
			{
				if (!nodes[i].isLeaf)
				{
					newtrees[i][positions[i]].children[o] = newnodes[i];
					if (newnodes[i].isLeaf)
					{
						if (newnodes[i].leaf.pattern == 0xff)
						{
							leafCnts[i]++;
							leafPatterns[i] |= (1 << o);
						}
						else if (newnodes[i].leaf.pattern == 0)
							leafCnts[i]++;
					}
					newtrees[i][positions[i]].vol += newnodes[i].volume(
							newtrees[i], bitvol);
				}
			}
		}

		//coalesce any all leaves
		for (unsigned i = 0; i < N; i++)
		{
			if (leafCnts[i] == 8)
			{
				newtrees[i].pop_back();
				newroots[i].isLeaf = true;
				newroots[i].leaf.pattern = leafPatterns[i];
				newroots[i].leaf.numbits = __builtin_popcount(leafPatterns[i]);
				//volume is correct
			}
		}

		return changed;
	}
}

//round trees in in as much as possible while maintaining local distiguishability
//and put the result in out
bool MappableOctTree::createRoundedSet(unsigned N, const MappableOctTree**in,
		bool roundUp, MappableOctTree** out)
{
	MChildNode roots[N];
	float dim = 0;
	for (unsigned i = 0; i < N; i++)
	{
		roots[i] = in[i]->root;
		if (i == 0)
			dim = in[i]->dimension;
		else
			assert(dim == in[i]->dimension);
	}
	vector<vector<MOctNode> > newtrees(N);
	vector<MChildNode> newroots(N);
	bool changed = createRoundedSet_r(N, roots, in, roundUp, newtrees, newroots,
			dim * dim * dim);

	for (unsigned i = 0; i < N; i++)
	{
		out[i] = newFromVector(newtrees[i], newroots[i], dim);
	}

	return changed;
}

MappableOctTree* MappableOctTree::createFromIntersection(unsigned N,
		const MappableOctTree** in)
{
	MChildNode roots[N];
	float dim = 0;
	for (unsigned i = 0; i < N; i++)
	{
		roots[i] = in[i]->root;
		if (i == 0)
			dim = in[i]->dimension;
		else
			assert(dim == in[i]->dimension);
	}
	vector<MOctNode> newtree;
	MChildNode newroot = createFrom_r(N, roots, in, newtree, false,
			dim * dim * dim);

	return newFromVector(newtree, newroot, dim);
}

MChildNode MappableOctTree::createTruncated_r(const MChildNode& node, float res,
		bool fillIn, vector<MOctNode>& newtree, float cubeVol) const
		{
	if (cubeVol <= res * res * res) //at level that needs to be truncated to a leaf
	{
		if (fillIn && node.volume(tree, cubeVol) != 0)
		{
			return MChildNode(0xff, 8);
		}
		else
		{
			return MChildNode(0, 0);
		}
	}
	else if (node.isLeaf)
	{
		return node;
	}
	else //recursively descend
	{
		unsigned pos = newtree.size();
		MChildNode ret(pos);
		newtree.push_back(MOctNode());

		unsigned numFullEmpty = 0;
		unsigned newPattern = 0;
		for (unsigned i = 0; i < 8; i++)
		{
			MChildNode nchild = createTruncated_r(
					tree[node.node.index].children[i], res, fillIn, newtree,
					cubeVol / 8);
			newtree[pos].children[i] = nchild;
			if (nchild.isLeaf
					&& (nchild.leaf.pattern == 0 || nchild.leaf.pattern == 0xff))
			{
				numFullEmpty++;
				if (nchild.leaf.pattern == 0xff)
					newPattern |= 1 << i;
			}

			newtree[pos].vol += nchild.volume(newtree, cubeVol / 8);
		}

		if (numFullEmpty == 8)
		{
			newtree.pop_back();
			ret.isLeaf = true;
			ret.leaf.pattern = newPattern;
			ret.leaf.numbits = __builtin_popcount(newPattern);
		}
		return ret;
	}
}

MappableOctTree* MappableOctTree::createTruncated(float res, bool fillIn) const
		{
	vector<MOctNode> newtree;
	MChildNode newroot = createTruncated_r(root, res, fillIn, newtree,
			dimension * dimension * dimension);

	return newFromVector(newtree, newroot, dimension);
}

MChildNode MappableOctTree::createRounded_r(const MChildNode& node,
		float volcutoff, bool fillIn, vector<MOctNode>& newtree,
		float maxvol) const
		{
	float vol = node.volume(tree, maxvol);
	float empty = maxvol - vol;

	if (fillIn && empty < volcutoff)
	{
		return MChildNode(0xff, 8);
	}
	else if (!fillIn && vol < volcutoff)
	{
		return MChildNode(0, 0);
	}
	else if (node.isLeaf)
	{
		return node;
	}
	else //recursively descend
	{
		unsigned pos = newtree.size();
		MChildNode ret(pos);
		newtree.push_back(MOctNode());

		unsigned numFullEmpty = 0;
		unsigned newPattern = 0;
		for (unsigned i = 0; i < 8; i++)
		{
			MChildNode nchild = createRounded_r(
					tree[node.node.index].children[i], volcutoff, fillIn,
					newtree, maxvol / 8);
			newtree[pos].children[i] = nchild;
			if (nchild.isLeaf
					&& (nchild.leaf.pattern == 0 || nchild.leaf.pattern == 0xff))
			{
				numFullEmpty++;
				if (nchild.leaf.pattern == 0xff)
					newPattern |= 1 << i;
			}

			newtree[pos].vol += nchild.volume(newtree, maxvol / 8);
		}

		if (numFullEmpty == 8)
		{
			newtree.pop_back();
			ret.isLeaf = true;
			ret.leaf.pattern = newPattern;
			ret.leaf.numbits = __builtin_popcount(newPattern);
			//volume is correct
		}
		return ret;
	}
}

MappableOctTree* MappableOctTree::createRounded(float vol, bool fillIn) const
		{
	vector<MOctNode> newtree;
	MChildNode newroot = createRounded_r(root, vol, fillIn, newtree,
			dimension * dimension * dimension);

	return newFromVector(newtree, newroot, dimension);
}

//volume calculations that don't require creating a tmp tree
float MappableOctTree::intersectVolume(const MappableOctTree * rhs) const
		{
	float ival = 0, uval = 0;
	assert(dimension == rhs->dimension);
	root.intersectUnionVolume(tree, rhs->root, rhs->tree,
			dimension * dimension * dimension, ival, uval);
	return ival;
}

float MappableOctTree::unionVolume(const MappableOctTree *rhs) const
		{
	float ival = 0, uval = 0;
	assert(dimension == rhs->dimension);
	root.intersectUnionVolume(tree, rhs->root, rhs->tree,
			dimension * dimension * dimension, ival, uval);
	return uval;
}

void MappableOctTree::intersectUnionVolume(const MappableOctTree *rhs,
		float& ival, float& uval) const
		{
	ival = 0;
	uval = 0;
	assert(dimension == rhs->dimension);

	root.intersectUnionVolume(tree, rhs->root, rhs->tree,
			dimension * dimension * dimension, ival, uval);
}

bool MappableOctTree::containedIn(const MappableOctTree *rhs) const
		{
	return root.containedIn(tree, rhs->root, rhs->tree);
}

//return total volume contained in octtree
float MappableOctTree::volume() const
{
	return root.volume(tree, dimension * dimension * dimension);
}

//return number of leaves
unsigned MappableOctTree::leaves() const
{
	unsigned cnt = 0;
	for (unsigned i = 0; i < N; i++)
	{
		for (unsigned j = 0; j < 8; j++)
		{
			if (tree[i].children[j].isLeaf)
				cnt++;
		}
	}
	return cnt;
}

void MappableOctTree::write(ostream& out) const
		{
	unsigned sz = sizeof(MappableOctTree) + N * sizeof(MOctNode);
	out.write((char*) this, sz);
}

float MappableOctTree::relativeVolumeDistance(const MappableOctTree * rhs) const
{
	float ival = 0, uval = 0;
	assert(dimension == rhs->dimension);

	root.intersectUnionVolume(tree, rhs->root, rhs->tree,
			dimension * dimension * dimension, ival, uval);
	if (uval == 0) //possible when don't have anchored shapes
		return 0.0; //two empty shapes..
	return 1 - ival / uval;
}

float MappableOctTree::absoluteVolumeDistance(const MappableOctTree * rhs) const
		{
	float ival = 0, uval = 0;
	assert(dimension == rhs->dimension);

	root.intersectUnionVolume(tree, rhs->root, rhs->tree,
			dimension * dimension * dimension, ival, uval);

	return uval - ival;
}

void MChildNode::invert(MOctNode* tree, float maxvol)
{
	if (isLeaf)
	{
		leaf.pattern = ~leaf.pattern;
		leaf.numbits = 8 - leaf.numbits;
	}
	else
	{
		tree[node.index].vol = 0;
		for (unsigned i = 0; i < 8; i++)
		{
			tree[node.index].children[i].invert(tree, maxvol / 8);
			tree[node.index].vol += tree[node.index].children[i].volume(tree,
					maxvol / 8);
		}
	}
}

//compute both the intersection and union volume at once
void MChildNode::intersectUnionVolume(const MOctNode* tree,
		const MChildNode& rhs, const MOctNode* rtree, float cubeVol,
		float& intersectval, float& unionval) const
		{
	float bitv = cubeVol / 8;
	if (rhs.isLeaf && isLeaf)
	{
		if (leaf.pattern == 0 || rhs.leaf.pattern == 0)
			intersectval += 0; //no intersect val
		else
		{
			unsigned intpat = leaf.pattern & rhs.leaf.pattern;
			intersectval += __builtin_popcount(intpat) * bitv;
		}
		if (leaf.pattern == 0)
			unionval += rhs.volume(rtree, cubeVol);
		else if (rhs.leaf.pattern == 0)
			unionval += volume(tree, cubeVol);
		else
		{
			unsigned orpat = leaf.pattern | rhs.leaf.pattern;
			unionval += __builtin_popcount(orpat) * bitv;
		}
	}
	else if (isLeaf && leaf.pattern == 0xff)
	{
		unionval += volume(tree, cubeVol);
		intersectval += rhs.volume(rtree, cubeVol);
	}
	else if (rhs.isLeaf && rhs.leaf.pattern == 0xff)
	{
		unionval += rhs.volume(rtree, cubeVol);
		intersectval += volume(tree, cubeVol);
	}
	else if (rhs.isLeaf && rhs.leaf.pattern == 0)
	{
		//vol of this tree
		unionval += volume(tree, cubeVol);
	}
	else if (isLeaf && leaf.pattern == 0)
	{
		// from rhs
		unionval += rhs.volume(rtree, cubeVol);
	}
	else if (isLeaf)
	{
		assert(!rhs.isLeaf);
		assert(leaf.numbits > 0);
		//rhs has children, have to compare bit by bit
		for (unsigned i = 0; i < 8; i++)
		{
			MChildNode rchild = rtree[rhs.node.index].children[i];
			if (leaf.pattern & (1 << i))
			{
				unionval += bitv;
				intersectval += rchild.volume(rtree, bitv);
			}
			else
			{
				unionval += rchild.volume(rtree, bitv);
			}
		}
	}
	else if (rhs.isLeaf)
	{
		assert(!isLeaf);
		assert(rhs.leaf.numbits > 0);
		//rhs has children, have to compare bit by bit
		for (unsigned i = 0; i < 8; i++)
		{
			MChildNode child = tree[node.index].children[i];
			if (rhs.leaf.pattern & (1 << i))
			{
				unionval += bitv;
				intersectval += child.volume(tree, bitv);
			}
			else
			{
				unionval += child.volume(tree, bitv);
			}
		}
	}
	else //both have children
	{
		for (unsigned i = 0; i < 8; i++)
		{
			tree[node.index].children[i].intersectUnionVolume(tree,
					rtree[rhs.node.index].children[i], rtree, bitv,
					intersectval, unionval);
		}
	}
}

//return true if this is contained in rhs - short circuit eval
bool MChildNode::containedIn(const MOctNode* tree, const MChildNode& rhs,
		const MOctNode* rtree) const
		{
	if (rhs.isLeaf && isLeaf)
	{
		return (leaf.pattern & rhs.leaf.pattern) == leaf.pattern;
	}
	else if (isLeaf && leaf.pattern == 0)
	{
		return true;
	}
	else if (rhs.isLeaf && rhs.leaf.pattern == 0)
	{
		return false; //somthing not in nothing (nothing handled above)
	}
	else if (rhs.isLeaf && rhs.leaf.pattern == 0xff)
	{
		return true;
	}
	else if (isLeaf && leaf.pattern == 0xff)
	{
		return false;
	}
	else if (isLeaf)
	{
		assert(!rhs.isLeaf);
		assert(leaf.numbits > 0);
		//rhs has children, have to compare bit by bit
		for (unsigned i = 0; i < 8; i++)
		{
			if (leaf.pattern & (1 << i))
			{
				MChildNode rchild = rtree[rhs.node.index].children[i];
				//rchild must be a leaf and full
				if (!rchild.isLeaf || rchild.leaf.pattern != 0xff)
					return false;
			}
		}
		return true;
	}
	else if (rhs.isLeaf)
	{
		assert(!isLeaf);
		assert(rhs.leaf.numbits > 0);
		//rhs has children, have to compare bit by bit
		for (unsigned i = 0; i < 8; i++)
		{
			if (!(rhs.leaf.pattern & (1 << i)))
			{
				MChildNode child = tree[node.index].children[i];
				//if rhs's bit is empty, then our child must be empty
				if (!child.isLeaf || child.leaf.pattern != 0)
					return false;
			}
		}
		return true;
	}
	else //both have children
	{
		for (unsigned i = 0; i < 8; i++)
		{
			if (!tree[node.index].children[i].containedIn(tree,
					rtree[rhs.node.index].children[i], rtree))
				return false;
		}
		return true;;
	}
	return false;
}

static void addVertices(vector<Vertex>& vertices, Cube box)
{
	for (unsigned i = 0; i < 8; i++)
	{
		float x = box.x;
		float y = box.y;
		float z = box.z;

		if (i & 1)
			x += box.dim;
		if (i & 2)
			y += box.dim;
		if (i & 4)
			z += box.dim;

		vertices.push_back(Vertex(x, y, z));
	}
}

//append vertices of this child to vertices, this child's box is provided
void MChildNode::collectVertices(vector<Vertex>& vertices, Cube box,
		const MOctNode* tree) const
		{
	if (isLeaf)
	{
		if (leaf.pattern == 0)
			return; //empty, no vertices
		else if (leaf.pattern == 0xff) //full, add this box's vertices
		{
			addVertices(vertices, box);
		}
		else
		{
			for (unsigned i = 0; i < 8; i++)
			{
				if (leaf.pattern & (1 << i))
				{
					addVertices(vertices, box.getOctant(i));
				}
			}
		}
	}
	else
	{
		for (unsigned i = 0; i < 8; i++)
		{
			tree[node.index].children[i].collectVertices(vertices,
					box.getOctant(i), tree);
		}
	}
}

//sort, then remove duplicates and remove vertices that appear 8 times altogether
static void uniqueRemove8ify(vector<Vertex>& V)
{
	sort(V.begin(), V.end());

	unsigned i = 1;
	unsigned n = V.size();
	unsigned cnt = 1;
	unsigned last = 0;
	while (i < n)
	{
		if (V[i] == V[last])
		{
			cnt++;
			i++;
		}
		else
		{
			if (cnt < 8)
				last++; //keep what we copied
			V[last] = V[i];
			cnt = 1;
		}
	}

	if (cnt < 8)
		V.resize(last + 1);
	else
		V.resize(last);
}

//the Hausdorff distance is the maximum minimum distance between the
//two shapes; it is not symmetric, but this computes the symmetric
//(max H(A,B),H(B,A)) distance
//
//TODO: make this efficient if it turns out to be worthwhile, right now
//do the easiest implimentation: collect all leaf corner coordinates, sort,
//remove completely buried (8 versions), compute all pairs distances
float MappableOctTree::hausdorffDistance(const MappableOctTree* B) const
{
	vector<Vertex> Av, Bv;
	float dim = dimension;

	Cube box(dim, -dim / 2, -dim / 2, -dim / 2);
	root.collectVertices(Av, box, tree);
	B->root.collectVertices(Bv, box, B->tree);

	uniqueRemove8ify(Av);
	uniqueRemove8ify(Bv);

	//compute distances once
	vector<float> Bmin(Bv.size(), HUGE_VAL);
	float Amax = 0;
	for (unsigned i = 0, n = Av.size(); i < n; i++)
	{
		float Amin = HUGE_VAL;
		for (unsigned j = 0, m = Bv.size(); j < m; j++)
		{
			float dist = Av[i].distanceSq(Bv[j]);
			if (dist < Amin)
				Amin = dist;

			if (dist < Bmin[j])
				Bmin[j] = dist;
		}
		if (Amin > Amax)
			Amax = Amin;
	}

	float Bmax = *std::max_element(Bmin.begin(), Bmin.end());
//cout << Av.size() << " " << oldA << " " << Bv.size() << " " << oldB << " " << Amax << " " << Bmax << "\n";
	return sqrt(std::max(Amax, Bmax));
}

bool MChildNode::equals(const MOctNode* tree, const MChildNode& rhs,
		const MOctNode* rtree) const
		{
	if (isLeaf)
	{
		if (rhs.isLeaf)
		{
			if (rhs.leaf.pattern == leaf.pattern)
				return true;
			else
				return false;
		}
		else
			return false;
	}
	else if (rhs.isLeaf)
	{
		return false;
	}
	else
	{
		for (unsigned i = 0; i < 8; i++)
		{
			if (!tree[node.index].children[i].equals(tree,
					rtree[rhs.node.index].children[i], rtree))
				return false;
		}
		return true;
	}
	return true;
}

bool MappableOctTree::equals(const MappableOctTree* rhs) const
		{
	return root.equals(tree, rhs->root, rhs->tree);
}

//recurse down and see if i,j,k is set
//max is the integral dimension of this node
bool MChildNode::checkCoord(const MOctNode *tree, unsigned i, unsigned j,
		unsigned k, unsigned max) const
		{
	//first compute the octant identified by i,j,k
	unsigned half = max / 2;
	unsigned oct = 0;
	if (i >= half)
	{
		oct |= 1;
		i -= half;
	}
	if (j >= half)
	{
		oct |= 2;
		j -= half;
	}
	if (k >= half)
	{
		oct |= 4;
		k -= half;
	}

	if (isLeaf)
	{
		return leaf.pattern & (1 << oct);
	}
	else
	{
		return tree[node.index].children[oct].checkCoord(tree, i, j, k, half);
	}
}

//dump an mmp formated grid
void MappableOctTree::dumpGrid(ostream& out, float res) const
		{
	//size
	float dim = dimension;
	out << "header\n";
	unsigned max = ceil(dim / res);
	out << max << " " << max << " " << max << "\n";
	//origin
	out << -dim / 2 << " " << -dim / 2 << " " << -dim / 2 << "\n";
	//resolution
	out << res << " " << res << " " << res << "\n";

	//now coordinates
	//definitely not the most efficient way
	for (unsigned i = 0; i < max; i++)
	{
		for (unsigned j = 0; j < max; j++)
		{
			for (unsigned k = 0; k < max; k++)
			{
				if (root.checkCoord(tree, i, j, k, max))
					out << "1 ";
				else
					out << "0 ";
			}
			out << "\n";
		}
		out << "\n";
	}
}

//no header, just binary output
void MappableOctTree::dumpRawGrid(ostream& out, float res) const
		{
	char one = 0xff;
	char zero = 0;
	unsigned max = ceil(dimension / res);

	//now coordinates
	//definitely not the most efficient way
	for (unsigned i = 0; i < max; i++)
	{
		for (unsigned j = 0; j < max; j++)
		{
			for (unsigned k = 0; k < max; k++)
			{
				if (root.checkCoord(tree, i, j, k, max))
					out.write(&one, 1);
				else
					out.write(&zero, 1);
			}
		}
	}
}

//autodock 4 grid
void MappableOctTree::dumpAD4Grid(ostream& out, float res) const
{

	unsigned max = ceil(dimension / res);

	out << "GRID_PARAMETER_FILE\nGRID_DATA_FILE\nMACROMOLECULE\n";
	out << "SPACING " << res << "\n";
	out << "NELEMENTS " << max-1 << " " << max-1 << " " << max-1 << "\n";
	out << "CENTER 0 0 0\n";
	//now coordinates - z,y,x
	//definitely not the most efficient way
	for (unsigned k = 0; k < max; k++)
	{
		for (unsigned j = 0; j < max; j++)
		{
			for (unsigned i = 0; i < max; i++)
			{
				if (root.checkCoord(tree, i, j, k, max))
					out << "1\n";
				else
					out << "0\n";
			}
		}
	}
}

#include <arpa/inet.h>

//sproxel csv format, voxelValue is the string to specify voxels occupied by the object
void MappableOctTree::dumpSproxelGrid(ostream& out, float res,
		const string& voxelValue) const
		{
	unsigned short max = ceil(dimension / res);
	out << max << "," << max << "," << max << "\n";

	string blank = "#00000000";
	for (unsigned j = 0; j < max; j++)
	{ //sproxel csv is in xz slices
		for (unsigned i = 0; i < max; i++)
		{
			for (unsigned k = 0; k < max; k++)
			{
				if (k != 0)
					out << ",";
				if (root.checkCoord(tree, i, j, k, max))
					out << voxelValue;
				else
					out << blank;
			}
			out << "\n";
		}
		out << "\n";
	}
}

//mira grid format
void MappableOctTree::dumpMiraGrid(ostream& out, float res) const
		{
	unsigned short max = ceil(dimension / res);

	MiraHeader fileheader =
			{
					{ 'V', 'O', 'X', 'E', 'L' }, 26, 1, htons(max), htons(max), htons(
							max), 0, 256, (unsigned)(256 + 3 * max * sizeof(double)) };

	//we store resolution information in unused part of header
	stringstream str;
	str << "resolution = " << res;
	strncpy(fileheader.unused, str.str().c_str(), sizeof(fileheader.unused)-1);

	char one = 0xff;
	char zero = 0;

	out.write((char*) &fileheader, sizeof(fileheader));
	for (unsigned i = 0; i < 3; i++)
	{
		double pos = 0;
		for (unsigned j = 0; j < max; j++)
		{
			out.write((char*) &pos, sizeof(double));
			pos += res;
		}
	}

	//now coordinates
	//definitely not the most efficient way
	for (unsigned i = 0; i < max; i++)
	{
		for (unsigned j = 0; j < max; j++)
		{
			for (unsigned k = 0; k < max; k++)
			{
				if (root.checkCoord(tree, i, j, k, max))
					out.write(&one, 1);
				else
					out.write(&zero, 1);
			}
		}
	}
}

void MChildNode::countLeavesAtDepths(const MOctNode* tree, unsigned depth,
		vector<unsigned>& counts) const
		{
	if (counts.size() <= depth)
		counts.resize(depth + 1, 0);

	if (isLeaf)
	{
		counts[depth]++;
	}
	else
	{
		for (unsigned i = 0; i < 8; i++)
		{
			tree[node.index].children[i].countLeavesAtDepths(tree, depth + 1,
					counts);
		}
	}
}

void MappableOctTree::countLeavesAtDepths(vector<unsigned>& counts) const
		{
	counts.clear();
	root.countLeavesAtDepths(tree, 0, counts);
}

void MappableOctTree::makeGrid(MGrid& grid, float res) const
		{
	//size
	float dim = dimension;
	unsigned max = ceil(dim / res);
	grid = MGrid(dim, res);
	float start = -dim / 2;

	//now coordinates
	//definitely not the most efficient way
	for (unsigned i = 0; i < max; i++)
	{
		for (unsigned j = 0; j < max; j++)
		{
			for (unsigned k = 0; k < max; k++)
			{
				if (root.checkCoord(tree, i, j, k, max))
				{
					float x = start + i * res;
					float y = start + j * res;
					float z = start + k * res;
					grid.setPoint(x, y, z);
				}
			}
		}
	}
}
