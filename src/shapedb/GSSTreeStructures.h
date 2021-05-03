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
 * GSSTreeStructures.h
 *
 *  Created on: Oct 14, 2011
 *      Author: dkoes
 *
 *  This contains data types/structures that are shared between the creator and the searcher
 *
 *  These can be mapped into a file.
 */

#ifndef GSSTREESTRUCTURES_H_
#define GSSTREESTRUCTURES_H_

#include "MappableOctTree.h"
#include "GSSTypes.h"
#include <cassert>


struct GSSDoubleTree
{
	unsigned msvOffset;
	unsigned char data[]; //miv than msv

	//write out tree data in a manner consistent with structure layout
	static void writeTreeS(ostream& out, MappableOctTree *miv, MappableOctTree *msv)
	{
		unsigned p = miv->bytes();
		out.write((char*)&p, sizeof(unsigned));
		miv->write(out);
		msv->write(out);
	}

	const MappableOctTree* getMIV() const
	{
		return (const MappableOctTree*)data;
	}

	const MappableOctTree* getMSV() const
	{
		return (const MappableOctTree*)&data[msvOffset];
	}
};

//header of leaf and internal nodes
struct GSSNodeCommon
{
	bool isLeaf: 1;
	unsigned N: 31;
};

//a GSSLeaf only needs to store a single tree for each object, and the positions
//are within the object file, not the nodes file
class GSSLeaf
{
public:
	struct Child
	{
		file_index object_pos;
		MappableOctTree tree;
	} __attribute__((__packed__));
private:
	GSSNodeCommon info;
	unsigned child_positions[]; //variable length, offset of trees within data section
	
public:
	static void writeLeaf(const DataViewer *data, const Cluster& cluster, ostream& outNodes, ostream& outTrees);
	const Child* getChild(unsigned i) const {
#pragma GCC diagnostic ignored "-Warray-bounds"
		unsigned char *d = data();
		return (Child*)&d[child_positions[i]];
	}
	unsigned char *data() const { return (unsigned char*)&child_positions; }
	unsigned size() const { return info.N; }
	unsigned bytes() const;
} __attribute__((__packed__));

//a GSSInternalNode stores both the MIV and MSV for each subnode, and points to their locations
class GSSInternalNode
{
public:
	struct Child
	{
		friend class GSSInternalNode;
	private:
		bool isLeaf: 1;
		file_index node_pos: 63;
		file_index leaves_start;
		file_index leaves_end;
		unsigned MSVindex; //where in data MSV starts (MIV at 0)
		unsigned char data[];
	public:
		Child(): isLeaf(false), node_pos(0), leaves_start(0), leaves_end(0), MSVindex(0) {}

		const MappableOctTree* getMIV() const { return (const MappableOctTree*)data; }
		const MappableOctTree* getMSV() const { return (MappableOctTree*)&data[MSVindex];}

		file_index position() const { return node_pos; }
		bool isLeafPosition() const { return isLeaf; }
		unsigned bytes() const;
	} __attribute__((__packed__));
private:
	GSSNodeCommon info;
	unsigned child_positions[]; //offset of children within data section

public:
	static void writeNode(const DataViewer *data, const Cluster& cluster, ostream& outNodes, ostream& outTrees);
	static GSSInternalNode* createMergedNode(const vector<GSSInternalNode*>& nodes);
	const Child* getChild(unsigned i) const { 
		unsigned char *d = data();
		return (Child*)&d[child_positions[i]]; 
	}
	unsigned size() const { return info.N; }
	unsigned bytes() const;
	unsigned char *data() const { return (unsigned char*)&child_positions; }

	//truncates MIV/MSV of children, returns result in malloced memory
	GSSInternalNode* createTruncated(float dimension, float resolution) const;

	void setChildPos(unsigned i, file_index newpos, bool isLeaf, file_index lstart, file_index lend);
} __attribute__((__packed__));






#endif /* GSSTREESTRUCTURES_H_ */
