#include "GSSTreeStructures.h"
#include "DataViewers.h"

//write out a leaf node to outNodes and the representative trees to outTrees
void GSSLeaf::writeLeaf(const DataViewer *data, const Cluster& cluster, ostream& outNodes, ostream& outTrees)
{
	GSSNodeCommon info;
	info.isLeaf = true;
	info.N = cluster.size();

	unsigned positions[info.N];

	outNodes.write((char*)&info, sizeof(info));

	//compute the position offsets
	unsigned curoffset = info.N*sizeof(unsigned); //offsets are from start of data section
	for(unsigned i = 0, n = cluster.size(); i < n; i++)
	{
		unsigned ind = cluster[i];
		positions[i] = curoffset;
		const MappableOctTree *tree = data->getMIV(ind);
		curoffset += sizeof(file_index);
		curoffset += tree->bytes();
	}

	//write positions
	outNodes.write((char*)positions, sizeof(unsigned)*info.N);

	//write tree info
	for(unsigned i = 0, n = cluster.size(); i < n; i++)
	{
		unsigned ind = cluster[i];
		//for leaf, write out object file index
		file_index index = data->getIndex(ind);
		outNodes.write((char*)&index, sizeof(file_index));

		//and single tree
		const MappableOctTree *tree = data->getMIV(ind);
		tree->write(outNodes);
	}

	//write out trees of cluster to a GSSDoubleTree
	unsigned offset = cluster.MIV->bytes();

	outTrees.write((char*)&offset, sizeof(offset));
	cluster.MIV->write(outTrees);
	cluster.MSV->write(outTrees);

}

//very similar to write leaf, but write out MIV/MSV of subnodes
void GSSInternalNode::writeNode(const DataViewer *data, const Cluster& cluster, ostream& outNodes, ostream& outTrees)
{
	GSSNodeCommon info;
	info.isLeaf = false;
	info.N = cluster.size();

	unsigned positions[info.N];
	outNodes.write((char*)&info, sizeof(info));

	//compute the position offsets
	unsigned curoffset = info.N*sizeof(unsigned); //offsets are from start of data section
	for(unsigned i = 0, n = cluster.size(); i < n; i++)
	{
		unsigned ind = cluster[i];
		positions[i] = curoffset;
		const MappableOctTree *miv = data->getMIV(ind);
		const MappableOctTree *msv = data->getMSV(ind);
		curoffset += sizeof(Child);
		curoffset += miv->bytes();
		curoffset += msv->bytes();
	}

	//write positions
	outNodes.write((char*)positions, sizeof(unsigned)*info.N);

	//write tree info
	for(unsigned i = 0, n = cluster.size(); i < n; i++)
	{
		unsigned ind = cluster[i];
		const MappableOctTree *miv = data->getMIV(ind);
		const MappableOctTree *msv = data->getMSV(ind);

		Child child;
		child.MSVindex = miv->bytes();
		child.node_pos = data->getIndex(ind);
		outNodes.write((char*)&child, sizeof(Child));

		miv->write(outNodes);
		msv->write(outNodes);
	}

	//write out trees of cluster to a GSSDoubleTree
	unsigned offset = cluster.MIV->bytes();

	outTrees.write((char*)&offset, sizeof(offset));
	cluster.MIV->write(outTrees);
	cluster.MSV->write(outTrees);

}

//create a single super node from nodes as a malloced node
GSSInternalNode* GSSInternalNode::createMergedNode(const vector<GSSInternalNode*>& nodes)
{
	unsigned numChildren = 0;
	vector<unsigned> positions;
	unsigned curoffset = 0;
	for(unsigned i = 0, n = nodes.size(); i < n; i++)
	{
		unsigned nc = nodes[i]->size();
		numChildren += nc;

		for(unsigned c = 0; c < nc; c++)
		{
			positions.push_back(curoffset);
			curoffset += nodes[i]->getChild(c)->bytes();
		}
	}

	//add in position offset
	unsigned posoff = positions.size()*sizeof(unsigned);
	for(unsigned i = 0, n = positions.size(); i < n; i++)
	{
		positions[i] += posoff;
	}
	GSSInternalNode *ret = (GSSInternalNode*)malloc(sizeof(GSSInternalNode)+curoffset+posoff);

	ret->info.isLeaf = false;
	ret->info.N = numChildren;
	unsigned offset = 0;
	memcpy(ret->data, &positions[0], posoff);
	offset += posoff;

	//copy in children
	for(unsigned i = 0, n = nodes.size(); i < n; i++)
	{
		unsigned nc = nodes[i]->size();
		for(unsigned c = 0; c < nc; c++)
		{
			const Child* ch = nodes[i]->getChild(c);
			unsigned nb = ch->bytes();
			memcpy(ret->data+offset, ch, nb);
			offset += nb;
		}
	}
	assert(offset == curoffset+posoff);

	return ret;
}


unsigned GSSLeaf::bytes() const
{
	unsigned ret = sizeof(GSSLeaf);
	if(info.N == 0) return ret;
	ret += child_positions[info.N-1]; //offset to last child
	const Child* child = (const Child*)&data[child_positions[info.N-1]];
	//add size of last child
	ret += sizeof(file_index);
	ret += child->tree.bytes();
	return ret;
}

unsigned GSSInternalNode::Child::bytes() const
{
	unsigned ret = sizeof(Child);
	ret += MSVindex;
	const MappableOctTree* msv = getMSV();
	ret += msv->bytes();

	return ret;
}

unsigned GSSInternalNode::bytes() const
{
	unsigned ret = sizeof(GSSInternalNode);
	if(info.N == 0) return ret;
	ret += child_positions[info.N-1]; //offset to last child
	const Child* child = (const Child*)&data[child_positions[info.N-1]];
	//add size of last child
	ret += child->bytes();
	return ret;
}

void GSSInternalNode::setChildPos(unsigned i, file_index newpos, bool isLeaf, file_index lstart, file_index lend)
{
	Child *child =  (Child*)&data[child_positions[i]];
	child->node_pos = newpos;
	child->isLeaf = isLeaf;
	child->leaves_start = lstart;
	child->leaves_end = lend;
}

//malloc a new node that has reduced sized trees
GSSInternalNode* GSSInternalNode::createTruncated(float dimension, float res) const
{
	unsigned nc = size();
	MappableOctTree *MIVs[nc];
	MappableOctTree *MSVs[nc];
	const MappableOctTree *oldMIVs[nc];
	const MappableOctTree *oldMSVs[nc];
	unsigned positions[nc];

	//reduce resolution to lowest discernable value
	/*
	float volcut = res*res*res/2.0;
	while(volcut < dimension*dimension*dimension)
	{
		volcut *= 2;

		for (unsigned c = 0; c < nc; c++)
		{
			const Child *child = getChild(c);
			MIVs[c] = child->getMIV()->createRounded(dimension, volcut,
					false);
			MSVs[c] = child->getMSV()->createRounded(dimension, volcut,
					true);
		}

		bool hasExactMatch = false;
		for(unsigned i = 0; i < nc && !hasExactMatch; i++)
		{
			for(unsigned j = 0; j < i && !hasExactMatch; j++)
			{
				if(MIVs[i]->equals(MIVs[j]))
					hasExactMatch = true;
				else if(MSVs[i]->equals(MSVs[j]))
					hasExactMatch = true;
			}
		}

		for (unsigned c = 0; c < nc; c++)
		{
			free(MIVs[c]);
			free(MSVs[c]);
		}

		if(hasExactMatch)
		{
			volcut /= 2;
			break;
		}
	}
*/
	unsigned curoffset = nc*sizeof(unsigned); //skip of positions

	for(unsigned c = 0; c < nc; c++)
	{
		const Child *child = getChild(c);
		oldMIVs[c] = child->getMIV();
		oldMSVs[c] = child->getMSV();
	}

	bool changed = false;

	changed |= MappableOctTree::createRoundedSet(nc, oldMIVs, false, MIVs);
	changed |= MappableOctTree::createRoundedSet(nc, oldMSVs, true, MSVs);

	while(changed)
	{
		changed = false;
		memcpy(oldMIVs,MIVs, sizeof(MIVs));
		memcpy(oldMSVs,MSVs, sizeof(MSVs));

		changed |= MappableOctTree::createRoundedSet(nc, oldMIVs, false, MIVs);
		changed |= MappableOctTree::createRoundedSet(nc, oldMSVs, true, MSVs);

		for(unsigned c = 0; c < nc; c++)
		{
			free((void*)oldMIVs[c]);
			free((void*)oldMSVs[c]);
		}
	}

	for(unsigned c = 0; c < nc; c++)
	{
		positions[c] = curoffset;

		curoffset += sizeof(Child);
		curoffset += MIVs[c]->bytes();
		curoffset += MSVs[c]->bytes();
	}

	unsigned char* buffer = (unsigned char*)malloc(curoffset+sizeof(GSSInternalNode));

	memcpy(buffer, this, sizeof(GSSInternalNode));
	unsigned offset = sizeof(GSSInternalNode);
	memcpy(buffer+offset, &positions, sizeof(unsigned)*nc);
	offset += sizeof(unsigned)*nc;

	for(unsigned c = 0; c < nc; c++)
	{
		Child child = *getChild(c);
		child.MSVindex = MIVs[c]->bytes();

		memcpy(buffer+offset, &child, sizeof(Child));
		offset += sizeof(Child);
		memcpy(buffer+offset, MIVs[c], MIVs[c]->bytes());
		offset += MIVs[c]->bytes();
		memcpy(buffer+offset, MSVs[c], MSVs[c]->bytes());
		offset += MSVs[c]->bytes();

		free(MIVs[c]);
		free(MSVs[c]);
	}

	assert(offset == curoffset+sizeof(GSSInternalNode));
	return (GSSInternalNode*)buffer;
}

