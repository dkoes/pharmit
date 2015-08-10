/*
 * GSSTreeCreator.cpp
 *
 *  Created on: Oct 13, 2011
 *      Author: dkoes
 */

#include "GSSTreeCreator.h"
#include "DataViewers.h"
#include "TopDownPartitioner.h"
#include <climits>

using namespace boost;

//convience function for creating an indexed path name
static string nextString(filesystem::path p, const char *base, unsigned i)
{
	stringstream str;
	str << base << i;
	return filesystem::path(p / str.str()).string();
}

//creates the shape index, give that the objs and first level of trees have been created
bool GSSTreeCreator::createIndex()
{
	Timer t;
	WorkFile nexttrees;
	string nexttreesfile = filesystem::path(dbpath / "nexttrees").string();

	//partition leaves into bottom level
	//setup level
	nodes.push_back(
			WorkFile(nextString(dbpath, "level", nodes.size()).c_str()));

	vector<file_index> nodeindices;
	nodeindices.reserve(objindices.size() / 2);

	//setup next trees
	nexttrees.set(nexttreesfile.c_str());

	//map the data
	currenttrees.switchToMap();
	//this clears the indices
	unsigned long numobjs = objindices.size();
	LeafViewer leafdata(currenttrees.map->get_address(), treeindices,
			objindices);
	leveler->createNextLevel(leafdata, nodes.back().file, nodeindices,
			nexttrees.file, treeindices);

	while (nodeindices.size() > 1)
	{
		//setup curr/next trees
		currenttrees.remove();
		swap(currenttrees, nexttrees);
		currenttrees.switchToMap();

		//this stores and resets the indices
		NodeViewer nodedata(currenttrees.map->get_address(), treeindices,
				nodeindices);

		//setup next level
		nodes.push_back(
				WorkFile(nextString(dbpath, "level", nodes.size()).c_str()));

		nexttrees.set(nextString(dbpath, "trees", nodes.size()).c_str());

		leveler->createNextLevel(nodedata, nodes.back().file, nodeindices,
				nexttrees.file, treeindices);
	}
	currenttrees.remove();
	nexttrees.remove();

	cout << "Create/write index\t" << t.elapsed() << "\n";
	t.restart();

	optimizeLevels();

	cout << "Optimized levels\t" << t.elapsed() << "\n";
	//output general info
	filesystem::path infoname = dbpath / "info";
	ofstream info(infoname.string().c_str());
	info << dimension << " " << resolution << " " << nodes.size() << " "
			<< numobjs << "\n";

	//clear workfile memory
	for (unsigned i = 0, n = nodes.size(); i < n; i++)
	{
		nodes[i].clear();
	}
	nexttrees.clear();
	currenttrees.clear();
	return true;
}

//create index from existing trees, which are copied over
bool GSSTreeCreator::create(filesystem::path dir, filesystem::path treedir,
		float dim,
		float res)
{
	initialize(dir, dim, res);

	Timer t;
	const int MAXBUF = 8192;
	char buffer[MAXBUF];
	//read in trees
	filesystem::path srctreesname = treedir / "trees";
	ifstream srctrees(srctreesname.string().c_str());
	if (!srctrees)
		return false;
	streamsize sz = 0;
	while ((sz = srctrees.readsome(buffer, MAXBUF)) > 0)
	{
		currenttrees.file->write(buffer, sz);
	}

	//read in objs
	filesystem::path srcobjsname = treedir / "objs";
	ifstream srcobjs(srcobjsname.string().c_str());
	if (!srcobjs)
		return false;
	while ((sz = srcobjs.readsome(buffer, MAXBUF)) > 0)
	{
		objects.file->write(buffer, sz);
	}

	//read in treeindices
	treeindices.clear();
	filesystem::path srctiname = treedir / "treeindices";
	ifstream srcti(srctiname.string().c_str());
	if (!srcti)
		return false;
	file_index ind;
	while (srcti.read((char*) &ind, sizeof(file_index)))
	{
		treeindices.push_back(ind);
	}

	//read in objindices
	objindices.clear();
	filesystem::path srcoiname = treedir / "objindices";
	ifstream srcoi(srcoiname.string().c_str());
	if (!srcoi)
		return false;
	while (srcoi.read((char*) &ind, sizeof(file_index)))
	{
		objindices.push_back(ind);
	}
	cout << "Read/write trees\t" << t.elapsed() << "\n";
	t.restart();

	return createIndex();

}

//recursive helper for optimizing level output
//return new position
file_index GSSTreeCreator::optimizeLevelsR(ostream& outnodes,
		ostream& outleaves, const GSSNodeCommon *n, unsigned level,
		file_index& lstart, file_index& lend)
{
	if (n->isLeaf)
	{
		const GSSLeaf* leaf = (const GSSLeaf*) n;
		file_index ret = outleaves.tellp();
		outleaves.write((const char*) leaf, leaf->bytes());
		lstart = ret;
		lend = outleaves.tellp();

		numLeaves++;
		if (leaf->size() >= leafContentDistribution.size())
			leafContentDistribution.resize(leaf->size() + 1);
		leafContentDistribution[leaf->size()]++;

		return ret; //leaves are neg
	}
	else
	{
		numNodes++;
		const GSSInternalNode* node = (const GSSInternalNode*) n;
		file_index ret = outnodes.tellp();

		if (node->size() >= nodeContentDistribution.size())
			nodeContentDistribution.resize(node->size() + 1);
		nodeContentDistribution[node->size()]++;

		GSSInternalNode* newnode = node->createTruncated(dimension, resolution);

		outnodes.write((const char*) newnode, newnode->bytes());
		unsigned nextlevel = level - 1;

		lstart = ULONG_MAX;
		lend = 0;
		for (unsigned i = 0, n = newnode->size(); i < n; i++)
		{
			file_index ls = ULONG_MAX, le = 0;
			const GSSInternalNode::Child *child = newnode->getChild(i);
			const GSSNodeCommon* next =
					(const GSSNodeCommon*) ((const char*) nodes[nextlevel].map->get_address()
							+ child->position());
			file_index newpos = optimizeLevelsR(outnodes, outleaves, next,
					nextlevel, ls, le);
			lstart = min(lstart, ls);
			lend = max(lend, le);
			newnode->setChildPos(i, newpos, next->isLeaf, ls, le);
		}

		outnodes.seekp(ret, ios_base::beg);
		outnodes.write((const char*) newnode, newnode->bytes());
		outnodes.seekp(0, ios_base::end);
		free(newnode);

		return ret;
	}
}

//this combines all the internal nodes at maxlevel into newroots (which contains malloced, truncated
//copies of the nodes).
//levels count down
void GSSTreeCreator::getNodesForSuperNode(const GSSInternalNode* node,
		vector<GSSInternalNode*>& newroots, unsigned curlevel,
		unsigned stoplevel)
{
	if (curlevel == stoplevel)
	{
		newroots.push_back(node->createTruncated(dimension, resolution));
	}
	else //explore children
	{
		for (unsigned i = 0, n = node->size(); i < n; i++)
		{
			const GSSInternalNode::Child *child = node->getChild(i);
			const GSSNodeCommon* next =
					(const GSSNodeCommon*) ((const char*) nodes[curlevel - 1].map->get_address()
							+ child->position());
			assert(!next->isLeaf); // for simplicity always have root be internal node
			getNodesForSuperNode((const GSSInternalNode*) next, newroots,
					curlevel - 1, stoplevel);
		}
	}
}

//re-write levels in depth first order; leaves get their own file
void GSSTreeCreator::optimizeLevels()
{
	for (unsigned i = 0, n = nodes.size(); i < n; i++)
	{
		nodes[i].switchToMap();
	}
	numNodes = 0;
	numLeaves = 0;
	nodeContentDistribution.clear();
	nodeContentDistribution.resize(leveler->getPack() + 1, 0);
	leafContentDistribution.clear();
	leafContentDistribution.resize(leveler->getPack() + 1, 0);

	filesystem::path npath = dbpath / "nodes";
	ofstream outnodes(npath.string().c_str());

	filesystem::path lpath = dbpath / "leaves";
	ofstream outleaves(lpath.string().c_str());

	const GSSNodeCommon* root =
			(GSSNodeCommon*) nodes.back().map->get_address();

	if (root->isLeaf)
	{
		file_index ls, le;
		optimizeLevelsR(outnodes, outleaves, root, nodes.size() - 1, ls, le);

		for (unsigned i = 0, n = nodes.size(); i < n; i++)
		{
			nodes[i].remove();
		}
	}
	else //construct super node
	{
		unsigned stopLevel = 1;
		if (nodes.size() > superNodeDepth + 1)
			stopLevel = nodes.size() - superNodeDepth - 1;

		vector<GSSInternalNode*> superroots;
		getNodesForSuperNode((const GSSInternalNode*) root, superroots,
				nodes.size() - 1, stopLevel);
		GSSInternalNode* newnode = GSSInternalNode::createMergedNode(
				superroots);
		//cout << "Supernode children " << newnode->size() << "\n";
		for (unsigned i = 0, n = superroots.size(); i < n; i++)
		{
			free(superroots[i]);
		}
		superroots.clear();

		file_index ret = outnodes.tellp();
		outnodes.write((const char*) newnode, newnode->bytes());
		unsigned nextlevel = stopLevel - 1;

		file_index lstart = ULONG_MAX;
		file_index lend = 0;
		for (unsigned i = 0, n = newnode->size(); i < n; i++)
		{
			file_index ls = ULONG_MAX, le = 0;
			const GSSInternalNode::Child *child = newnode->getChild(i);
			const GSSNodeCommon* next =
					(const GSSNodeCommon*) ((const char*) nodes[nextlevel].map->get_address()
							+ child->position());
			file_index newpos = optimizeLevelsR(outnodes, outleaves, next,
					nextlevel, ls, le);
			lstart = min(lstart, ls);
			lend = max(lend, le);
			newnode->setChildPos(i, newpos, next->isLeaf, ls, le);
		}

		outnodes.seekp(ret, ios_base::beg);
		outnodes.write((const char*) newnode, newnode->bytes());
		outnodes.seekp(0, ios_base::end);
		free(newnode);

		for (unsigned i = 0, n = nodes.size(); i < n; i++)
		{
			nodes[i].remove();
		}
	}

}

//print out some distributions
void GSSTreeCreator::printStats(ostream& out) const
		{
	out << "Nodes " << numNodes << "  Leaves " << numLeaves << "\n";

	unsigned mind = min(nodeContentDistribution.size(),
			leafContentDistribution.size());
	unsigned nsum = 0, lsum = 0;
	for (unsigned i = 0; i < mind; i++)
	{
		out << i << ": " << nodeContentDistribution[i] << "\t"
				<< leafContentDistribution[i] << "\n";
		nsum += nodeContentDistribution[i];
		lsum += leafContentDistribution[i];
	}
	cout << nodeContentDistribution.size() << "\t"
			<< leafContentDistribution.size() << "\n";
}

//top down partition
void GSSLevelCreator::createNextLevel(DataViewer& data, ostream* nodefile,
		vector<file_index>& nodeindices, ostream* treefile,
		vector<file_index>& treeindices)
{
	if (data.size() == 0)
		return;

	TopDownPartitioner *thispart = partitioner->create(&data);

	packingSize = nodePack;
	if (data.isLeaf()) //making leaf nodes
		packingSize = leafPack;

	outNodes = nodefile;
	outTrees = treefile;
	nodeIndices = &nodeindices;
	treeIndices = &treeindices;
	//recursively partition
	createNextLevelR(thispart);
	delete thispart;
}

void GSSLevelCreator::createNextLevelR(TopDownPartitioner *P)
{
	if (P->size() == 0)
		return;

	if (P->size() <= packingSize)
	{
		//bottom up pack
		const DataViewer* dv = P->getData();
		vector<Cluster> clusters;
		vector<unsigned> dvindices;
		P->extractIndicies(dvindices);

		const DataViewer* data = dv->createSlice(dvindices); //reindex for better caching

		packer->pack(data, clusters);

		//each cluster creates a new node
		for (unsigned c = 0, nc = clusters.size(); c < nc; c++)
		{
			//write out node
			nodeIndices->push_back((file_index) outNodes->tellp());
			treeIndices->push_back((file_index) outTrees->tellp());
			if (data->isLeaf())
			{
				//the children are all single trees
				GSSLeaf::writeLeaf(data, clusters[c], *outNodes, *outTrees);
			}
			else
			{
				GSSInternalNode::writeNode(data, clusters[c], *outNodes,
						*outTrees);
			}
		}

		delete data;
	}
	else
	{
		vector<TopDownPartitioner*> parts;
		P->partition(parts);

		for (unsigned i = 0, n = parts.size(); i < n; i++)
		{
			createNextLevelR(parts[i]);
			delete parts[i];
		}
	}
}

