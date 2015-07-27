/*
 * MatcherPacker.cpp
 *
 *  Created on: Nov 2, 2011
 *      Author: dkoes
 */

#include "MatcherPacker.h"
#include <boost/unordered_set.hpp>
#include <lemon/full_graph.h>
#include <lemon/matching.h>
#include <lemon/list_graph.h>
#include <lemon/connectivity.h>
#include <lemon/graph_to_eps.h>
#include <lemon/lp.h>
#include <ANN/ANN.h>
#include "Timer.h"
using namespace lemon;

void MatcherPacker::pack(const DataViewer* dv, vector<Cluster>& clusters) const
{
	if (dv->size() == 0)
		return;

	vector<unsigned> indices(dv->size());
	for (unsigned i = 0, n = indices.size(); i < n; i++)
		indices[i] = i;

	//start with each object in its own cluster
	clusters.clear();
	clusters.resize(indices.size());
	for (unsigned i = 0, n = indices.size(); i < n; i++)
	{
		unsigned index = indices[i];
		clusters[i].setToSingleton(index, dv->getMIV(index), dv->getMSV(index));
	}

	//combine everything as much as possible
	unsigned curSz = 1;

	if (K == 0)
	{
		FullCache dcache(dv);

		while (curSz < packSize
				&& fullMergeClusters(dv, clusters, curSz, dcache))
			curSz *= 2;
	}
	else //use knn
	{
		OnDemandCache dcache(dv);

		while (curSz < packSize && knnMergeClusters(dv, clusters, curSz, dcache))
			curSz *= 2;

		//cout << "dcache " << dcache.size() << "\n";
	}

}


/* find the pseudo-optimal matching of every pair of clusters
 * use the knn graph instead of the full graph
 * maxSz is what the largest current cluster should be
 */
bool MatcherPacker::knnMergeClusters(const DataViewer *D,
		vector<Cluster>& clusters, unsigned maxSz, DCache& dcache) const
{
	unsigned N = clusters.size();
	SmartGraph SG;
	SmartGraph::EdgeMap<double> Sweights(SG);
	SmartGraph::NodeMap<unsigned> Sindex(SG);
	vector<SmartGraph::Node> Snodes;

	ClusterCache ccache(clusters.size());
	makeKNNGraph(D, clusters, maxSz, dcache, ccache, SG, Sweights, Sindex,
			Snodes);

	Timer t;
	MaxWeightedMatching<SmartGraph, SmartGraph::EdgeMap<double> > matcher(SG,
			Sweights);
	matcher.run();
	//merge clusters that were matched
	vector<Cluster> newclusters;
	newclusters.reserve(N / 2 + 1);
	vector<bool> packed(N, false);
	bool didmerge = false;
	for (unsigned i = 0; i < N; i++)
	{
		if (!packed[i])
		{
			SmartGraph::Node m = matcher.mate(Snodes[i]);
			if (!SG.valid(m))
			{
				newclusters.push_back(Cluster());
				newclusters.back().moveInto(clusters[i]);
			}
			else
			{
				unsigned j = Sindex[m];

				if (clusters[i].size() + clusters[j].size() <= packSize) //merge
				{
					newclusters.push_back(Cluster());
					newclusters.back().mergeInto(clusters[i], clusters[j]);
					didmerge = true;
					packed[i] = true;
					packed[j] = true;
				}
				else
				{
					//didn't merge, keep clusters as is
					newclusters.push_back(Cluster());
					newclusters.back().moveInto(clusters[i]);
					newclusters.push_back(Cluster());
					newclusters.back().moveInto(clusters[j]);
					packed[i] = true;
					packed[j] = true;
				}
			}
		}
	}

	swap(clusters, newclusters);
	return didmerge;
}

//write dot file from graph
void writeDot(SmartGraph& g, SmartGraph::EdgeMap<double>& length,
		SmartGraph::NodeMap<unsigned>& index,
		ostream& out)
{
	out << "graph name {" << endl;
	out << "  node [ shape=ellipse, fontname=Helvetica, fontsize=10 ];" << endl;
	for (SmartGraph::NodeIt n(g); n != INVALID; ++n)
	{
		out << "  n" << g.id(n) << " [ label=\"" << index[n] << "\" ]; " << endl;
	}
	out << "  edge [ shape=ellipse, fontname=Helvetica, fontsize=10 ];" << endl;
	for (SmartGraph::EdgeIt e(g); e != INVALID; ++e)
	{
		out << "  n" << g.id(g.u(e)) << " -- " << " n" << g.id(g.v(e))
				<< " [ label=\"" << length[e] << "\" ]; " << endl;
	}
	out << "}" << endl;

}

typedef pair<unsigned, unsigned> IndexPair;

void writeEGDot(SmartGraph& g, SmartGraph::EdgeMap<double>& length,
		ostream& out, SmartGraph::NodeMap<IndexPair>& edgePairs,
		boost::unordered_set<unsigned>& matched)
{
	out << "graph name {" << endl;
	out << "  node [ shape=ellipse, fontname=Helvetica, fontsize=10 ];" << endl;
	for (SmartGraph::NodeIt n(g); n != INVALID; ++n)
	{
		IndexPair p = edgePairs[n];
		if(p.first != p.second)
			out << "  n" << g.id(n) << " [ label=\"" << p.first << ":" << p.second << "\" ]; " << endl;
		else
			out << "  n" << g.id(n) << " [ shape=triangle, label=\"" << p.first << "\" ]; \n";
	}
	out << "  edge [ shape=ellipse, fontname=Helvetica, fontsize=10 ];" << endl;
	for (SmartGraph::EdgeIt e(g); e != INVALID; ++e)
	{
		string extra = "";
		if(matched.count(g.id(e)))
			extra = "style = bold, ";
		out << "  n" << g.id(g.u(e)) << " -- " << " n" << g.id(g.v(e))
				<< " [ " << extra << "label=\"" << length[e] << "\" ]; " << endl;
	}
	out << "}" << endl;
}

/* find the optimal matching of every pair of clusters
 * maxSz is what the largest current cluster should be
 */
bool MatcherPacker::fullMergeClusters(const DataViewer *D,
		vector<Cluster>& clusters, unsigned maxSz, DCache& dcache) const
{
	unsigned C = clusters.size();
	if (C == 1)
		return false;

	//create full graph with distances; enforce an even number of nodes if
	//necessary by adding a dummy node
	unsigned N = C;
	if (C % 2)
		N++;

	FullGraph G(N);
	FullGraph::EdgeMap<double> weights(G);
	for (unsigned i = 0; i < C; i++)
	{
		for (unsigned j = 0; j < i; j++)
		{
			float dist = clusterDistance(D, clusters[i], clusters[j], dcache);
			weights[G.arc(G(i), G(j))] = -dist; //because the algo maximizes
		}
	}

	//init dummy node
	for (unsigned i = C; i < N; i++)
	{
		for (unsigned j = 0; j < C; j++)
		{
			if (clusters[j].size() < maxSz)
			{
				//force smaller clusters to merge
				weights[G.arc(G(i), G(j))] = -10000;
				weights[G.arc(G(j), G(i))] = -10000;
			}
			else
			{
				weights[G.arc(G(i), G(j))] = 0;
				weights[G.arc(G(j), G(i))] = 0;
			}
		}
	}

	MaxWeightedPerfectMatching<FullGraph, FullGraph::EdgeMap<double> > matcher(
			G, weights);
	matcher.run();

	//merge clusters that were matched
	vector<Cluster> newclusters;
	newclusters.reserve(N / 2 + 1);
	vector<bool> packed(N, false);
	bool didmerge = false;
	for (unsigned i = 0; i < C; i++)
	{
		if (!packed[i])
		{
			unsigned j = G.index(matcher.mate(G(i)));
			if (j == C)
			{
				newclusters.push_back(Cluster());
				newclusters.back().moveInto(clusters[i]);
			}
			else
			{
				if (clusters[i].size() + clusters[j].size() <= packSize) //merge
				{
					newclusters.push_back(Cluster());
					newclusters.back().mergeInto(clusters[i], clusters[j]);
					didmerge = true;
					packed[i] = true;
					packed[j] = true;
				}
				else
				{
					//didn't merge, keep clusters as is
					newclusters.push_back(Cluster());
					newclusters.back().moveInto(clusters[i]);
					newclusters.push_back(Cluster());
					newclusters.back().moveInto(clusters[j]);
					packed[i] = true;
					packed[j] = true;
				}
			}
		}
	}

	swap(clusters, newclusters);
	return didmerge;
}
