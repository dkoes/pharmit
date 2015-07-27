/*
 * GreedyPacker.cpp
 *
 *  Created on: Nov 2, 2011
 *      Author: dkoes
 */

#include "GreedyPacker.h"
#include <boost/multi_array.hpp>

using namespace lemon;


typedef bool (*ICDComp)(const IntraClusterDist& lhs,
		const IntraClusterDist& rhs);
static bool reverseDist(const IntraClusterDist& lhs,
		const IntraClusterDist& rhs)
{
	return lhs.dist > rhs.dist;
}

//a priority q for distances, roll our own to support more efficient
//removal of clusters and re-use of distances array
class DistancePQ
{
	vector<IntraClusterDist>& Q;

public:
	DistancePQ(vector<IntraClusterDist>& q) :
			Q(q)
	{
		make_heap(Q.begin(), Q.end(), reverseDist);
	}

	bool empty() const
	{
		return Q.size() == 0;
	}

	const IntraClusterDist& top() const
	{
		return Q[0];
	}

	void push(const IntraClusterDist& c)
	{
		Q.push_back(c);
		push_heap(Q.begin(), Q.end(), reverseDist);
	}

	void pop()
	{
		pop_heap(Q.begin(), Q.end(), reverseDist);
		Q.pop_back();
	}

	//remove all references to clusters i and j; this is linear
	void removeIJ(unsigned i, unsigned j)
	{
		unsigned c = 0;
		unsigned end = Q.size();
		while (c < end)
		{
			if (Q[c].i != i && Q[c].i != j && Q[c].j != i && Q[c].j != j)
			{
				c++;
			}
			else
			{
				end--;
				swap(Q[c], Q[end]);
			}
		}

		Q.erase(Q.begin() + end, Q.end());
		make_heap(Q.begin(), Q.end(), reverseDist);
	}

};

//keep track of what clusters are connected to what
class ClusterConnections
{
	vector< vector<unsigned> > table;
public:

	//assumes distances are distinct edges
	ClusterConnections(const vector<IntraClusterDist>& distances)
	{
		for(unsigned i = 0, n = distances.size(); i < n; i++)
		{
			unsigned I = distances[i].i;
			unsigned J = distances[i].j;

			unsigned m = max(I,J);
			if(table.size() <= m)
				table.resize(m+1);

			table[I].push_back(J);
			table[J].push_back(I);
		}

		for(unsigned i = 0, n = table.size(); i < n; i++)
		{
			sort(table[i].begin(), table[i].end());
		}
	}

	//return all clusters connected to c
	const vector<unsigned>& neighbors(unsigned c) const
	{
		assert(c < table.size());
		return table[c];
	}

	//merge i and j into a new cluster c
	//i and j go away, c is connected to everything they are connected to
	void merge(unsigned i, unsigned j, unsigned c)
	{
		if(table.size() <= c)
			table.resize(c+1);

		table[c].insert(table[c].end(), table[i].begin(), table[i].end());
		table[c].insert(table[c].end(), table[j].begin(), table[j].end());

		table[i].clear();
		table[j].clear();

		//uniquify
		sort(table[c].begin(), table[c].end());
		vector<unsigned>::iterator newend = unique(table[c].begin(), table[c].end());
		table[c].erase(newend, table[c].end());
	}
};

void GreedyPacker::pack(const DataViewer* dv,
		vector<Cluster>& finalclusters) const
{
	if (dv->size() == 0)
		return;

	vector<unsigned> indices(dv->size());
	for (unsigned i = 0, n = indices.size(); i < n; i++)
		indices[i] = i;

	//start with each object in its own cluster
	unsigned N = indices.size();
	vector<Cluster> clusters;
	clusters.reserve(N * 2);
	clusters.resize(N);
	for (unsigned i = 0, n = indices.size(); i < n; i++)
	{
		unsigned index = indices[i];
		clusters[i].setToSingleton(index, dv->getMIV(index), dv->getMSV(index));
	}

	vector<IntraClusterDist> distances;
	DCache *dcacheptr = NULL;

	computeDistanceVector(dv, clusters, distances, dcacheptr);

	DistancePQ pQ(distances); //manipulates a reference of distances

	//keep track of what clusters are connected for the knn case
	//this isn't really necessary for non-KNN, but that's a curiousity anyway
	ClusterConnections connections(distances);

	unsigned max = packSize;
	while (!pQ.empty())
	{
		//get the smallest distance
		unsigned i = pQ.top().i;
		unsigned j = pQ.top().j;
		pQ.pop();

		if (clusters[i].isValid() && clusters[j].isValid())
		{
			if (clusters[i].size() + clusters[j].size() <= max)
			{
				unsigned c = clusters.size();
				//merge, create a new cluster
				clusters.push_back(Cluster());
				clusters.back().mergeInto(clusters[i], clusters[j]);

				connections.merge(i, j, c);
				pQ.removeIJ(i, j);
				//now compute the distance between this cluster and all remaining clusters
				//and add to priority Q
				const vector<unsigned>& neighs = connections.neighbors(c);
				for (unsigned nbr = 0, n = neighs.size(); nbr < n; nbr++)
				{
					unsigned d = neighs[nbr];
					if (clusters[d].isValid())
					{
						float dist = clusterDistance(dv, clusters.back(),
								clusters[d], *dcacheptr);
						pQ.push(IntraClusterDist(c, d, dist));
					}
				}
			}
		}
	}

	delete dcacheptr;
	dcacheptr = NULL;

	//remove empty clusters
	finalclusters.clear();
	finalclusters.reserve(ceil(N / (double) packSize));
	for (unsigned i = 0, n = clusters.size(); i < n; i++)
	{
		if (clusters[i].isValid())
		{
			finalclusters.push_back(Cluster());
			finalclusters.back().moveInto(clusters[i]);
		}
	}
}
