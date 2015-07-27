/*
 * FullMergePacker.cpp
 *
 *  Created on: Oct 18, 2011
 *      Author: dkoes
 */

#include "FullMergePacker.h"
#include "ShapeDistance.h"

//bottom up even clustering (merge all pairs, repeat)
void FullMergePacker::pack(const DataViewer* dv, vector<Cluster>& clusters) const
{
	if (dv->size() == 0)
		return;

	vector<unsigned> indices(dv->size());
	for(unsigned i = 0, n = indices.size(); i < n; i++)
		indices[i] = i;

	//start with each object in its own cluster
	clusters.clear();
	clusters.resize(indices.size());
	for (unsigned i = 0, n = indices.size(); i < n; i++)
	{
		unsigned index = indices[i];
		clusters[i].setToSingleton(index, dv->getMIV(index), dv->getMSV(index));
	}

	DCache *dcache = NULL;
	//combine everything as much as possible
	while (fullMergeClusters(dv, clusters, dcache))
		;

	if(dcache) delete dcache;

}


bool FullMergePacker::fullMergeClusters(const DataViewer *D,
		vector<Cluster>& clusters, DCache *& dcache) const
{
	unsigned N = clusters.size();
	if (N == 1)
		return false;

	vector<IntraClusterDist> distances;
	computeDistanceVector(D, clusters, distances, dcache);

	unsigned maxclust = 0;

	sort(distances.begin(), distances.end());

	vector<Cluster> newclusters;
	newclusters.reserve(N / 2 + 1);
	unsigned merged = 0;
	bool didmerge = false;
	for (unsigned d = 0, maxdists = distances.size(); d < maxdists; d++)
	{
		//merge if not already merged
		unsigned i = distances[d].i;
		unsigned j = distances[d].j;
		if (clusters[i].isValid() && clusters[j].isValid()
				&& clusters[i].size() + clusters[j].size() <= packSize)
		{
			newclusters.push_back(Cluster());
			newclusters.back().mergeInto(clusters[i], clusters[j]);
			merged += 2;
			didmerge = true;
			if (newclusters.back().size() > maxclust)
				maxclust = newclusters.back().size();
		}

		if (merged >= N - 1)
			break;
	}

	if (merged != N) //find the loner(s), keep it as a singleton (may also be too big)
	{
		for (unsigned i = 0; i < N; i++)
		{
			if (clusters[i].isValid())
			{
				newclusters.push_back(Cluster());
				newclusters.back().moveInto(clusters[i]);
			}
		}
	}

	//this method may leave a singleton at then end, merge it
	if (!didmerge && newclusters.size() > 1 && newclusters.back().size() == 1)
	{
		//push the lone singleton into some other cluster
		double minval = HUGE_VAL;
		unsigned best = 0;
		for (unsigned i = 0, n = newclusters.size() - 1; i < n; i++)
		{
			float dist = clusterDistance(D, newclusters.back(), newclusters[i], *dcache);
			if (dist < minval)
			{
				minval = dist;
				best = i;
			}
		}

		newclusters[best].addInto(newclusters.back());
		newclusters.pop_back();
	}

	swap(newclusters, clusters);

	return didmerge;
}

