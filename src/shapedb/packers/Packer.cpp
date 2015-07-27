/*
 * Packer.cpp
 *
 *  Created on: Nov 2, 2011
 *      Author: dkoes
 */

#include "Packer.h"
#include "ShapeDistance.h"

#include <boost/unordered_set.hpp>
#include <lemon/full_graph.h>
#include <lemon/matching.h>
#include <lemon/list_graph.h>
#include <lemon/connectivity.h>
#include <lemon/graph_to_eps.h>
#include <lemon/lp.h>
#include <ANN/ANN.h>

using namespace lemon;

//return the distance between two clusters, this may be configurable to other metrics
float Packer::clusterDistance(const DataViewer* D, const Cluster& a,
		const Cluster& b, DCache& dcache) const
{
	switch (distMetric)
	{
	case AverageLink:
		if(a.size() == 1 && b.size() == 1)
		{
			return dcache.get(a[0],b[0]);
		}
		else
		{
			return shapeDistance(a.MIV, a.MSV, b.MIV, b.MSV);
		}
	case CompleteLink:
	{
		//this is the maximum o distance between any two cluster members
		float max = 0;
		for (unsigned i = 0, ni = a.size(); i < ni; i++)
		{
			for (unsigned j = 0, nj = b.size(); j < nj; j++)
			{
				unsigned l = a[i];
				unsigned r = b[j];
				float dist = dcache.get(l,r);

				if (dist > max)
					max = dist;
			}
		}
		return max;
	}
	case TotalLink:
	{
		//this is the total sum of all the linkages, distiguished from average link
		//in that smaller clusters end up with smaller values
		float sum = 0;
		for (unsigned i = 0, ni = a.size(); i < ni; i++)
		{
			for (unsigned j = 0, nj = b.size(); j < nj; j++)
			{
				unsigned l = a[i];
				unsigned r = b[j];
				float dist = dcache.get(l,r);

				sum += dist;
			}
		}
		return sum;
	}
	case SingleLink:
	{
		//the minimum distance overall
		float min = HUGE_VAL;
		for (unsigned i = 0, ni = a.size(); i < ni; i++)
		{
			for (unsigned j = 0, nj = b.size(); j < nj; j++)
			{
				unsigned l = a[i];
				unsigned r = b[j];
				float dist = dcache.get(l,r);

				if (dist < min)
					min = dist;
			}
		}
		return min;
	}
	default:
		abort();
		break;
	}
	return 0;
}


//create reverse mapping from forward mapping
static void reverse(const vector<vector<unsigned> >& F,
		vector<vector<unsigned> >& R)
{
	R.clear();
	R.resize(F.size());
	for (unsigned v = 0, n = F.size(); v < n; v++)
	{
		for (unsigned i = 0, m = F[v].size(); i < m; i++)
		{
			unsigned neigh = F[v][i];
			R[neigh].push_back(v);
		}
	}

	for (unsigned v = 0, n = R.size(); v < n; v++)
	{
		sort(R[v].begin(), R[v].end());
	}
}

//insert in sorted order, make K items, no limit if 0
bool Packer::KNNSlice::update(const IndDist& origitem, unsigned k)
{
	IndDist item(origitem);

	for (unsigned i = 0, n = neighbors.size(); i < n; i++)
	{
		if(!isfinite(item.dist))
			abort();
		if (neighbors[i].j == item.j && neighbors[i].dist != item.dist)
		{
			//honestly, have no idea how this happens, only seems to happen with
			//very small (more than 10 sig figs) differences in distances
			//but in attempt at a slap-dash bugfix, always return false
			fprintf(stderr,"Inconsistent distances: %.12f %.12f\n",neighbors[i].dist, item.dist);
			return false;
//			abort();
		}
	}

	vector<IndDist>::iterator pos = lower_bound(neighbors.begin(),
			neighbors.end(), item);

	if (k > 0 && (pos - neighbors.begin()) >= k)
		return false;

	for (vector<IndDist>::iterator look = pos;
			look != neighbors.end() && look->dist == pos->dist; look++)
	{
		if (look->j == item.j)
			return false; //already in here
		if (look->dist == item.dist) //not actually better
			return false;
	}

	neighbors.insert(pos, item);
	if (k > 0)
		neighbors.resize(k);

	return true;
}

//processed and unprocessed; once extracted becomes processed
void Packer::KNNSlice::getOldNew(vector<unsigned>& O,
		vector<unsigned>& N)
{
	for (unsigned i = 0, n = neighbors.size(); i < n; i++)
	{
		if (neighbors[i].unprocessed)
		{
			N.push_back(neighbors[i].j);
			neighbors[i].unprocessed = false;
		}
		else
			O.push_back(neighbors[i].j);
	}
}

//create an initial knn graph quickly
void Packer::initialKNNSample(const DataViewer *D,
		vector<Cluster>& clusters, unsigned maxSz, DCache& dcache,
		ClusterCache& ccache, vector<KNNSlice>& V, float& maxdist) const
{
	if (S == 0) //random sampling
	{
		unsigned N = clusters.size();
		V.resize(N);

		vector<unsigned> indices(N);
		for (unsigned i = 0; i < N; i++)
		{
			indices[i] = i;
		}

		random_shuffle(indices.begin(), indices.end());

		unsigned pos = 0;
		for (unsigned i = 0; i < N; i++)
		{
			V[i].neighbors.reserve(K);
			for (unsigned j = 0; j < K; j++)
			{
				//get random neighbor
				unsigned neigh = indices[pos];
				pos = (pos + 1) % N;
				if (neigh == i)
				{
					neigh = indices[pos];
					pos = (pos + 1) % N;
				}

				float dist = 0;
				if (!ccache.isCached(i, neigh, dist))
				{
					dist = clusterDistance(D, clusters[i], clusters[neigh],
							dcache);
					ccache.set(i, neigh, dist);
				}

				if (dist > maxdist)
					maxdist = dist;
				V[i].neighbors.push_back(IndDist(neigh, dist));
			}

			sort(V[i].neighbors.begin(), V[i].neighbors.end());
		}
	}
	else
	{
		//pick S sentinals
		unsigned N = clusters.size();
		V.resize(N);

		ANNpointArray points = annAllocPts(N, S);

		//farthest first, initialize point data as we go
		vector<unsigned> indices;
		indices.push_back(0);
		for (unsigned s = 1; s < S; s++)
		{
			float max = 0;
			unsigned besti = 0;
			for (unsigned i = 0; i < N; i++)
			{
				float min = HUGE_VAL;
				unsigned send = indices.size() - 1; //all but last have distances already computed
				for (unsigned j = 0; j < send; j++)
				{
					float dist = points[i][j];
					if (dist < min)
						min = dist;
				}
				float dist = 0;
				if (i != indices.back())
				{
					if (!ccache.isCached(i, indices.back(), dist))
					{
						dist = clusterDistance(D, clusters[i],
								clusters[indices.back()], dcache);
						ccache.set(i, indices.back(), dist);
					}
				}
				points[i][send] = dist;
				if (dist < min)
					min = dist;

				if (min > max)
				{
					besti = i;
					max = min;
				}
			}
			indices.push_back(besti);
		}assert(indices.size() == S);
		//get distances for last chosen sentinal
		for (unsigned i = 0; i < N; i++)
		{
			float dist = 0;
			if (i != indices.back())
			{
				if (!ccache.isCached(i, indices.back(), dist))
				{
					dist = clusterDistance(D, clusters[i],
							clusters[indices.back()], dcache);
					ccache.set(i, indices.back(), dist);
				}
			}
			points[i][S - 1] = dist;
		}

		ANNkd_tree searcher(points, N, S);

		ANNdistArray dists = new ANNdist[K + 1];
		vector<int> nnIdx(K + 1, 0);
		for (unsigned i = 0; i < N; i++)
		{
			V[i].neighbors.reserve(K);
			searcher.annkSearch(points[i], K + 1, &nnIdx[0], dists);

			for (unsigned j = 0; j <= K; j++)
			{
				if (nnIdx[j] != (int) i && nnIdx[j] >= 0)
				{
					unsigned neigh = nnIdx[j];
					float dist = 0;
					if (!ccache.isCached(i, neigh, dist))
					{
						dist = clusterDistance(D, clusters[i], clusters[neigh],
								dcache);
						ccache.set(i, neigh, dist);
					}

					V[i].neighbors.push_back(IndDist(nnIdx[j], dist));

					if (dist > maxdist)
						maxdist = dist;
				}
			}
		}

		delete[] dists;
		annDeallocPts(points);
	}
}

//create a knn graph, return max value of any edge since we reverse the meanings of distances
float Packer::makeKNNGraph(const DataViewer *D, vector<Cluster>& clusters,
		unsigned maxsz, DCache& dcache, ClusterCache& ccache, SmartGraph& G,
		SmartGraph::EdgeMap<double>& E, SmartGraph::NodeMap<unsigned>& Sindex,
		vector<SmartGraph::Node>& nodes) const
{
	G.clear();
	unsigned N = clusters.size();
	float max = 0;

	if (K >= clusters.size()) //full knn
	{
		vector<vector<IndDist> > dists(N, vector<IndDist>(N));
		for (unsigned i = 0; i < N; i++)
		{
			for (unsigned j = 0; j < i; j++)
			{
				float dist = 0;

				if (!ccache.isCached(i, j, dist))
				{
					dist = clusterDistance(D, clusters[i], clusters[j], dcache);
					ccache.set(i, j, dist);
				}
				unsigned minsz = min(clusters[i].size(), clusters[j].size());
				while (minsz < maxsz)
				{
					dist /= 2;
					minsz *= 2;
				}
				dists[i][j] = IndDist(j, dist);
				dists[j][i] = IndDist(i, dist);

				if (dist > max)
					max = dist;
			}
			dists[i][i] = IndDist(i, HUGE_VAL);
		}

		G.reserveNode(N);
		G.reserveEdge(N * K);

		nodes.resize(N);
		for (unsigned i = 0; i < N; i++)
		{
			nodes[i] = G.addNode();
			Sindex[nodes[i]] = i;
		}

		for (unsigned i = 0; i < N; i++)
		{
			sort(dists[i].begin(), dists[i].end());

			for (unsigned j = 0; j < K && j < N; j++)
			{
				if (dists[i][j].j != i)
				{
					SmartGraph::Edge e = G.addEdge(nodes[i],
							nodes[dists[i][j].j]);
					E[e] = max - dists[i][j].dist; //compute max weighting
				}
			}
		}
	}
	else //approximate knn
	{
		vector<KNNSlice> B(N);

		initialKNNSample(D, clusters, maxsz, dcache, ccache, B, max);

		bool keepgoing = true;
		while (keepgoing)
		{
			keepgoing = false;
			unsigned changed = 0;
			vector<vector<unsigned> > oldf(N), newf(N);
			vector<vector<unsigned> > oldr(N), newr(N);
			for (unsigned v = 0; v < N; v++)
			{
				B[v].getOldNew(oldf[v], newf[v]);
			}

			reverse(oldf, oldr);
			reverse(newf, newr);

			for (unsigned v = 0; v < N; v++)
			{
				//merge forward and reverse neighbors
				vector<unsigned> oldn, newn;
				sort(oldf[v].begin(), oldf[v].end());
				sort(oldr[v].begin(), oldr[v].end());
				sort(newf[v].begin(), newf[v].end());
				sort(newr[v].begin(), newr[v].end());

				merge(oldf[v].begin(), oldf[v].end(), oldr[v].begin(),
						oldr[v].end(), back_inserter(oldn));
				merge(newf[v].begin(), newf[v].end(), newr[v].begin(),
						newr[v].end(), back_inserter(newn));

				//compare new to new and new to old
				for (unsigned i = 0, n = newn.size(); i < n; i++)
				{
					unsigned u1 = newn[i];
					for (unsigned j = 0; j < n; j++)
					{
						unsigned u2 = newn[j];
						if (u1 < u2)
						{
							float dist = 0;
							if (!ccache.isCached(u1, u2, dist))
							{
								dist = clusterDistance(D, clusters[u1],
										clusters[u2], dcache);
								ccache.set(u1, u2, dist);
							}

							changed += B[u1].update(IndDist(u2, dist), K);
							changed += B[u2].update(IndDist(u1, dist), K);

							if (dist > max)
								max = dist;
						}
					}

					for (unsigned j = 0, m = oldn.size(); j < m; j++)
					{
						unsigned u2 = oldn[j];
						if (u2 != u1)
						{
							float dist = 0;
							if (!ccache.isCached(u1, u2, dist))
							{
								dist = clusterDistance(D, clusters[u1],
										clusters[u2], dcache);
								ccache.set(u1, u2, dist);
							}

							if (B[u1].update(IndDist(u2, dist), K))
								changed++;
							if (B[u2].update(IndDist(u1, dist), K))
								changed++;

							if (dist > max)
								max = dist;
						}
					}
				}
			}

			if (changed > 0)
				keepgoing = true;
		}

		G.reserveNode(N);
		G.reserveEdge(N * K);

		nodes.resize(N);
		for (unsigned i = 0; i < N; i++)
		{
			nodes[i] = G.addNode();
			Sindex[nodes[i]] = i;
		}

		typedef pair<unsigned, unsigned> UPair;
		boost::unordered_set<UPair> seen;
		for (unsigned i = 0; i < N; i++)
		{
			for (unsigned j = 0, nn = B[i].neighbors.size(); j < nn; j++)
			{
				unsigned neigh = B[i].neighbors[j].j;

				if (seen.count(UPair(i, neigh)) == 0)
				{
					seen.insert(UPair(i, neigh));
					seen.insert(UPair(neigh, i));
					float dist = B[i].neighbors[j].dist;
					SmartGraph::Edge e = G.addEdge(nodes[i], nodes[neigh]);
					E[e] = max - dist; //compute max weighting
				}
			}
		}
	}
	return max;
}

//compute a vector of all the edge/distances we care about between these clusters
//will be the full n^2 list if Knn is off, otherwise constructs the Knn graph
//this abstracts away the knn graph creation for those packers who don't need to
//work on the actual graph
//initializes dcache if necessary
void Packer::computeDistanceVector(const DataViewer* dv, vector<Cluster>& clusters, vector<IntraClusterDist>& distances, DCache *& dcache) const
{
	distances.clear();
	unsigned N = clusters.size();
	if (K == 0)
	{
		//compute all intra cluster distances
		distances.reserve(N * N / 2);
		boost::multi_array<float, 2> darray(boost::extents[N][N]);
		if(dcache == NULL)
			dcache = new FullCache(dv);

		for (unsigned i = 0; i < N; i++)
		{
			darray[i][i] = 0;
			for (unsigned j = 0; j < i; j++)
			{
				float dist = clusterDistance(dv, clusters[i], clusters[j],
						*dcache);
				darray[i][j] = darray[j][i] = dist;
				distances.push_back(IntraClusterDist(i, j, dist));
			}
		}

	}
	else
	{
		if(dcache == NULL)
			dcache = new OnDemandCache(dv);

		ClusterCache ccache(clusters.size());
		SmartGraph SG;
		SmartGraph::EdgeMap<double> Sweights(SG);
		SmartGraph::NodeMap<unsigned> Sindex(SG);
		vector<SmartGraph::Node> Snodes;

		makeKNNGraph(dv, clusters, packSize, *dcache, ccache, SG, Sweights, Sindex,
				Snodes);

		for (SmartGraph::EdgeIt e(SG); e != INVALID; ++e)
		{
			SmartGraph::Node u = SG.u(e);
			SmartGraph::Node v = SG.v(e);

			unsigned i = Sindex[u];
			unsigned j = Sindex[v];
			float dist = Sweights[e];

			distances.push_back(IntraClusterDist(i,j,dist));
		}

	}
}

