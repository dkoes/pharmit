/*
 * Packer.h
 *
 *  Created on: Oct 18, 2011
 *      Author: dkoes
 */

#ifndef PACKER_H_
#define PACKER_H_

#include "DataViewers.h"
#include "ShapeDistance.h"
#include <boost/multi_array.hpp>
#include <lemon/smart_graph.h>

//distances between i and j
struct IntraClusterDist
{
	unsigned i;
	unsigned j;
	float dist;

	IntraClusterDist() :
			i(0), j(0), dist(HUGE_VAL)
	{
	}

	IntraClusterDist(unsigned I, unsigned J, float d) :
			i(I), j(J), dist(d)
	{
	}

	bool operator<(const IntraClusterDist& rhs) const
	{
		return dist < rhs.dist;
	}
};

//abstract class for a bottom up packer, can run in quadratic time
class Packer
{
public:
	enum ClusterDistance
	{
		CompleteLink, AverageLink, SingleLink, TotalLink, NotApplicable
	};
protected:
	unsigned packSize;
	ClusterDistance distMetric;

	//caches intercluster distances
	class ClusterCache
	{
		struct JDist
		{
			unsigned j;
			float dist;

			JDist() :
					j(0), dist(0)
			{
			}
			JDist(unsigned J, float D) :
					j(J), dist(D)
			{
			}

			bool operator<(const JDist& rhs) const
			{
				return j < rhs.j;
			}

			bool operator==(const JDist& rhs) const
			{
				return j == rhs.j;
			}
		};

		vector<vector<JDist> > cache;

	public:
		//indices are assumed to be consecutive and sorted
		ClusterCache(unsigned n) :
				cache(n)
		{

		}

		bool isCached(unsigned i, unsigned j, float& ret)
		{
			if (j > i) //normalize order - store only once
				swap(i, j);
			else if (i == j)
			{
				ret = 0;
				return true;
			}

			if (i >= cache.size())
				return false;

			vector<JDist>& row = cache[i];
			JDist val(j, 0);
			vector<JDist>::iterator pos = lower_bound(row.begin(), row.end(),
					val);

			if (pos != row.end() && pos->j == j)
			{
				//already exists
				ret = pos->dist;
				return true;
			}
			else
			{
				return false;
			}
		}

		void set(unsigned i, unsigned j, float v)
		{
			if (j > i) //normalize order - store only once
				swap(i, j);
			else if (i == j)
				return;

			if (i >= cache.size())
				cache.resize(i + 1);

			vector<JDist>& row = cache[i];
			JDist val(j, v);
			vector<JDist>::iterator pos = lower_bound(row.begin(), row.end(),
					val);

			if (pos != row.end() && pos->j == j)
			{
				//already exists
				pos->dist = v;
			}
			else
			{
				row.insert(pos, val);
			}
		}
	};

	class DCache //interface for caching distance information between shapes
	{
	public:
		virtual float get(unsigned i, unsigned j) = 0;
		virtual ~DCache() {}

	};

	//stores and precomputes everything
	class FullCache: public DCache
	{
		boost::multi_array<float, 2> cache;
		const DataViewer* dv;

	public:
		//indices are assumed to be consecutive and sorted
		FullCache(const DataViewer* D) :
				cache(boost::extents[D->size()][D->size()]), dv(D)
		{
			unsigned N = D->size();
			for (unsigned i = 0; i < N; i++)
			{
				for (unsigned j = 0; j < i; j++)
				{
					float dist = shapeDistance(D->getMIV(i), D->getMSV(i),
							D->getMIV(j), D->getMSV(j));
					cache[i][j] = dist;
					cache[j][i] = dist;
				}
				cache[i][i] = 0;
			}
		}

		float get(unsigned i, unsigned j)
		{
			return cache[i][j];
		}

		unsigned size() const
		{
			return cache.size() * cache.size();
		}
	};

	//does not precompute, uses much less memory
	class OnDemandCache: public DCache
	{
		struct JDist
		{
			unsigned j;
			float dist;

			JDist() :
					j(0), dist(0)
			{
			}
			JDist(unsigned J, float D) :
					j(J), dist(D)
			{
			}

			bool operator<(const JDist& rhs) const
			{
				return j < rhs.j;
			}

			bool operator==(const JDist& rhs) const
			{
				return j == rhs.j;
			}
		};

		vector<vector<JDist> > cache;
		const DataViewer* dv;

	public:
		//indices are assumed to be consecutive and sorted
		OnDemandCache(const DataViewer* D) :
				cache(D->size()), dv(D)
		{

		}

		float get(unsigned i, unsigned j)
		{
			if (j > i) //normalize order - store only once
				swap(i, j);
			else if (i == j)
				return 0;

			vector<JDist>& row = cache[i];
			JDist val(j, 0);
			vector<JDist>::iterator pos = lower_bound(row.begin(), row.end(),
					val);

			if (pos != row.end() && pos->j == j)
			{
				//already exists
				return pos->dist;
			}
			else
			{
				//must compute and insert
				val.dist = shapeDistance(dv->getMIV(i), dv->getMSV(i),
						dv->getMIV(j), dv->getMSV(j));
				row.insert(pos, val);
				return val.dist;
			}
		}

		unsigned size() const
		{
			unsigned ret = 0;
			for (unsigned i = 0, n = cache.size(); i < n; i++)
			{
				ret += cache[i].size();
			}
			return ret;
		}
	};

	struct IndDist
	{
		unsigned j;
		bool unprocessed;
		double dist;

		IndDist() :
				j(0), unprocessed(true), dist(HUGE_VAL)
		{
		}
		IndDist(unsigned J, double D) :
				j(J), unprocessed(true), dist(D)
		{
		}

		bool operator<(const IndDist& rhs) const
		{
			return dist < rhs.dist;
		}
	};

	struct KNNSlice
	{
		vector<IndDist> neighbors;

		bool update(const IndDist& item, unsigned k); //insert in sorted order, make K items, no limit if 0
		void getOldNew(vector<unsigned>& o, vector<unsigned>& n); //processed and unprocessed; once extracted becomes processed

		void getCommonNeighbors(const vector<unsigned>& rev,
				vector<unsigned>& neighs) const;
	};

	void initialKNNSample(const DataViewer *D, vector<Cluster>& clusters,
			unsigned maxSz, DCache& dcache, ClusterCache& ccache,
			vector<KNNSlice>& V, float& maxdist) const;

	float makeKNNGraph(const DataViewer *D, vector<Cluster>& clusters,
			unsigned maxSz, DCache& dcache, ClusterCache& ccache,
			lemon::SmartGraph& G, lemon::SmartGraph::EdgeMap<double>& E,
			lemon::SmartGraph::NodeMap<unsigned>& Sindex,
			vector<lemon::SmartGraph::Node>& nodes) const;

	float clusterDistance(const DataViewer* D, const Cluster& a,
			const Cluster& b, DCache& dcache) const;

	void computeDistanceVector(const DataViewer* D, vector<Cluster>& clusters,
			vector<IntraClusterDist>& distances, DCache *& dcache) const;

	unsigned K; //for knn, if zero use full
	unsigned S; //number of sentinals for knn building, if zero, do random
public:

	Packer(unsigned ps, ClusterDistance metric, unsigned k = 0, unsigned s = 0) :
			packSize(ps), distMetric(metric), K(k), S(s)
	{
	}
	virtual ~Packer()
	{
	}
	//assumes dv has been sliced into sequential indices
	virtual void pack(const DataViewer* dv,
			vector<Cluster>& clusters) const = 0;

	unsigned getPack() const
	{
		return packSize;
	}
};

#endif /* PACKER_H_ */
