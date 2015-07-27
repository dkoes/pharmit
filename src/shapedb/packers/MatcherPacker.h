/*
 * MatcherPacker.h
 *
 *  Created on: Nov 2, 2011
 *      Author: dkoes
 *
 *      This packer is similar to the FullMerge packer, but computes the optimal
 *      edge matching at each merge.
 */

#ifndef MATCHERPACKER_H_
#define MATCHERPACKER_H_

#include "Packer.h"


class MatcherPacker: public Packer
{
	bool fullMergeClusters(const DataViewer* D, vector<Cluster>& clusters,
			unsigned maxSz, DCache& dcache) const;

	bool knnMergeClusters(const DataViewer* D, vector<Cluster>& clusters,
			unsigned maxSz, DCache& dcache) const;


public:
	MatcherPacker(unsigned ps, unsigned k = 0, unsigned s = 0,
			ClusterDistance metric = AverageLink) :
			Packer(ps, metric, k, s)
	{
	}
	~MatcherPacker()
	{
	}

	virtual void pack(const DataViewer* dv, vector<Cluster>& clusters) const;

};

#endif /* MATCHERPACKER_H_ */
