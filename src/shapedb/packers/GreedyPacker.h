/*
 * GreedyPacker.h
 *
 *  Created on: Nov 2, 2011
 *      Author: dkoes
 *
 *      A greedy packer that pulls cluster distances off a priority queue.
 *      The two closest clusters are merged, their distances are recomputed,
 *      and then put back on the Q.
 */

#ifndef GREEDYPACKER_H_
#define GREEDYPACKER_H_

#include "Packer.h"

class GreedyPacker: public Packer
{
public:
	GreedyPacker(unsigned ps, ClusterDistance metric=AverageLink, unsigned k = 0, unsigned s = 0): Packer(ps, metric, k, s) {}
	~GreedyPacker() {}

	virtual void pack(const DataViewer* dv, vector<Cluster>& clusters) const;

};

#endif /* GREEDYPACKER_H_ */
