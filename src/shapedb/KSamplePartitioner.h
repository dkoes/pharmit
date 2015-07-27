/*
 * KSamplePartitioner.h
 *
 *  Created on: Oct 17, 2011
 *      Author: dkoes
 */

#ifndef KSAMPLEPARTITIONER_H_
#define KSAMPLEPARTITIONER_H_

#include "TopDownPartitioner.h"
#include "MappableOctTree.h"

class KSamplePartitioner: public TopDownPartitioner
{
public:
	enum CenterChoice {AveCenter, MinMaxCenter};
private:
	unsigned kcenters; //number of centers
	unsigned ksamples; //sample ksamples*kcenters
	CenterChoice centerFind; //algorithm for choosing center of a sampled cluster
	unsigned stopPartitionSize; // size where we stop partitioning
	void kCluster(const vector<unsigned>& indices, vector< vector<unsigned> >& clusters);

	void getCenter(const vector<unsigned>& cluster, const MappableOctTree *& MIV, const MappableOctTree *& MSV) const;

	unsigned fitKCenterToSize(unsigned  n) const;
	virtual TopDownPartitioner* create(const DataViewer* dv, vector<unsigned>& ind) const;

public:
	KSamplePartitioner(unsigned kc, unsigned ks, CenterChoice ch, unsigned stop): kcenters(kc), ksamples(ks), stopPartitionSize(stop) {}
	KSamplePartitioner(const DataViewer *dv, unsigned kc, unsigned ks, CenterChoice ch, unsigned stop): TopDownPartitioner(dv), kcenters(kc), ksamples(ks), centerFind(ch), stopPartitionSize(stop) {}
	~KSamplePartitioner() {}

	virtual TopDownPartitioner* create(const DataViewer* dv) const;
	virtual void partition(vector<TopDownPartitioner*>& parts);
};

#endif /* KSAMPLEPARTITIONER_H_ */
