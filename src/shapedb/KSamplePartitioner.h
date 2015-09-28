/*
Pharmit
Copyright (c) David Ryan Koes, University of Pittsburgh and contributors.
All rights reserved.

Pharmit is licensed under both the BSD 3-clause license and the GNU
Public License version 2. Any use of the code that retains its reliance
on the GPL-licensed OpenBabel library is subject to the terms of the GPL2.

Use of the Pharmit code independently of OpenBabel (or any other
GPL2 licensed software) may choose between the BSD or GPL licenses.

See the LICENSE file provided with the distribution for more information.

*/
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
	KSamplePartitioner(unsigned kc=8, unsigned ks=5, CenterChoice ch=AveCenter, unsigned stop=32768): kcenters(kc), ksamples(ks), stopPartitionSize(stop) {}
	KSamplePartitioner(const DataViewer *dv, unsigned kc, unsigned ks, CenterChoice ch, unsigned stop): TopDownPartitioner(dv), kcenters(kc), ksamples(ks), centerFind(ch), stopPartitionSize(stop) {}
	~KSamplePartitioner() {}

	virtual TopDownPartitioner* create(const DataViewer* dv) const;
	virtual void partition(vector<TopDownPartitioner*>& parts);
};

#endif /* KSAMPLEPARTITIONER_H_ */
