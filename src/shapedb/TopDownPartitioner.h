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
 * TopDownPartitioner.h
 *
 *  Created on: Oct 18, 2011
 *      Author: dkoes
 */

#ifndef TOPDOWNPARTITIONER_H_
#define TOPDOWNPARTITIONER_H_

#include "GSSTypes.h"
#include "DataViewers.h"


//abstract class for a top down partitioner, is expected to run in linear time
//can maintain state as the partitions are refined
class TopDownPartitioner
{
protected:
	const DataViewer *data;
	vector<unsigned> indices; //indexes into data
public:
	TopDownPartitioner(): data(NULL) {}
	TopDownPartitioner(const DataViewer *d): data(d) {}
	virtual ~TopDownPartitioner() {}
	virtual TopDownPartitioner* create(const DataViewer* dv) const = 0;
	virtual void partition(vector<TopDownPartitioner*>& parts) = 0;

	virtual unsigned size() const { return indices.size(); }
	virtual const DataViewer * getData() const { return data; }
	virtual void extractIndicies(vector<unsigned>& ind) { swap(ind, indices); }

	void initFromData()
	{
		indices.clear(); indices.resize(data->size());
		for(unsigned i = 0, n = data->size(); i < n; i++)
			indices[i] = i;
	}
};

#endif /* TOPDOWNPARTITIONER_H_ */
