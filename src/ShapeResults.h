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
 * ShapeResults.h
 *
 *  Created on: Aug 10, 2015
 *      Author: dkoes
 *
 * Implements shapedb results interface for putting hits into a corresponder queue
 */

#ifndef SHAPERESULTS_H_
#define SHAPERESULTS_H_

#include "Results.h"
#include "MTQueue.h"
#include "cors.h"
#include "ShapeConstraints.h"
#include "params.h"

class ShapeResults: public Results
{
	std::shared_ptr<PharmerDatabaseSearcher> dbptr;
	MTQueue<CorrespondenceResult*>& resultQ;
	CorAllocator& alloc;
	const QueryParameters& qparams;
	vector<PharmaPoint> points; //pharmacophore query points after alignment to grid

	unsigned db;
	unsigned numdb;
	RMSDResult defaultR; //all molecules have same transformation
	bool& stop;

public:
	ShapeResults(std::shared_ptr<PharmerDatabaseSearcher>& dptr, const vector<PharmaPoint>& querypoints,
			MTQueue<CorrespondenceResult*> & Q, CorAllocator& ca,
			const QueryParameters& qp, const ShapeConstraints& cons, unsigned whichdb, unsigned totaldb, bool& stopEarly);
	virtual ~ShapeResults() {}

	virtual void clear() {} //meaningless
	virtual void add(const char *data, double score);

	virtual void reserve(unsigned n) {}

	virtual unsigned size() const;

	virtual bool stopEarly() const { return stop; }
};

#endif /* SHAPERESULTS_H_ */
