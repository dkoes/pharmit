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
	boost::shared_ptr<PharmerDatabaseSearcher> dbptr;
	MTQueue<CorrespondenceResult*>& resultQ;
	CorAllocator& alloc;
	const QueryParameters& qparams;
	vector<PharmaPoint> points; //pharmacophore query points after alignment to grid

	unsigned db;
	unsigned numdb;
	RMSDResult defaultR; //all molecules have same transformation
	bool& stop;

public:
	ShapeResults(boost::shared_ptr<PharmerDatabaseSearcher>& dptr, const vector<PharmaPoint>& querypoints,
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
