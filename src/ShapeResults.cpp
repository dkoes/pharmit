/*
 * ShapeResults.cpp
 *
 *  Created on: Aug 10, 2015
 *      Author: dkoes
 */

#include <ShapeResults.h>
#include "ShapeObj.h"

ShapeResults::ShapeResults(boost::shared_ptr<PharmerDatabaseSearcher>& dptr,
		const vector<PharmaPoint>& querypoints,
		MTQueue<CorrespondenceResult*>& Q, CorAllocator& ca,
		const QueryParameters& qp,
		const ShapeConstraints& cons, unsigned whichdb, unsigned totaldb, bool& stopEarly) :
		dbptr(dptr), resultQ(Q), alloc(ca), qparams(qp), db(whichdb), numdb(totaldb), stop(stopEarly)
{
	Affine3d transform = cons.getGridTransform();
	Affine3d itransform = transform.inverse();
	defaultR = RMSDResult(0, itransform.translation(), itransform.rotation());

	//create query points transformed to grid space
	points = querypoints;

	for(unsigned i = 0, n = points.size(); i < n; i++)
	{
		PharmaPoint& pt = points[i];
		Vector3d xyz(pt.x, pt.y, pt.z);
		xyz = transform*xyz;
		pt.x = xyz.x();
		pt.y = xyz.y();
		pt.z = xyz.z();

		//also the vectors!
		for(unsigned v = 0, nv = pt.vecs.size(); v < nv; v++)
		{
			vector3& orig = pt.vecs[v];
			Vector3d vec(orig.x(), orig.y(), orig.z());
			vec = transform.linear()*vec;
			orig.Set(vec.x(), vec.y(), vec.z());
		}
	}
}

//add to queue
void ShapeResults::add(const char *data, double score)
{
	ShapeObj::MolInfo minfo;
	memcpy(&minfo, data, sizeof(ShapeObj::MolInfo));

	unsigned long loc = minfo.molPos;
	unsigned mid = ThreePointData::unpackMolID(loc);

	//filter out unsavory characters
	//check against query params
	if(minfo.nrot < qparams.minRot)
		return;
	if(minfo.nrot > qparams.maxRot)
		return;

	if(minfo.weight < qparams.reducedMinWeight)
		return;
	if(minfo.weight > qparams.reducedMaxWeight)
		return;

	for(unsigned i = 0, n = qparams.propfilters.size(); i < n; i++)
	{
		const PropFilter& prop = qparams.propfilters[i];
		//look up value for property
		double val = dbptr->getMolProp(prop.kind, mid);
		if(val < prop.min || val > prop.max)
		{
			return;
		}
	}

	if(points.size() > 0)
	{
		//pharmacophore filter
		if(!dbptr->alignedPharmasMatch(minfo.pharmPos, points))
			return;
	}

	CorrespondenceResult *res = alloc.newCorResult();

	//fill in cr with minfo
	res->location = loc * numdb + db;
	res->molid = mid * numdb + db;
	res->dbid = mid;
	res->val = 1.0-score;
	res->weight = ThreePointData::unreduceWeight(minfo.weight);
	res->nRBnds = minfo.nrot;
	res->rmsd = defaultR;
	res->rmsd.setValue(1.0-score); //overloading

	resultQ.push(res);
}

unsigned ShapeResults::size() const //this isn't really meaningful
{
	return resultQ.size();
}
