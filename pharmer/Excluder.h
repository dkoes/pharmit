/*
 * Excluder.h
 *
 *  Created on: Jun 11, 2012
 *      Author: dkoes
 *
 *      This class is used to define an excluded space and is used to
 *      check to see if any points fall within this space.
 */

#ifndef EXCLUDER_H_
#define EXCLUDER_H_

#include "FloatCoord.h"
#include <json/json.h>
#include <vector>
using namespace std;

class Excluder
{
	struct ExclusionSphere
	{
		float x, y, z;
		float rSq;

		bool contains(float a, float b, float c) const
		{
			float d1 = a-x;
			d1 *= d1;
			float d2 = b-y;
			d2 *= d2;
			float d3 = c-z;
			d3 *= d3;
			float distSq = d1 + d2 + d3;
			return distSq <= rSq;
		}

		ExclusionSphere(): x(0),y(0),z(0),rSq(0) {}
		ExclusionSphere(float _x, float _y, float _z, float _r): x(_x), y(_y), z(_z), rSq(_r*_r) {}
	};

	vector<ExclusionSphere> spheres;

public:
	Excluder() {}
	~Excluder() {}

	//read exclusion spheres from json
	bool addJSONPoints(Json::Value& root);

	//write exclusion information to root
	void addToJSON(Json::Value& root) const;

	//return true if in the exclusion zone
	bool isExcluded(const FloatCoord& pnt) const;

	bool isDefined() const { return spheres.size() > 0; }

	void clear() { spheres.clear(); }

	void addExclusionSphere(float x, float y, float z, float r)
	{
		spheres.push_back(ExclusionSphere(x,y,z,r));
	}
};

#endif /* EXCLUDER_H_ */
