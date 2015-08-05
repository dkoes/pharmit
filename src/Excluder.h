/*
 * Excluder.h
 *
 *  Created on: Jun 11, 2012
 *      Author: dkoes
 *
 *      This class is used to define an excluded space and is used to
 *      check to see if any points fall within this space.
 */

#ifndef PHARMITSERVER_EXCLUDER_H_
#define PHARMITSERVER_EXCLUDER_H_

#include "FloatCoord.h"
#include "PMol.h"
#include <json/json.h>
#include <vector>
using namespace std;

class Excluder
{
	struct Sphere
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

		Sphere(): x(0),y(0),z(0),rSq(0) {}
		Sphere(float _x, float _y, float _z, float _r): x(_x), y(_y), z(_z), rSq(_r*_r) {}
	};

	vector<Sphere> exspheres; //exclusion spheres - can't overlap any
	vector<Sphere> inspheres; //inclusion spheres - must overlap all

public:
	Excluder() {}
	~Excluder() {}

	//read steric constraints from json
	bool readJSONExclusion(Json::Value& root);

	//write exclusion information to root
	void addToJSON(Json::Value& root) const;

	//return true if coordinates are excluded by spatial constraints
	//provide the molecule and any transformation to the coordinates
	bool isExcluded(PMol *mol, const RMSDResult& res) const;

	bool isDefined() const { return exspheres.size() > 0 || inspheres.size() > 0; }

	void clear() { exspheres.clear(); inspheres.clear(); }

	void addExclusionSphere(float x, float y, float z, float r)
	{
		exspheres.push_back(Sphere(x,y,z,r));
	}

	void addInclusionSphere(float x, float y, float z, float r)
	{
		inspheres.push_back(Sphere(x,y,z,r));
	}
};

#endif /* PHARMITSERVER_EXCLUDER_H_ */
