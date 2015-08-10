/*
 * ShapeConstraints.h
 *
 *  Created on: Jun 11, 2012
 *      Author: dkoes
 *
 *      This class is used to define an excluded space and is used to
 *      check to see if any points fall within this space.
 */

#ifndef PHARMITSERVER_SHAPECONSTRAINTS_H_
#define PHARMITSERVER_SHAPECONSTRAINTS_H_

#include "FloatCoord.h"
#include "pharmarec.h"
#include "PMol.h"
#include "MGrid.h"
#include <json/json.h>
#include <vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include "ShapeObj.h"

using namespace std;

class ShapeConstraints
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

	MGrid excludeGrid;
	MGrid includeGrid;
	MGrid ligandGrid;
	Eigen::Affine3d gridtransform;//transformation to be applied to points to move into grid space

	Eigen::Affine3d computeTransform(Json::Value& root);
	void makeGrid(MGrid& grid, OpenBabel::OBMol& mol, const Eigen::Affine3d& transform, double tolerance);

	static const double probeRadius;

	void getMesh(MGrid& grid, Json::Value& mesh);

	enum Kind {None, Shape, Spheres};

	Kind inclusiveKind;
	Kind exclusiveKind;
public:
	ShapeConstraints(): excludeGrid(PHARMIT_DIMENSION,PHARMIT_RESOLUTION),
	includeGrid(PHARMIT_DIMENSION,PHARMIT_RESOLUTION),
	ligandGrid(PHARMIT_DIMENSION,PHARMIT_RESOLUTION),
	inclusiveKind(None),exclusiveKind(None) {}
	~ShapeConstraints() {}

	//read steric constraints from json
	bool readJSONExclusion(Json::Value& root);

	//write exclusion information to root
	void addToJSON(Json::Value& root) const;

	//return true if coordinates are excluded by spatial constraints
	//provide the molecule and any transformation to the coordinates
	bool isExcluded(PMol *mol, const RMSDResult& res) const;

	bool isDefined() const { return inclusiveKind != None || exclusiveKind != None; }

	void clear() { exspheres.clear(); inspheres.clear(); }

	void addExclusionSphere(float x, float y, float z, float r)
	{
		exspheres.push_back(Sphere(x,y,z,r));
	}

	void addInclusionSphere(float x, float y, float z, float r)
	{
		inspheres.push_back(Sphere(x,y,z,r));
	}

	void getExclusiveMesh(Json::Value& mesh);
	void getInclusiveMesh(Json::Value& mesh);

	const MGrid& getExclusiveGrid() const { return excludeGrid; }
	const MGrid& getInclusiveGrid() const { return includeGrid; }
	const MGrid& getLigandGrid() const { return ligandGrid; }

	Eigen::Affine3d getGridTransform() const { return gridtransform; }
	static void computeInteractionPoints(OpenBabel::OBMol& ligand, OpenBabel::OBMol& receptor, vector<Eigen::Vector3d>& points);
};

#endif /* PHARMITSERVER_SHAPECONSTRAINTS_H_ */
