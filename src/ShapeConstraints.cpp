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
 * Excluder.cpp
 *
 *  Created on: Jun 11, 2012
 *      Author: dkoes
 */

#include "ShapeObj.h"
#include <string>
#include <cmath>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <ShapeConstraints.h>
#include "MappableOctTree.h"
using namespace OpenBabel;
using namespace Eigen;


const double ShapeConstraints::probeRadius = 1.4; //radius of water

//return a "good" transformation to put the query coordinate system into
//a centered coordinate system; use the ligand coordinates if available in query,
//otherwise use pharmacophore points
Eigen::Affine3d ShapeConstraints::computeTransform(Json::Value& root)
{
	vector<Vector3d> coords;


	if(root["ligand"].isString() && root["ligandFormat"].isString())
	{
		OBConversion conv;
		string lname = root["ligandFormat"].asString();
		conv.SetInFormat(OBConversion::FormatFromExt(lname.c_str()));
		OBMol lig;
		conv.ReadString(&lig, root["ligand"].asString());
		for (OBAtomIterator aitr = lig.BeginAtoms(); aitr != lig.EndAtoms(); ++aitr)
		{
			OBAtom* atom = *aitr;
			if(!atom->IsHydrogen())
			{
				Vector3d c(atom->x(), atom->y(), atom->z());
				coords.push_back(c);
			}
		}
	}

	if(coords.size() == 0) //fallback on pharmacophroe featuers
	{
		Json::Value jpoints = root["points"];
		for (unsigned i = 0, n = jpoints.size(); i < n; i++)
		{
			Json::Value jpnt = jpoints[i];
			if (jpnt.isMember("enabled") && !jpnt["enabled"].asBool())
				continue;
			string name = jpnt["name"].asString();
			if(name == "ExclusionSphere" || name == "InclusionSphere")
				continue; //not "real" pharmacophores

			double x = jpnt["x"].asDouble();
			double y = jpnt["y"].asDouble();
			double z = jpnt["z"].asDouble();

			coords.push_back(Vector3d(x,y,z));

		}
	}

	Vector3d translate = Vector3d::Zero();
	Matrix3d rotate = Matrix3d::Identity();

	ShapeObj::computeAndApplyNormalization(coords, translate, rotate);

	//grid must also be aligned to moments
	Affine3d ret = Affine3d::Identity();
	ret.rotate(rotate);
	ret.translate(translate);

	return ret;
}

void ShapeConstraints::makeGrid(MGrid& grid, OBMol& mol, const Affine3d& transform, double tolerance)
{
	//set adjustments based on tolerance
	double grow = 0.0;
	double shrink = 0.0;
	if(tolerance < 0) grow = -tolerance;
	else shrink = tolerance;

	grid.clear();
	//create vdw grid
	for (OBAtomIterator aitr = mol.BeginAtoms(); aitr != mol.EndAtoms(); ++aitr)
	{
		OBAtom* atom = *aitr;
		Vector3d c(atom->x(), atom->y(), atom->z());
		c = transform*c;
		grid.markXYZSphere(c.x(), c.y(), c.z(), etab.GetVdwRad(atom->GetAtomicNum()));
	}

	//make solvent accessible grid
	MGrid sagrid = grid;
	sagrid.grow(probeRadius+grow);
	grid.makeSurface(sagrid,probeRadius);
	grid.shrink(shrink);
}

//return true if the shape constraints actually specify something
bool ShapeConstraints::isMeaningful() const
{
	bool ingood = inclusiveKind != None;
	if(inclusiveKind == Spheres && inspheres.size() == 0)
		ingood = false;
	if(inclusiveKind == Shape && includeGrid.numSet() == 0)
		ingood = false;

	bool exgood = exclusiveKind != None;
	if(exclusiveKind == Spheres && exspheres.size() == 0)
		exgood = false;
	if(exclusiveKind == Shape && excludeGrid.numSet() == 0)
		exgood = false;

	return ingood || exgood;
}


//read exclusion sphere points from a json formatted stream
//and add to excluder
bool ShapeConstraints::readJSONExclusion(Json::Value& root)
{
	try
	{
		gridtransform = computeTransform(root);
		excludeGrid.clear();
		includeGrid.clear();
		ligandGrid.clear();

		//read setting for type of shape constraints
		Json::Value exselect = root["exselect"];
		Json::Value inselect = root["inselect"];

		//set what sort of shape constraints these are
		if(exselect.isString())
		{
			if(exselect.asString() == "points")
				exclusiveKind = Spheres;
			else if(exselect.asString() == "receptor")
				exclusiveKind = Shape;
			else
				exclusiveKind = None;
		}
		else
			exclusiveKind = None;

		if(inselect.isString())
		{
			if(inselect.asString() == "points")
				inclusiveKind = Spheres;
			else if(inselect.asString() == "ligand")
				inclusiveKind = Shape;
			else
				inclusiveKind = None;
		}
		else
			inclusiveKind = None;


		if(exclusiveKind == Spheres || inclusiveKind == Spheres)
		{ //scan through points
			Json::Value jpoints = root["points"];
			for (unsigned i = 0, n = jpoints.size(); i < n; i++)
			{
				Json::Value jpnt = jpoints[i];
				if (jpnt.isMember("enabled") && !jpnt["enabled"].asBool())
					continue;
				string name = jpnt["name"].asString();
				if(exclusiveKind == Spheres && name == "ExclusionSphere")
				{
					double radius = 0;
					if (jpnt.isMember("radius"))
						radius = jpnt["radius"].asDouble();
					if(radius > 0)
					{
						double x = jpnt["x"].asDouble();
						double y = jpnt["y"].asDouble();
						double z = jpnt["z"].asDouble();
						exspheres.push_back(Sphere(x,y,z,radius));

						//for grid, convert to grid space
						Vector3d c(x,y,z);
						c = gridtransform*c;
						excludeGrid.markXYZSphere(c.x(),c.y(),c.z(),radius);
					}
				}
				else if(inclusiveKind == Spheres && name == "InclusionSphere")
				{
					double radius = 0;
					if (jpnt.isMember("radius"))
						radius = jpnt["radius"].asDouble();
					if(radius > 0)
					{
						double x = jpnt["x"].asDouble();
						double y = jpnt["y"].asDouble();
						double z = jpnt["z"].asDouble();

						inspheres.push_back(Sphere(x,y,z,radius));
						Vector3d c(x,y,z);
						c = gridtransform*c;
						includeGrid.markXYZSphere(c.x(),c.y(),c.z(),radius);
					}
				}
			}
		}

		//check for full shape constraints
		if(exclusiveKind == Shape)
		{
			//get tolerance,
			double tolerance = 0.0;

			if(root["extolerance"].isNumeric())
			{
				tolerance = root["extolerance"].asDouble();
			}

			//parse receptor
			OBConversion conv;
			string rname = root["recname"].asString();
			conv.SetInFormat(OBConversion::FormatFromExt(rname.c_str()));
			//ignore bonds since we don't need them and openbabel likes to crash perceiving them
			conv.AddOption("b",OBConversion::INOPTIONS);

			OBMol rec;
			conv.ReadString(&rec, root["receptor"].asString());

			makeGrid(excludeGrid, rec, gridtransform, tolerance);
		}

		OBMol lig;
		if(root["ligand"].isString())
		{
			ligandGrid.clear();
			//parse ligand if available and create grid
			OBConversion conv;
			string lname = root["ligandFormat"].asString();
			conv.SetInFormat(OBConversion::FormatFromExt(lname.c_str()));
			OBMol lig;
			conv.ReadString(&lig, root["ligand"].asString());

			double resolution = ligandGrid.getResolution();
			double dimension = ligandGrid.getDimension();

			//transform to grid space
			for (OBAtomIterator aitr = lig.BeginAtoms(); aitr != lig.EndAtoms(); ++aitr)
			{
				OBAtom* atom = *aitr;
				Vector3d c(atom->x(), atom->y(), atom->z());
				c = gridtransform*c;
				atom->SetVector(c.x(), c.y(), c.z());
			}

			//compute grid using obanalytic for closest fidelty to database shapes
			ShapeObj::MolInfo minfo;
			ShapeObj obj(lig, Vector3d::Zero(), Matrix3d::Identity(), minfo, dimension, resolution);
			MappableOctTree *tree = MappableOctTree::create(dimension, resolution, obj);
			tree->makeGrid(ligandGrid, 0.5);
			delete tree;
		}
		else
		{
			ligandGrid = includeGrid; //let ligand be spheres if it doesn't exist
		}

		if(inclusiveKind == Shape)
		{
			//get tolerance,
			double tolerance = 0.0;
			if(root["intolerance"].isNumeric())
			{
				tolerance = root["intolerance"].asDouble();
			}

			includeGrid = ligandGrid;

			if(tolerance > 0)
				includeGrid.shrink(tolerance);
			else if(tolerance < 0)
				includeGrid.grow(-tolerance);
		}

	} catch (std::exception& e)
	{
		//poorly formated json
		cerr << "Parse " << e.what() << "\n";
		return false;
	}
	return true;
}

bool ShapeConstraints::isExcluded(PMol *mol, const RMSDResult& res) const
{
	using namespace Eigen;
	vector<Vector3f> coords;

	Transform<double, 3, Affine>  transform = Transform<double, 3, Affine>::Identity();
	Vector3f trans = res.translationVector();
	Matrix3f rot = res.rotationMatrix();
	transform.translate(trans.cast<double>());
	transform.rotate(rot.cast<double>());

	mol->getCoords(coords, transform);

	if(exclusiveKind == Shape)
	{
		for (unsigned c = 0, nc = coords.size(); c < nc; c++)
		{
			const Vector3d& pnt = gridtransform*coords[c].cast<double>();
			float x = pnt.x();
			float y = pnt.y();
			float z = pnt.z();
			if(excludeGrid.inGrid(x,y,z) && excludeGrid.test(x,y,z))
				return true;
		}
	}
	else if(exclusiveKind == Spheres)
	{
		//todo: accelerate with grid?  but it isn't perfect.. evaluate accuracy/performance tradeoff

		//if any coordinate overlaps with any exclusion sphere, no good
		for(unsigned i = 0, n = exspheres.size(); i < n; i++)
		{
			for (unsigned c = 0, nc = coords.size(); c < nc; c++)
			{
				const Vector3f& pnt = coords[c];
				if(exspheres[i].contains(pnt.x(),pnt.y(),pnt.z()))
					return true;
			}
		}
	}

	if(inclusiveKind == Shape)
	{
		//as long as a single heavy atom is included, we are good
		bool nogood = true;
		for (unsigned c = 0, nc = coords.size(); c < nc && nogood; c++)
		{
			const Vector3d& pnt = gridtransform*coords[c].cast<double>();
			float x = pnt.x();
			float y = pnt.y();
			float z = pnt.z();
			if(includeGrid.inGrid(x,y,z) && includeGrid.test(x,y,z))
				nogood = false;
		}
		if(nogood) //nothing overlapped
			return true;
	}
	else if(inclusiveKind == Spheres)
	{
		//at least one atom must overlap inclusion sphere
		for(unsigned i = 0, n = inspheres.size(); i < n; i++)
		{
			unsigned c = 0, nc = coords.size();
			for ( ; c < nc; c++)
			{
				const Vector3f& pnt = coords[c];
				if(inspheres[i].contains(pnt.x(),pnt.y(),pnt.z()))
					break;
			}
			if(c == nc) //nothing overlapped
				return true;
		}
	}
	return false;
}

void ShapeConstraints::addToJSON(Json::Value& root) const
{
	Json::Value jpoints = root["points"];
	unsigned start = jpoints.size();
	for(unsigned i = 0, n = exspheres.size(); i < n; i++)
	{
		Json::Value& pt = root["points"][start+i];
		pt["name"] = "ExclusionSphere";
		pt["x"] = exspheres[i].x;
		pt["y"] = exspheres[i].y;
		pt["z"] = exspheres[i].z;
		pt["radius"] = sqrt(exspheres[i].rSq);
	}
	start = jpoints.size();
	for(unsigned i = 0, n = inspheres.size(); i < n; i++)
	{
		Json::Value& pt = root["points"][start+i];
		pt["name"] = "InclusionSphere";
		pt["x"] = inspheres[i].x;
		pt["y"] = inspheres[i].y;
		pt["z"] = inspheres[i].z;
		pt["radius"] = sqrt(inspheres[i].rSq);
	}
}

//truncate floats to reduce the number of digits
static double Round2(double num)
{
    return round(num * 100.0)/100.0;
}

//set json formatted mesh of provided grid
void ShapeConstraints::getMesh(MGrid& grid, Json::Value& mesh)
{
	Json::Value& verts = mesh["vertexArr"] =  Json::arrayValue;
	Json::Value& norms =mesh["normalArr"] =  Json::arrayValue;
	Json::Value& jfaces = mesh["faceArr"] =  Json::arrayValue;

	vector<Vector3f> vertices;
	vector<Vector3f> normals;
	vector<int> faces;

	//create mesh from mgrid
	grid.makeMesh(vertices, normals, faces);

	//copy face indices
	for(unsigned i = 0, n = faces.size(); i < n; i++)
	{
		jfaces[i] = faces[i];
	}

	Affine3d trans = gridtransform.inverse();

	//have to transform vertices and normals
	for(unsigned i = 0, n = vertices.size(); i < n; i++)
	{
		Vector3d v = trans*vertices[i].cast<double>();
		verts[i]["x"] = Round2(v.x());
		verts[i]["y"] = Round2(v.y());
		verts[i]["z"] = Round2(v.z());
	}

	return;
	for(unsigned i = 0, n = normals.size(); i < n; i++)
	{
		Vector3d norm = normals[i].cast<double>();
		norm.normalize();
		Vector3d v = trans.linear()*norm;
		norms[i]["x"] = Round2(v.x());
		norms[i]["y"] = Round2(v.y());
		norms[i]["z"] = Round2(v.z());
	}
}

//set json formated mesh of exclusive grid
void ShapeConstraints::getExclusiveMesh(Json::Value& mesh)
{
	getMesh(excludeGrid, mesh);
}

void ShapeConstraints::getInclusiveMesh(Json::Value& mesh)
{
	getMesh(includeGrid, mesh);
}

void ShapeConstraints::computeInteractionPoints(OBMol& ligand, OBMol& receptor, vector<Vector3d>& points)
{
	points.clear();
	//compute steric points
	ShapeObj::MolInfo minfo;
	ShapeObj obj(ligand, Vector3d::Zero(), Matrix3d::Identity(), minfo, PHARMIT_DIMENSION, PHARMIT_RESOLUTION);
	obj.computeInteractionPoints(receptor, points);
}
