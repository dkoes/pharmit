/*
 * Excluder.cpp
 *
 *  Created on: Jun 11, 2012
 *      Author: dkoes
 */

#include "Excluder.h"
#include "ShapeObj.h"
#include <string>
#include <cmath>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
using namespace OpenBabel;
using namespace Eigen;


const double Excluder::probeRadius = 1.4; //radius of water

//return a "good" transformation to put the query coordinate system into
//a centered coordinate system; use the ligand coordinates if available in query,
//otherwise use pharmacophore points
Eigen::Affine3d Excluder::computeTransform(Json::Value& root)
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
			Vector3d c(atom->x(), atom->y(), atom->z());
			coords.push_back(c);
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

	Affine3d ret = Affine3d::Identity();
	ret.rotate(rotate);
	ret.translate(translate);

	return ret;
}

void Excluder::makeGrid(MGrid& grid, OBMol& mol, const Affine3d& transform, double tolerance)
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


//read exclusion sphere points from a json formatted stream
//and add to excluder
bool Excluder::readJSONExclusion(Json::Value& root)
{
	try
	{
		bool getexpoints = true, getinpoints = true;

		//read setting for type of shape constraints
		Json::Value exselect = root["exselect"];
		if(exselect.isString() && exselect.asString() != "points")
			getexpoints = false;

		Json::Value inselect = root["inselect"];
		if(inselect.isString() && inselect.asString() != "points")
			getinpoints = false;

		if(getinpoints || getexpoints)
		{ //scan through points
			Json::Value jpoints = root["points"];
			for (unsigned i = 0, n = jpoints.size(); i < n; i++)
			{
				Json::Value jpnt = jpoints[i];
				if (jpnt.isMember("enabled") && !jpnt["enabled"].asBool())
					continue;
				string name = jpnt["name"].asString();
				if(getexpoints && name == "ExclusionSphere")
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
					}
				}
				else if(getinpoints && name == "InclusionSphere")
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
					}
				}
			}
		}

		//check for full shape constraints
		gridtransform = computeTransform(root);
		excludeGrid.clear();
		includeGrid.clear();
		if(exselect.isString() && exselect.asString() == "receptor" && root["receptor"].isString())
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
			OBMol rec;
			conv.ReadString(&rec, root["receptor"].asString());

			makeGrid(excludeGrid, rec, gridtransform, tolerance);
		}


	} catch (std::exception& e)
	{
		//poorly formated json
		cerr << "Parse " << e.what() << "\n";
		return false;
	}
	return true;
}

bool Excluder::isExcluded(PMol *mol, const RMSDResult& res) const
{
	using namespace Eigen;
	vector<Vector3f> coords;

	Transform<double, 3, Affine>  transform = Transform<double, 3, Affine>::Identity();
	Vector3f trans = res.translationVector();
	Matrix3f rot = res.rotationMatrix();
	transform.translate(trans.cast<double>());
	transform.rotate(rot.cast<double>());

	mol->getCoords(coords, transform);

	if(excludeGrid.numSet() > 0)
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

	return false;
}

void Excluder::addToJSON(Json::Value& root) const
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

//set json formated mesh of exclusive grid
void Excluder::getExclusiveMesh(Json::Value& mesh)
{
	Json::Value& verts = mesh["vertexArr"] =  Json::arrayValue;
	Json::Value& norms =mesh["normalArr"] =  Json::arrayValue;
	Json::Value& jfaces = mesh["faceArr"] =  Json::arrayValue;

	vector<Vector3f> vertices;
	vector<Vector3f> normals;
	vector<int> faces;

	//create mesh from mgrid
	excludeGrid.makeMesh(vertices, normals, faces);

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
		verts[i]["x"] = v.x();
		verts[i]["y"] = v.y();
		verts[i]["z"] = v.z();
	}

	for(unsigned i = 0, n = normals.size(); i < n; i++)
	{
		Vector3d v = trans.linear()*vertices[i].cast<double>();
		v.normalize();
		norms[i]["x"] = v.x();
		norms[i]["y"] = v.y();
		norms[i]["z"] = v.z();
	}
}

