/*
 * Excluder.cpp
 *
 *  Created on: Jun 11, 2012
 *      Author: dkoes
 */

#include "Excluder.h"
#include <string>
#include <cmath>

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
