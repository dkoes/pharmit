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
bool Excluder::addJSONPoints(Json::Value& root)
{
	try
	{
		Json::Value jpoints = root["points"];
		for (unsigned i = 0, n = jpoints.size(); i < n; i++)
		{
			Json::Value jpnt = jpoints[i];
			if (jpnt.isMember("enabled") && !jpnt["enabled"].asBool())
				continue;
			string name = jpnt["name"].asString();
			if(name == "ExclusionSphere")
			{
				double radius = 0;
				if (jpnt.isMember("radius"))
					radius = jpnt["radius"].asDouble();
				if(radius > 0)
				{
					double x = jpnt["x"].asDouble();
					double y = jpnt["y"].asDouble();
					double z = jpnt["z"].asDouble();
					spheres.push_back(ExclusionSphere(x,y,z,radius));
				}
			}
		}

	} catch (std::exception& e)
	{
		//poorly formated json
		cerr << "Parse " << e.what() << "\n";
		return false;
	}
	return true;
}

bool Excluder::isExcluded(const FloatCoord& pnt) const
{
	for(unsigned i = 0, n = spheres.size(); i < n; i++)
	{
		if(spheres[i].contains(pnt.x,pnt.y,pnt.z))
			return true;
	}
	return false;
}

void Excluder::addToJSON(Json::Value& root) const
{
	Json::Value jpoints = root["points"];
	unsigned start = jpoints.size();
	for(unsigned i = 0, n = spheres.size(); i < n; i++)
	{
		Json::Value& pt = root["points"][start+i];
		pt["name"] = "ExclusionSphere";
		pt["x"] = spheres[i].x;
		pt["y"] = spheres[i].y;
		pt["z"] = spheres[i].z;
		pt["radius"] = sqrt(spheres[i].rSq);
	}
}
