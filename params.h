/*
Pharmer: Efficient and Exact 3D Pharmacophore Search
Copyright (C) 2011  David Ryan Koes and the University of Pittsburgh

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

/*
 * params.h
 *
 *  Created on: Oct 19, 2010
 *      Author: dkoes
 *
 *      Simple parameter storing classes.
 */

#ifndef PARAMS_H_
#define PARAMS_H_
#include <iostream>
#include <boost/tuple/tuple.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/unordered_set.hpp>
#include <vector>
#include <json/json.h>

#include "ThreePointData.h"

using namespace std;

namespace SortType {
enum SortType {Undefined, RMSD, MolWeight, NRBnds};
}
typedef SortType::SortType SortTyp;

// Isolates parameters of a query.
struct QueryParameters
{
	double maxRMSD;
	unsigned reduceConfs;
	unsigned orientationsPerConf;
	unsigned maxHits; //total returned hits
	SortTyp sort; //for determining how to truncate
	//add support for more later

	double minWeight;
	double maxWeight;

	unsigned reducedMinWeight;
	unsigned reducedMaxWeight;
	unsigned minRot;
	unsigned maxRot;

	QueryParameters() :
		maxRMSD(HUGE_VAL), reduceConfs(UINT_MAX), orientationsPerConf(UINT_MAX), maxHits(UINT_MAX),
		sort(SortType::Undefined), minWeight(0), maxWeight(UINT_MAX), reducedMinWeight(0), reducedMaxWeight(UINT_MAX), minRot(0), maxRot(UINT_MAX)
	{

	}

	//extract parameters from json
	QueryParameters(Json::Value& data) :
		maxRMSD(HUGE_VAL), reduceConfs(UINT_MAX), orientationsPerConf(UINT_MAX), maxHits(UINT_MAX),
		minWeight(0), maxWeight(HUGE_VAL),reducedMinWeight(0), reducedMaxWeight(UINT_MAX), minRot(0), maxRot(UINT_MAX)
	{
		if (data["maxRMSD"].isNumeric())
			maxRMSD = data["maxRMSD"].asDouble();
		if (data["reduceConfs"].isNumeric())
			reduceConfs = data["reduceConfs"].asUInt();

		if(data["max-hits"].isNumeric())
			maxHits = data["max-hits"].asUInt();
		if(data["max-orient"].isNumeric())
			orientationsPerConf = data["max-orient"].asUInt();

		if(data["minMolWeight"].isNumeric())
		{
			minWeight = data["minMolWeight"].asDouble();
			reducedMinWeight = ThreePointData::reduceWeight(minWeight);
		}
		if(data["maxMolWeight"].isNumeric())
		{
			maxWeight = data["maxMolWeight"].asDouble();
			reducedMaxWeight = ThreePointData::reduceWeight(maxWeight);
		}

		if(data["minrotbonds"].isNumeric())
		{
			minRot = ThreePointData::reduceRotatable(data["minrotbonds"].asUInt());
		}
		if(data["maxrotbonds"].isNumeric())
			maxRot = data["maxrotbonds"].asUInt();

	}
};

//parameters for retrieving data
struct DataParameters
{
	unsigned start; //starting index
	unsigned num; //number of results to return from start; 0 means all

	SortTyp sort;
	bool reverseSort;
	bool extraInfo; //get mol info (slower)
	//extract parameters from posted data
	DataParameters() :
		start(0), num(0), sort(SortType::Undefined), reverseSort(false), extraInfo(false)
	{

	}

};


#endif /* PARAMS_H_ */
