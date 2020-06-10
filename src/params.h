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
 * params.h
 *
 *  Created on: Oct 19, 2010
 *      Author: dkoes
 *
 *      Simple parameter storing classes.
 */

#ifndef PHARMITSERVER_PARAMS_H_
#define PHARMITSERVER_PARAMS_H_
#include <iostream>
#include <boost/tuple/tuple.hpp>
#include <boost/unordered_set.hpp>
#include <vector>
#include <json/json.h>

#include "ThreePointData.h"
#include "MolProperties.h"

using namespace std;

namespace SortType {
enum SortType {Undefined, RMSD, MolWeight, NRBnds};
}
typedef SortType::SortType SortTyp;

//mol property min/max specification
struct PropFilter
{
	MolProperties::PropIDs kind;
	double min;
	double max;

	PropFilter(): kind(MolProperties::None), min(-HUGE_VAL), max(HUGE_VAL) {}
	PropFilter(MolProperties::PropIDs k, double low, double high): kind(k), min(low), max(high) {}
};

// Isolates parameters of a query.
struct QueryParameters
{
	double maxRMSD;
	unsigned reduceConfs;
	unsigned orientationsPerConf;
	unsigned maxHits; //total returned hits
	SortTyp sort; //for determining how to truncate
	bool reverseSort;
	//add support for more later

	double minWeight;
	double maxWeight;

	unsigned reducedMinWeight;
	unsigned reducedMaxWeight;
	unsigned minRot;
	unsigned maxRot;

	bool isshape;
	string subset;

	vector<PropFilter> propfilters;

	QueryParameters() :
		maxRMSD(HUGE_VAL), reduceConfs(UINT_MAX), orientationsPerConf(UINT_MAX), maxHits(UINT_MAX),
		sort(SortType::Undefined), reverseSort(false), minWeight(0), maxWeight(UINT_MAX), reducedMinWeight(0), reducedMaxWeight(UINT_MAX), minRot(0), maxRot(UINT_MAX), isshape(false)
	{

	}


	void addPropFilter(MolProperties::PropIDs p, string name, Json::Value& data)
	{
		string min("min");
		min += name;
		string max("max");
		max += name;

		if(data[min].isNumeric() || data[max].isNumeric())
		{
			PropFilter f;
			f.kind = p;
			if(data[min].isNumeric())
				f.min = data[min].asDouble();
			if(data[max].isNumeric())
				f.max = data[max].asDouble();
			propfilters.push_back(f);
		}
	}

	//extract parameters from json
	QueryParameters(Json::Value& data) :
		maxRMSD(HUGE_VAL), reduceConfs(UINT_MAX), orientationsPerConf(UINT_MAX), maxHits(UINT_MAX),sort(SortType::Undefined), reverseSort(false),
		minWeight(0), maxWeight(HUGE_VAL),reducedMinWeight(0), reducedMaxWeight(UINT_MAX), minRot(0), maxRot(UINT_MAX), isshape(false)
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

		if(data["subset"].isString())
			subset = data["subset"].asString();

		if(data["ShapeModeSelect"].isString() && data["ShapeModeSelect"].asString() == "search")
			isshape = true; //shape search

		//this sort if for truncating, do something reasonable (close to query)
		//TODO: specify in query object
		sort = SortType::RMSD;
		reverseSort = isshape;


		addPropFilter(MolProperties::LogP, "logp", data);
		addPropFilter(MolProperties::PSA, "psa", data);
		addPropFilter(MolProperties::NAromatics, "aromatics", data);
		addPropFilter(MolProperties::HBA, "hba", data);
		addPropFilter(MolProperties::HBD, "hbd", data);
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

	unsigned drawCode; //for datatables
	//extract parameters from posted data
	DataParameters() :
		start(0), num(0), sort(SortType::Undefined), reverseSort(false), extraInfo(false), drawCode(0)
	{

	}

};


#endif /* PHARMITSERVER_PARAMS_H_ */
