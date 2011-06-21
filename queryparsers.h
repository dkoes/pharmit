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
 * queryparsers.h
 *
 *  Created on: Aug 17, 2010
 *      Author: dkoes
 */

#ifndef QUERYPARSERS_H_
#define QUERYPARSERS_H_

#include "pharmarec.h"
#include "tinyxml/tinyxml.h"

//base class of query file parsers
class QueryParser
{
public:
	virtual bool parse(const Pharmas& pharmas, istream& in, vector<PharmaPoint>& points) = 0;
};

//parses our own text based format
class TextQueryParser: public QueryParser
{
public:
	virtual bool parse(const Pharmas& pharmas, istream& in, vector<PharmaPoint>& points)
	{
		points.clear();
		PharmaPoint p;
		while(p.read(pharmas, in))
		{
			points.push_back(p);
		}
		return points.size() > 0;
	}

};

//parses our own json based format
class JSonQueryParser: public QueryParser
{
public:
	virtual bool parse(const Pharmas& pharmas, istream& in, vector<PharmaPoint>& points)
	{
		try {
			points.clear();
			Json::Value root; // will contains the root value after parsing.
			in >> root;
			return readPharmaPointsJSON(pharmas, root, points);
		}
		catch(std::exception& e)
		{
			return false;
		}
	}

};

//attempts to parse MOE's ph4 format
//does not do a lot of checking, since I couldn't find an official description
//of the format
class PH4Parser : public QueryParser
{
	//return pharma corresponding to a single moe name (does not handle |)
	const Pharma* getPharma(const char *moename, const Pharmas& pharmas)
	{
		unsigned len = strlen(moename);
		if(len < 3)
			return NULL;
		switch (moename[0]) {
			case 'D': //donor
				if(moename[1] != 'o' || moename[2] != 'n')
					return NULL;
				if(len > 3 && moename[3] == '2') //projected point, not yet supported
					return false;
				return pharmas.pharmaFromName("HydrogenDonor");
				break;
			case 'A': //acceptor or aromatic or anion
				if(moename[1] == 'c')
				{
					//acceptor
					if(moename[2] != 'c')
						return false;
					if(len > 3 && moename[3] == '2') //projected point, not yet supported
						return false;
					return pharmas.pharmaFromName("HydrogenAcceptor");
				}
				else if(moename[1] == 'r')
				{
					//aromatic
					if(moename[2] != 'o')
						return false;
					return pharmas.pharmaFromName("Aromatic");

				}
				else if(moename[1] == 'n')
				{
					//anion (neg)
					if(moename[2] != 'i')
						return false;
					return pharmas.pharmaFromName("NegativeIon");
				}
				else
					return NULL;
				break;
			case 'H':
				if(moename[1] != 'y' || moename[2] != 'd')
					return NULL;
				//hydrophobe
				return pharmas.pharmaFromName("Hydrophobic");
				break;
			case 'C': //cation (pos)
				if(moename[1] != 'a' || moename[2] != 't')
					return NULL;
				return pharmas.pharmaFromName("PositiveIon");
				break;
			default:
				break;
		}
		return NULL;
	}
public:
	virtual bool parse(const Pharmas& pharmas, istream& in, vector<PharmaPoint>& points)
	{
		points.clear();
		//first #moe
		string line;
		if(!getline(in, line))
			return false;
		if(line.substr(0,4) != "#moe")
			return false;
		//then #pharmacophore
		if(!getline(in, line))
			return false;
		//then scheme
		if(!getline(in, line))
			return false;
		//then look for #feature, number of points
		while(getline(in, line))
		{
			if(line.substr(0,8) == "#feature")
				break;
		}
		if(!in)
			return false;
		//features are whitespace separated, read until we find a # or eof
		while(in)
		{
			string feature;
			in >> feature;
			if(feature.size() == 0 || feature[0] == '#')
				break;

			vector<const Pharma*> types;
			//read all feature types within this feature (separated by |)
			types.push_back(getPharma(feature.c_str(), pharmas));
			size_t pos = feature.find_first_of('|');
			while(pos != string::npos)
			{
				const Pharma* p = getPharma(&feature[pos+1], pharmas);
				if(p != NULL && p != types.back())
					types.push_back(p);
				pos = feature.find_first_of('|',pos+1);
			}

			//now read point data
			string dummy;
			PharmaPoint point;

			in >> dummy;
			in >> point.x;
			in >> point.y;
			in >> point.z;
			in >> point.radius;
			in >> dummy;
			in >> dummy;
			if(!in)
				return false;

			for(unsigned i = 0, n = types.size(); i < n; i++)
			{
				if(types[i] != NULL)
				{
					point.pharma = types[i];
					points.push_back(point);
				}
			}
		}

		return true;
	}
};


//attempts to parse LigandScount's xml-based pml format
class PMLParser : public QueryParser
{
public:
	virtual bool parse(const Pharmas& pharmas, istream& in, vector<PharmaPoint>& points)
	{
		points.clear();
		TiXmlDocument doc;
		in >> doc;

		TiXmlElement *base = doc.FirstChildElement("MolecularEnvironment");
		if(base == NULL)
			base = doc.ToElement();
		if(base == NULL)
			return false;
		TiXmlElement *ph = base->FirstChildElement("pharmacophore");
		if(ph == NULL)
			return false;

		TiXmlNode *pt = NULL;
		while((pt = ph->IterateChildren(pt)) != NULL)
		{
			TiXmlElement *el = pt->ToElement();
			if(el == NULL)
				continue;
			const char *n = el->Attribute("name"); //short name
			if(n)
			{
				PharmaPoint point;
				if(el->Attribute("optional") && strcmp(el->Attribute("optional"), "true") == 0)
					point.requirements = PharmaPoint::Optional;

				switch(n[0])
				{
				case 'H': //HBD, HBA, H
					if(n[1] == 'B') //hydrogen bond
					{
						if(n[2] == 'A')
							point.pharma = pharmas.pharmaFromName("HydrogenAcceptor");
						else if(n[2] == 'D')
							point.pharma = pharmas.pharmaFromName("HydrogenDonor");
						else
							break;

						//if it points to the ligand, the target is the ligand point
						bool toligand = false;
						if(el->Attribute("pointsToLigand") && strcmp(el->Attribute("pointsToLigand"), "true") == 0)
							toligand = true;

						TiXmlElement *orig = toligand ? el->FirstChildElement("target") : el->FirstChildElement("origin");
						if(orig == NULL)
							break;
						orig->Attribute("x3", &point.x);
						orig->Attribute("y3", &point.y);
						orig->Attribute("z3", &point.z);
						orig->Attribute("tolerance",&point.radius);

						TiXmlElement *targ = toligand ? el->FirstChildElement("origin") : el->FirstChildElement("target");
						if (targ != NULL)
						{
							vector3 vec;
							targ->Attribute("x3", &vec.x());
							targ->Attribute("y3", &vec.y());
							targ->Attribute("z3", &vec.z());
							vec -= vector3(point.x, point.y, point.z);
							double t = 0;
							targ->Attribute("tolerance",&t);
							//use tolerance to compute vec pivot
							point.vecpivot = atan2(t, vec.length());
							point.vecs.push_back(vec);
						}
						points.push_back(point);
					}
					else if(n[1] == 0) //hydrophobic
					{
						TiXmlElement *pos = el->FirstChildElement("position");
						if(pos)
						{
							point.pharma = pharmas.pharmaFromName("Hydrophobic");
							pos->Attribute("x3", &point.x);
							pos->Attribute("y3", &point.y);
							pos->Attribute("z3", &point.z);
							pos->Attribute("tolerance",&point.radius);

							points.push_back(point);
						}
					}
					break;
				case 'P': //PI
				case 'N': //NI
					if(n[1] == 'I')
					{
						TiXmlElement *pos = el->FirstChildElement("position");
						if(pos)
						{
							point.pharma = pharmas.pharmaFromName(n[0] == 'P' ? "PositiveIon" : "NegativeIon");
							pos->Attribute("x3", &point.x);
							pos->Attribute("y3", &point.y);
							pos->Attribute("z3", &point.z);
							pos->Attribute("tolerance",&point.radius);

							points.push_back(point);
						}
					}
					break;
				case 'A':
					if(n[1] == 'R')
					{
						TiXmlElement *pos = el->FirstChildElement("position");
						if(pos)
						{
							point.pharma = pharmas.pharmaFromName("Aromatic");
							pos->Attribute("x3", &point.x);
							pos->Attribute("y3", &point.y);
							pos->Attribute("z3", &point.z);
							pos->Attribute("tolerance",&point.radius);

							TiXmlElement *norm = el->FirstChildElement("normal");
							if (norm != NULL)
							{
								vector3 vec;
								norm->Attribute("x3", &vec.x());
								norm->Attribute("y3", &vec.y());
								norm->Attribute("z3", &vec.z());
								vec -= vector3(point.x, point.y, point.z);
								double t = 0;
								norm->Attribute("tolerance",&t);
								//use tolerance to compute vec pivot
								point.vecpivot = atan2(t, vec.length());
								point.vecs.push_back(vec);
							}
							if(el->Attribute("optional") && strcmp(el->Attribute("optional"), "true") == 0)
								point.requirements = PharmaPoint::Optional;

							points.push_back(point);
						}
					}
					break;
				}
			}
		}

		//we store hydrogen bond points with a single combined vector, so merge any
		//hbond features at the same location with different vectors
		vector<PharmaPoint> newpoints;
		vector<bool> merged(points.size(), false);
		for(unsigned i = 0, n = points.size(); i < n; i++)
		{
			if(merged[i])
				continue;
			if(points[i].pharma->name == "HydrogenAcceptor" || points[i].pharma->name == "HydrogenDonor")
			{
				for(unsigned j = i+1; j < n; j++)
				{
					vector<vector3> vecs = points[i].vecs;
					if(points[i].x == points[j].x && points[i].y == points[j].y &&
							points[i].radius == points[j].radius
							&& points[i].z == points[j].z && points[i].pharma->name == points[j].pharma->name)
					{
						merged[j] = true;
						//grab vectors
						vecs.insert(vecs.end(), points[j].vecs.begin(), points[j].vecs.end());
					}

					if(vecs.size() > 0)
					{
						vector3 vec(0, 0, 0);
						for (unsigned v = 0, nv = vecs.size(); v < nv; v++)
						{
							vec += vecs[v];
						}
						vec /= vecs.size();
						points[i].vecs.clear();
						points[i].vecs.push_back(vec);
					}
				}
			}
			newpoints.push_back(points[i]);
		}
		swap(points, newpoints);
		return true;
	}
};


#endif /* QUERYPARSERS_H_ */
