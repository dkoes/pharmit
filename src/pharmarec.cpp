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
 * pharmarec.cpp
 *
 *  Created on: Dec 1-, 2009
 *      Author: dkoes
 *
 *      Read in a file and identify all the pharma points.
 *      Can either use smarts or built in recognition routines.
 *      Can either print points to stdout or modify the file with representative atoms.
 *      This version uses OEChem instead of OpenBabel
 */

#include <getopt.h>
#include <iostream>
#include <cassert>
#include <string>
#include "pharmarec.h"
#include "Timer.h"
#include "ShapeConstraints.h"
#include <boost/lexical_cast.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/unordered_set.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/foreach.hpp>
#include <eigen3/Eigen/Core>

using namespace boost;
using namespace OpenBabel;

cl::opt<bool> DKoesTest("dkoes", cl::Hidden);

//default pharmacophore definitions
//The original definitions are due to Lidio Meireles and have been subsequently modified by me (dkoes)
const char *aromatic[] =
{ "a1aaaaa1", "a1aaaa1", NULL };

const char * hydrogen_donor[] =
{ "[#7!H0&!$(N-[SX4](=O)(=O)[CX4](F)(F)F)]", "[#8!H0&!$([OH][C,S,P]=O)]",
		"[#16!H0]", NULL };

const char * hydrogen_acceptor[] =
{ "[#7&!$([nX3])&!$([NX3]-*=[!#6])&!$([NX3]-[a])&!$([NX4])&!$(N=C([C,N])N)]",
		"[$([O])&!$([OX2](C)C=O)&!$(*(~a)~a)]", NULL };

const char * positive_ion[] =
{ "[+,+2,+3,+4]",
//amidine
		"[$(CC)](=N)N",
		//guanidine
		"[$(C(N)(N)=N)]", "[$(n1cc[nH]c1)]", NULL };

const char * negative_ion[] =
{ "[-,-2,-3,-4]", "C(=O)[O-,OH,OX1]", "[$([S,P](=O)[O-,OH,OX1])]",
		"c1[nH1]nnn1", "c1nn[nH1]n1", "C(=O)N[OH1,O-,OX1]", "C(=O)N[OH1,O-]",
		"CO(=N[OH1,O-])",
		//trifluoromethyl sulfonamide
		"[$(N-[SX4](=O)(=O)[CX4](F)(F)F)]", NULL };

const char *hydrophobic[] =
		{
				"a1aaaaa1",
				"a1aaaa1",
				//branched terminals as one point
				"[$([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(**[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]",
				"[$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
				"*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
				//simple rings only; need to combine points to get good results for 3d structures
				"[C&r3]1~[C&r3]~[C&r3]1",
				"[C&r4]1~[C&r4]~[C&r4]~[C&r4]1",
				"[C&r5]1~[C&r5]~[C&r5]~[C&r5]~[C&r5]1",
				"[C&r6]1~[C&r6]~[C&r6]~[C&r6]~[C&r6]~[C&r6]1",
				"[C&r7]1~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]1",
				"[C&r8]1~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]1",
				//aliphatic chains
				"[CH2X4,CH1X3,CH0X2]~[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
				"[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]",
				"[$([CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]",
				// sulfur (apparently)
				"[$([S]~[#6])&!$(S~[!#6])]", NULL };

const vector<Pharma> defaultPharmaVec = assign::list_of(
		Pharma(0, "Aromatic", aromatic, 18, 1.0, 0.1))(
		Pharma(1, "HydrogenDonor", hydrogen_donor, 1, 1.0, 0.1))(
		Pharma(2, "HydrogenAcceptor", hydrogen_acceptor, 89, 1.0, 0.1))(
		Pharma(3, "PositiveIon", positive_ion, 7, 1.0, 0.1))(
		Pharma(4, "NegativeIon", negative_ion, 8, 1.0, 0.1))(
		Pharma(5, "Hydrophobic", hydrophobic, 6, 1.0, 2.0));

//TODO: make this user-configurable
static boost::unordered_map<string, PharmaInteract> pharmaInteractions =
		assign::map_list_of("Aromatic", PharmaInteract(0, 5, 1))(
				"HydrogenDonor", PharmaInteract(2, 4, 1))("HydrogenAcceptor",
				PharmaInteract(1, 4, 1))("PositiveIon", PharmaInteract(4, 5, 1))(
				"NegativeIon", PharmaInteract(3, 5, 1))("Hydrophobic",
				PharmaInteract(5, 6, 3));
//reduced set of pharmacophore definitions for proteins
const char * positive_ion_protein[] =
{ "[+,+2,+3,+4]",
//amidine
//guanidine
		"[$(C(N)(N)=N)]", "[$(n1cc[nH]c1)]", NULL };

const char * negative_ion_protein[] =
{ "[-,-2,-3,-4]", "C(=O)[O-,OH,OX1]", NULL };

const char *hydrophobic_protein[] =
		{
				"a1aaaaa1",
				"a1aaaa1",
				//branched terminals as one point
				"[$([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(**[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]",
				"[$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",

				"[CH2X4,CH1X3,CH0X2]~[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
				"[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]",
				"[$([S]~[#6])&!$(S~[!#6])]", NULL };

static const vector<Pharma> proteinPharmaVec = assign::list_of(
		Pharma(0, "Aromatic", aromatic, 18, 1.0, 0.1))(
		Pharma(1, "HydrogenDonor", hydrogen_donor, 1, 1.0, 0.1))(
		Pharma(2, "HydrogenAcceptor", hydrogen_acceptor, 89, 1.0, 0.1))(
		Pharma(3, "PositiveIon", positive_ion_protein, 7, 1.0, 0.1))(
		Pharma(4, "NegativeIon", negative_ion_protein, 8, 1.0, 0.1))(
		Pharma(5, "Hydrophobic", hydrophobic_protein, 6, 1.0, 2.0));
static Pharmas proteinPharmas(proteinPharmaVec);

//setup lookup tables, pharmas must be initialized
//use a manually managed array to make sure pointers stay valid
void Pharmas::initialize(const vector<Pharma>& ps)
{
	if (pharmas != NULL)
		delete[] pharmas;
	nameLookup.clear();
	pharmas = new Pharma[ps.size()];
	numPharmas = ps.size();
	for (unsigned i = 0; i < numPharmas; i++)
	{
		pharmas[i] = ps[i];
		nameLookup[pharmas[i].name] = i;
	}
}

bool Pharma::operator==(const Pharma& rhs) const
{
	if (name != rhs.name)
		return false;
	if (atomic_number_label != rhs.atomic_number_label)
		return false;
	if (index != rhs.index)
		return false;
	if (defaultSearchRadius != rhs.defaultSearchRadius)
		return false;
	if (getVectors != rhs.getVectors)
		return false;

	if (smarts.size() != rhs.smarts.size())
		return false;

	for (unsigned i = 0, n = smarts.size(); i < n; i++)
	{
		if (smarts[i].GetSMARTS() != rhs.smarts[i].GetSMARTS())
			return false;
	}

	return true;
}

//return pharma from name
const Pharma* Pharmas::pharmaFromName(const string& name) const
{
	boost::unordered_map<string, unsigned>::const_iterator pos = nameLookup.find(name);
	if (pos == nameLookup.end())
	{
		return NULL;
	}
	return &pharmas[pos->second];
}

bool Pharmas::operator==(const Pharmas& rhs) const
{
	if (numPharmas != rhs.numPharmas)
		return false;
	for (unsigned i = 0; i < numPharmas; i++)
	{
		if (pharmas[i] != rhs.pharmas[i])
			return false;
	}
	return true;
}

//override default search radius to a single specified constant
//because I think uniform 1.0 is the best default but made the mistake
//of embedding this information in the database
void Pharmas::setDefaultSearchRadius(double val = 1.0)
{
	for (unsigned i = 0; i < numPharmas; i++)
	{
		pharmas[i].defaultSearchRadius = val;
	}
}

static void readVectors(stringstream& str, vector<vector3>& vecs)
{
	while (str && ::isblank(str.peek()))
		str.get();

	while (str && str.peek() == '[')
	{
		double x, y, z;
		str >> x;
		str >> y;
		str >> z;
		vecs.push_back(vector3(x, y, z));
		char c;
		str >> c; //end bracket

		while (str && ::isblank(str.peek()))
			str.get();
	}

}

//read in query pharma
bool PharmaPoint::read(const Pharmas& pharmas, istream &in)
{
	string line;
	getline(in, line);
	stringstream str(line);

	//init
	size = 0;
	vecpivot = 0;
	requirements = Required;
	radius = 0;

	string name;
	str >> name;
	pharma = pharmas.pharmaFromName(name);
	if (pharma == NULL)
		return false;

	str >> x;
	str >> y;
	str >> z;

	if (!str)
		return false;
	readVectors(str, vecs);

	str >> radius;
	char r = 0;
	str >> r;
	switch (r)
	{
	case 'r':
	case 'R':
		requirements = Required;
		break;
	case 'n':
	case 'N':
		requirements = NotPresent;
		break;
	case 'o':
	case 'O':
		requirements = Optional;
		break;
	}
	str >> vecpivot;

	return true;
}

//write out textual description of pharma point
ostream &operator<<(ostream &stream, const PharmaPoint& obj)
{
	stream << obj.pharma->name << " ";
	stream << obj.x << ' ' << obj.y << ' ' << obj.z;

	for (unsigned i = 0, n = obj.vecs.size(); i < n; i++)
	{
		stream << " [ " << obj.vecs[i].x() << " " << obj.vecs[i].y() << " "
				<< obj.vecs[i].z() << " ]";
	}
	return stream;
}

//read from a json formatted stream, return true if successfull
bool readPharmaPointsJSON(const Pharmas& pharmas, Json::Value& root,
		vector<PharmaPoint>& points)
{
	try
	{
		points.clear();
		Json::Value jpoints = root["points"];
		for (unsigned i = 0, n = jpoints.size(); i < n; i++)
		{
			Json::Value jpnt = jpoints[i];
			if (jpnt.isMember("enabled") && !jpnt["enabled"].asBool())
				continue;
			string name = jpnt["name"].asString();
			double radius = .1;
			unsigned minsize = 0;
			unsigned maxsize = UINT_MAX;
			PharmaPoint::PointRequirements required = PharmaPoint::Required;
			if (jpnt.isMember("radius"))
				radius = jpnt["radius"].asDouble();
			if (jpnt.isMember("requirement") && jpnt["required"].isString())
			{
				string req = jpnt["requirement"].asString();
				algorithm::to_lower(req);
				if (req == "optional")
					required = PharmaPoint::Optional;
				else if (req == "notpresent")
					required = PharmaPoint::NotPresent;
			}

			if (jpnt.isMember("minsize") && !jpnt["minsize"].isNull())
			{
				if (jpnt["minsize"].isString())
				{
					if (jpnt["minsize"].asString().length() > 0)
					{
						try
						{
							minsize = lexical_cast<unsigned>(
									jpnt["minsize"].asString());
						} catch (std::exception& e)
						{
						}
					}
				}
				else
					minsize = jpnt["minsize"].asUInt();
			}

			if (jpnt.isMember("maxsize") && !jpnt["maxsize"].isNull())
			{
				if (jpnt["maxsize"].isString())
				{
					if (jpnt["maxsize"].asString().length() > 0)
					{
						try
						{
							maxsize = lexical_cast<unsigned>(
									jpnt["maxsize"].asString());
						} catch (std::exception& e)
						{
						}
					}
				}
				else
					maxsize = jpnt["maxsize"].asUInt();
			}

			double x = jpnt["x"].asDouble();
			double y = jpnt["y"].asDouble();
			double z = jpnt["z"].asDouble();
			const Pharma* p = pharmas.pharmaFromName(name);
			if (p == NULL)
				continue;
			unsigned size = 1;
			if (jpnt.isMember("size"))
				size = jpnt["size"].asUInt();

			PharmaPoint pnt;
			pnt.pharma = p;
			pnt.x = x;
			pnt.y = y;
			pnt.z = z;
			pnt.size = size;
			pnt.maxSize = maxsize;
			pnt.minSize = minsize;

			if (!jpnt.isMember("vector_on") || jpnt["vector_on"].asBool())
			{
				if (jpnt.isMember("svector") && !jpnt["svector"].isNull())
				{
					pnt.vecs.push_back(
							vector3(jpnt["svector"]["x"].asDouble(),
									jpnt["svector"]["y"].asDouble(),
									jpnt["svector"]["z"].asDouble()));
				}
				else if (jpnt.isMember("vector"))
				{ //queries typically have only one vector
					for (unsigned v = 0, nv = jpnt["vector"].size(); v < nv;
							v++)
					{
						pnt.vecs.push_back(
								vector3(jpnt["vector"][v]["x"].asDouble(),
										jpnt["vector"][v]["y"].asDouble(),
										jpnt["vector"][v]["z"].asDouble()));
					}
				}
			}

			pnt.radius = radius;
			pnt.requirements = required;
			points.push_back(pnt);
		}

	} catch (std::exception& e)
	{
		//poorly formated json
		cerr << "Parse " << e.what() << "\n";
		return false;
	}
	return true;
}

//return location of max degree, return false if no nonzero degrees
static bool findMax(vector<unsigned>& degrees, unsigned n, unsigned& maxi)
{
	unsigned m = 0;
	for (unsigned i = 0; i < n; i++)
	{
		if (degrees[i] > m)
		{
			m = degrees[i];
			maxi = i;
		}
	}
	return m > 0;
}

struct Coord
{
	double x;
	double y;
	double z;

	Coord(double X, double Y, double Z) :
			x(X), y(Y), z(Z)
	{
	}
	Coord() :
			x(0), y(0), z(0)
	{
	}

	bool operator==(const Coord& rhs) const
	{
		return x == rhs.x && y == rhs.y && z == rhs.z;
	}
};

static size_t hash_value(const Coord & t)
{
	size_t seed = 0;
	hash_combine(seed, t.x);
	hash_combine(seed, t.y);
	hash_combine(seed, t.z);

	return seed;
}

//sort an indexing array into points by distance to index point
class DistanceSorter
{
	unsigned index;
	const vector<PharmaPoint>& points;
public:
	DistanceSorter(unsigned i, const vector<PharmaPoint>& pts) :
			index(i), points(pts)
	{
	}

	bool operator()(unsigned i, unsigned j) const
	{
		double d1 = PharmaPoint::pharmaDist(points[index], points[i]);
		double d2 = PharmaPoint::pharmaDist(points[index], points[j]);
		return d1 < d2;
	}

};

// greedily create clusters of points that are all separated by dist
//and replace with the cluster centers
//the common cases if for all points we want to be clustered to be in a clique
//so hopefully no need for more complex algorithms
//mutates points
static void clusterPoints(const Pharma *pharma, vector<PharmaPoint>& points,
		const vector<vector<Coord> >& ptatoms, double dist)
{
	unsigned n = points.size();
	vector<PharmaPoint> npts;
	npts.reserve(n);
	vector<vector<double> > distance(n, vector<double>(n, 0));
	vector<unsigned> degrees(n, 0);

	//compute distances and degrees (number of close neighbors)
	for (unsigned i = 0; i < n; i++)
	{
		for (unsigned j = 0; j < i; j++)
		{
			double d = PharmaPoint::pharmaDist(points[i], points[j]);
			distance[j][i] = distance[i][j] = d;

			if (d <= dist)
			{
				degrees[i]++;
				degrees[j]++;
			}
		}
	}

	unsigned order[n];
	//add all zero degree points
	for (unsigned i = 0; i < n; i++)
	{
		order[i] = i;
		if (degrees[i] == 0)
			npts.push_back(points[i]);
	}

	unsigned maxi = 0;
	while (findMax(degrees, n, maxi))
	{
		vector<unsigned> cluster;
		cluster.push_back(maxi);
		degrees[maxi] = 0;

		DistanceSorter sorter(maxi, points);
		sort(order, order + n, sorter);
		//add points that are within threshold of all member of the cluster
		for (unsigned o = 0; o < n; o++)
		{
			unsigned i = order[o];
			if (degrees[i] > 0)
			{
				unsigned j = 0;
				for (unsigned csz = cluster.size(); j < csz; j++)
				{
					if (distance[i][cluster[j]] > dist)
						break;
				}
				if (j == cluster.size()) //close to all current cluster members
				{
					degrees[i] = 0;
					cluster.push_back(i);
				}
			}
		}

		//create point with average of atom coords
		unordered_set<Coord> atoms;
		for (unsigned i = 0, csz = cluster.size(); i < csz; i++)
		{
			unsigned index = cluster[i];
			for (unsigned j = 0, na = ptatoms[index].size(); j < na; j++)
			{
				atoms.insert(ptatoms[index][j]);
			}
		}
		PharmaPoint pt(pharma);
		BOOST_FOREACH(const Coord& a, atoms)
				{
					pt.x += a.x;
					pt.y += a.y;
					pt.z += a.z;
				}
		pt.x /= atoms.size();
		pt.y /= atoms.size();
		pt.z /= atoms.size();
		pt.size = atoms.size();

		npts.push_back(pt);
	}

	points.swap(npts);
}

/* Identify the pharma points of mol and put them in points.
 * Note that this is very inefficient in processing sequences of
 * multi-conformer molecules.
 * */
void getPharmaPoints(const Pharmas& pharmas, OBMol& mol,
		vector<PharmaPoint>& points)
{
	vector<vector<PharmaPoint> > pts;
	points.clear();
	getPharmaPointsMC(pharmas, mol, pts);
	if (pts.size() > 0)
		points.swap(pts[0]);
}

//return true if mol is a PharmaGist mol2 and fill in points with first model
bool isPharmaGist(const Pharmas& pharmas, const string& mol,
		vector<PharmaPoint>& points)
{
	if (mol.substr(0, 17) != "@<TRIPOS>MOLECULE")
		return false;

	istringstream str(mol);

	bool readatoms = false;
	string line;
	while (getline(str, line))
	{
		if (line == "@<TRIPOS>ATOM")
			readatoms = true;
		else if (line == "@<TRIPOS>BOND")
			return points.size() > 0;
		else if (readatoms)
		{
			stringstream l(line);
			unsigned num;
			string type;
			double x, y, z;
			l >> num >> type >> x >> y >> z;

			if(type == "DON")
			{
				points.push_back(PharmaPoint(x,y,z,pharmas.pharmaFromName("HydrogenDonor")));
			}
			else if(type == "ACC")
			{
				points.push_back(PharmaPoint(x,y,z,pharmas.pharmaFromName("HydrogenAcceptor")));
			}
			else if(type == "ANI")
			{
				points.push_back(PharmaPoint(x,y,z,pharmas.pharmaFromName("NegativeIon")));
			}
			else if(type == "CAT")
			{
				points.push_back(PharmaPoint(x,y,z,pharmas.pharmaFromName("PositiveIon")));
			}
			else if(type == "HYD")
			{
				points.push_back(PharmaPoint(x,y,z,pharmas.pharmaFromName("Hydrophobic")));
			}
			else if(type == "AR")
			{
				points.push_back(PharmaPoint(x,y,z,pharmas.pharmaFromName("Aromatic")));
			}
			else
			{
				points.clear();
				return false;
			}

		}
	}
	points.clear();
	return false;
}

void getPharmaPointsMC(const Pharmas& pharmas, OBMol& mol,
		vector<vector<PharmaPoint> >& points)
{
	points.clear();
	points.resize(mol.NumConformers());
	//for each kind of pharma
	for (int p = 0, np = pharmas.size(); p < np; p++)
	{
		const Pharma* pharma = pharmas[p];

		vector<vector<PharmaPoint> > cpoints(mol.NumConformers());
		vector<vector<vector<Coord> > > ccoords(mol.NumConformers());
		//foreach smart for the pharma
		for (int s = 0, ns = pharma->smarts.size(); s < ns; s++)
		{
			const OBSmartsPattern& matcher = pharma->smarts[s];
			vector<vector<int> > maplist;
			if (matcher.Match(mol, maplist, OBSmartsPattern::AllUnique))
			{
				for (unsigned c = 0, nc = mol.NumConformers(); c < nc; c++)
				{
					//foreach map, record the point as the average of the matched atoms
					for (int i = 0, n = maplist.size(); i < n; i++)
					{
						mol.SetConformer(c);
						ccoords[c].resize(ccoords[c].size() + 1);

						PharmaPoint point(pharma);
						point.x = point.y = point.z = 0;

						//get each atom
						int numatoms = maplist[i].size();
						point.size = numatoms;
						for (int a = 0; a < numatoms; a++)
						{
							OBAtom *atom = mol.GetAtom(maplist[i][a]);
							point.x += atom->x();
							point.y += atom->y();
							point.z += atom->z();
							if (pharma->clusterLimit > 0)
								ccoords[c].back().push_back(
										Coord(atom->x(), atom->y(), atom->z()));
						}
						//take average
						if (numatoms > 0)
						{
							point.x /= numatoms;
							point.y /= numatoms;
							point.z /= numatoms;
						}

						if (pharma->getVectors != NULL)
							pharma->getVectors(maplist[i], mol, point);
						cpoints[c].push_back(point);
					}
				}
			}
		}
		//now add all points for each conformation, clustering as necessary
		for (unsigned c = 0, nc = cpoints.size(); c < nc; c++)
		{
			if (pharma->clusterLimit > 0)
			{
				clusterPoints(pharma, cpoints[c], ccoords[c],
						pharma->clusterLimit);
			}

			//add (possibly clustered) conformer points
			points[c].insert(points[c].end(), cpoints[c].begin(),
					cpoints[c].end());
		}
	}
}

/* add the passed pharma points to the molecule model */
void addPharmaPoints(OBMol& mol, vector<PharmaPoint>& points)
{
	for (int p = 0, np = points.size(); p < np; p++)
	{
		OBAtom atom;
		atom.SetVector(points[p].x, points[p].y, points[p].z);
		atom.SetAtomicNum(points[p].pharma->atomic_number_label);
		mol.AddAtom(atom);
	}
}

//convert from vector version to json version
bool convertPharmaJson(Json::Value& root, const vector<PharmaPoint>& points)
{
	//add all points
	for (unsigned i = 0, n = points.size(); i < n; i++)
	{
		Json::Value& pt = root["points"][i];
		pt["name"] = points[i].pharma->name;
		pt["x"] = points[i].x;
		pt["y"] = points[i].y;
		pt["z"] = points[i].z;
		pt["size"] = points[i].size;

		vector3 sum;
		unsigned cnt = 0;
		for (unsigned j = 0, nv = points[i].vecs.size(); j < nv; j++)
		{
			if(isfinite(points[i].vecs[j].length())) //really do not want Nan - results in bad json
			{
				vector3 v = points[i].vecs[j];
				pt["vector"][j]["x"] = v.x();
				pt["vector"][j]["y"] = v.y();
				pt["vector"][j]["z"] = v.z();
				sum += v.normalize();
				cnt++;
			}
		}

		if(sum.length_2() > 0.005) //didn't just cancel out (e.g. aromatics)
		{
			//take average
			sum /= cnt;
			pt["svector"]["x"] = sum.x();
			pt["svector"]["y"] = sum.y();
			pt["svector"]["z"] = sum.z();
		}
		else if(cnt > 0) //take first
		{
			pt["svector"]["x"] = points[i].vecs[0].x();
			pt["svector"]["y"] = points[i].vecs[0].y();
			pt["svector"]["z"] = points[i].vecs[0].z();
		}

		pt["radius"] = points[i].radius;
		pt["enabled"] = points[i].requirements != PharmaPoint::NotPresent;
	}

	return true;
}

//look for hydrogens bonded to first atom, take bond to H as vector
static void genHDonorPointVector(const vector<int>& atom_indexes,
		const OBMol& mol, PharmaPoint& pnt)
{
	if (atom_indexes.size() == 0)
	{
		return;
	}
	OBAtom *atom = mol.GetAtom(atom_indexes[0]);
	OBBondIterator i;
	for (OBAtom* nbor = atom->BeginNbrAtom(i); nbor;
			nbor = atom->NextNbrAtom(i))
	{
		if (nbor->GetAtomicNum() == 1)
		{
			vector3 vec1 = atom->GetVector();
			vector3 vec2 = nbor->GetVector();
			vector3 vec = vec2 - vec1;

			if (!isfinite(vec.x()))
				continue; //sometimes, despite our best efforts, H has same coords as N

			vec.normalize();
			pnt.vecs.push_back(vec);
		}
	}
}

//acceptors are trickier.. if a hydrogen exists, use that direction,
//otherwise use opposite of average of existing bonds
static void genHAcceptorPointVector(const vector<int>& atom_indexes,
		const OBMol& mol, PharmaPoint& pnt)
{
	if (atom_indexes.size() == 0)
	{
		return;
	}
	OBAtom *atom = mol.GetAtom(atom_indexes[0]);
	OBBondIterator i;
	bool foundvec = false;
	for (OBAtom* nbor = atom->BeginNbrAtom(i); nbor;
			nbor = atom->NextNbrAtom(i))
	{
		if (nbor->GetAtomicNum() == 1)
		{
			vector3 vec1;
			vector3 vec2;
			vec1 = atom->GetVector();
			vec2 = nbor->GetVector();

			vector3 vec = vec2 - vec1;
			if (!isfinite(vec.x()))
				continue; //sometimes, despite our best efforts, H has same coords as N

			vec.normalize();

			pnt.vecs.push_back(vec);
			foundvec = true;
		}
	}

	if (!foundvec)
	{
		vector3 avebond(0, 0, 0);
		vector3 avec = atom->GetVector();
		unsigned cnt = 0;
		OBBondIterator i;
		for (OBAtom* nbor = atom->BeginNbrAtom(i); nbor;
				nbor = atom->NextNbrAtom(i))
		{
			avebond += nbor->GetVector() - avec;
			cnt++;
		}
		avebond /= cnt;

		avebond *= -1;
		avebond.normalize();

		if(isfinite(avebond.x()))
			pnt.vecs.push_back(avebond);
	}
}

//just pick some normal to the ring
static void genAromaticPointVector(const vector<int>& atom_indexes,
		const OBMol& mol, PharmaPoint& pnt)
{
	if (atom_indexes.size() < 3)
	{
		return;
	}

	vector3 base = mol.GetAtom(atom_indexes[0])->GetVector();
	vector3 pt1 = mol.GetAtom(atom_indexes[1])->GetVector();
	vector3 pt2 = mol.GetAtom(atom_indexes[2])->GetVector();

	vector3 norm = cross((pt2 - base), (pt1 - base));
	norm.normalize();

	pnt.vecs.push_back(norm);
	pnt.vecs.push_back(-norm);
}

void Pharma::setVectorFn(genPointVectorFn fn)
{
	if (fn != NULL)
	{
		getVectors = fn;
	}
	else
	{
		if (boost::starts_with(name, "Aromatic"))
			getVectors = genAromaticPointVector;
		else if (boost::starts_with(name,"HydrogenDonor"))
			getVectors = genHDonorPointVector;
		else if (boost::starts_with(name, "HydrogenAcceptor"))
			getVectors = genHAcceptorPointVector;
		else
			getVectors = NULL;
	}
}

//write out pharma definitions
void Pharmas::write(ostream& out) const
{
	for (unsigned i = 0; i < numPharmas; i++)
	{
		out << pharmas[i].name << " ";
		out << pharmas[i].atomic_number_label << " ";
		assert(i == pharmas[i].index);
		out << i << " ";
		out << pharmas[i].defaultSearchRadius << " ";
		//0, no vectors
		//1, default
		//anything else has to be done programatically
		if (pharmas[i].getVectors != NULL)
			out << "1 ";
		else
			out << "0 ";

		out << pharmas[i].clusterLimit << " ";
		out << "\n";
		//now list all smarts
		for (unsigned s = 0, ns = pharmas[i].smarts.size(); s < ns; s++)
		{
			out << pharmas[i].smarts[s].GetSMARTS() << "\n";
		}

		//blank line between pharmas
		out << "\n";
	}
}

//read in pharma definitions
bool Pharmas::read(istream& in)
{
	string line;
	vector<Pharma> pharmavec;
	while (getline(in, line))
	{
		//first line has all info
		Pharma p;
		stringstream str(line);
		str >> p.name;
		str >> p.atomic_number_label;
		str >> p.index;
		if (p.index != pharmavec.size())
			p.index = pharmavec.size(); //making this explicit was really stupid, just use order from file
		str >> p.defaultSearchRadius;
		int vecid = 0;
		str >> vecid;
		if (vecid == 1)
			p.setVectorFn();
		else if (vecid > 1)
			return false; //unknown vector generator

		str >> p.clusterLimit;

		if (!str) //fail bit got set
			return false;

		//now read in smarts
		while (getline(in, line))
		{
			if (line.length() == 0)
				break; //empty line

			p.smarts.push_back(OBSmartsPattern());
			p.smarts.back().Init(line);
		}
		if (p.smarts.size() == 0)
			return false;
		pharmavec.push_back(p);
	}

	initialize(pharmavec);
	return true;
}

//output a set of points in json format that is suitable for query based on
//the data in moldata
bool jsonPharmaQuery(const Pharmas& pharmas, Json::Value& root,
		const string& moldata, OBFormat *format, const string& recdata,
		OBFormat *rformat)
{
	using namespace Eigen;
	vector<PharmaPoint> points;
	vector<PharmaPoint> disabled;
	vector<Vector3d> ipoints;

	if (isPharmaGist(pharmas, moldata, points))
		; //got points from the file
	else
	{
		OBMol mol;
		OBConversion conv;
		conv.SetInFormat(format);

		if (!conv.ReadString(&mol, moldata))
			return false;

		if (rformat == NULL || recdata.length() == 0)
		{
			getPharmaPoints(pharmas, mol, points);
		}
		else //compute interaction pharma
		{
			conv.SetInFormat(rformat);
			OBMol rec;
			if (!conv.ReadString(&rec, recdata))
			{
				getPharmaPoints(pharmas, mol, points);
			}
			else
			{
				getInteractionPoints(pharmas, rec, mol, points, disabled);
				points.insert(points.end(), disabled.begin(), disabled.end());
				ShapeConstraints::computeInteractionPoints(mol, rec, ipoints);
			}
		}
	}

	if (!convertPharmaJson(root, points))
		return false;

	//set disabled points
	unsigned cutoff = points.size() - disabled.size();
	for(unsigned i = 0, n = points.size(); i < n; i++)
	{
		root["points"][i]["enabled"] = (i < cutoff);

	}

	//add inclusion spheres
	for(unsigned i = 0, n = ipoints.size(); i < n; i++)
	{
		const Vector3d& ipt = ipoints[i];
		unsigned indx = points.size() + i;
		Json::Value& pt = root["points"][indx];
		pt["name"] = "InclusionSphere";
		pt["x"] = ipt.x();
		pt["y"] = ipt.y();
		pt["z"] = ipt.z();
		pt["radius"] = 1.0; //default radius should be what?
		pt["enabled"] = false;
	}

	return true;
}

//OpenBabel SMARTs matching is a little slow, so for when we're matching large
//proteins use an optimizes set of pharmas, and then convert
void getProteinPharmaPoints(const Pharmas& pharmas, OBMol& mol,
		vector<PharmaPoint>& points)
{
	getPharmaPoints(proteinPharmas, mol, points);
	//now convert to requested pharmas based on name equality
	for (unsigned i = 0, n = points.size(); i < n; i++)
	{
		points[i].pharma = pharmas.pharmaFromName(points[i].pharma->name);
	}
}

//prune atoms from rec that are more than dist from ligand
//this is approximate since I'm just using the bounding box
static OBMol getPrunedReceptor(OBMol& rec, OBMol& ligand, double dist)
{
	if(ligand.NumAtoms() == 0) return rec;

	double min[3] = {HUGE_VAL, HUGE_VAL, HUGE_VAL};
	double max[3] = {-HUGE_VAL, -HUGE_VAL, -HUGE_VAL};

	//get bounding box of ligand
	FOR_ATOMS_OF_MOL(a, ligand)
	{
		vector3 pt = a->GetVector();
		for(unsigned i = 0; i < 3; i++) {
			if(pt[i] > max[i]) max[i] = pt[i];
			if(pt[i] < min[i]) min[i] = pt[i];
		}
	}

	//now expand by dist
	for(unsigned i = 0; i < 3; i++)
	{
		max[i] += dist;
		min[i] -= dist;
	}

	//identify receptor atoms outside of box
	vector<OBAtom*> toremove;
	OBMol ret = rec;
	FOR_ATOMS_OF_MOL(a, ret)
	{
		vector3 pt = a->GetVector();
		if(pt[0] > max[0] || pt[1] > max[1] || pt[2] > max[2] ||
				pt[0] < min[0] || pt[1] < min[1] || pt[2] < min[2])
		{
			toremove.push_back(&(*a));
		}
	}

	ret.BeginModify();
	for(unsigned i = 0, n = toremove.size(); i < n; i++)
	{
		ret.DeleteAtom(toremove[i],true);
	}
	ret.EndModify();
	return ret;
}

//compute pharmacophore of the interaction
//calculate both ligand and receptor points and then
//enable only those ligand points that are close to complimentary receptor points
void getInteractionPoints(const Pharmas& pharmas, OBMol& receptor,
		OBMol& ligand, vector<PharmaPoint>& points,
		vector<PharmaPoint>& screenedout)
{
	points.clear();
	screenedout.clear();

	vector<PharmaPoint> ligandpoints;
	vector<PharmaPoint> receptorpoints;

	getPharmaPoints(pharmas, ligand, ligandpoints);
	OBMol prunedrec = getPrunedReceptor(receptor, ligand, 10); //prune receptor
	getProteinPharmaPoints(pharmas, prunedrec, receptorpoints); //could potentially prune a very large molecule..

	//remove anything without interacting info
	vector<PharmaPoint> interactpoints;
	for (unsigned i = 0, n = ligandpoints.size(); i < n; i++)
	{
		const PharmaInteract& interact =
				pharmaInteractions[ligandpoints[i].pharma->name];
		if (interact.maxDist > 0)
		{
			interactpoints.push_back(ligandpoints[i]);
		}
		else
		{
			screenedout.push_back(ligandpoints[i]);
		}
	}
	//collate receptor points
	vector<vector<PharmaPoint> > rinteractpoints(pharmas.size());

	for (unsigned i = 0, n = receptorpoints.size(); i < n; i++)
	{
		const PharmaInteract& interact =
				pharmaInteractions[receptorpoints[i].pharma->name];
		if (interact.maxDist > 0)
		{
			rinteractpoints[receptorpoints[i].pharma->index].push_back(
					receptorpoints[i]);
		}
	}

	//screen ligand points
	for (unsigned i = 0, n = interactpoints.size(); i < n; i++)
	{
		const PharmaPoint& l = interactpoints[i];
		unsigned cnt = 0;
		const PharmaInteract& I = pharmaInteractions[l.pharma->name];
		const vector<PharmaPoint>& rvec = rinteractpoints[I.complement];
		unsigned j, m;
		for (j = 0, m = rvec.size(); j < m; j++)
		{
			const PharmaPoint& rp = rvec[j];
			double d = PharmaPoint::pharmaDist(l, rp);
			if (d <= I.maxDist)
				cnt++;
			if (cnt >= I.minMatch)
			{
				//just hydrogen bond features and set vector
				if (l.pharma->name == "HydrogenAcceptor"
						|| l.pharma->name == "HydrogenDonor")
				{
					PharmaPoint lp = l;
					vector3 rv(rp.x, rp.y, rp.z);
					vector3 lv(l.x, l.y, l.z);
					lp.vecs.clear();
					lp.vecs.push_back(rv - lv);
					//assume identical points with different vectors are right next to each other
					//and only do one
					if (points.size() == 0 || points.back().pharma != l.pharma
							|| points.back().x != l.x || points.back().y != l.y
							|| points.back().z != l.z)
						points.push_back(lp);
					else if(points.size() > 0 && points.back().vecs.size() > 0)//identical points, add vector
					{
						//alternatively, we could store all vectors..
						double len = lp.vecs[0].length_2();
						double oldlen = points.back().vecs[0].length_2();
						if(len < oldlen)
						{
							points.back().vecs[0] = lp.vecs[0];
						}
					}
				}
				else
				{
					points.push_back(l);
				}
				break;
			}
		}
		if (j == m)
			screenedout.push_back(l);
	}
}
