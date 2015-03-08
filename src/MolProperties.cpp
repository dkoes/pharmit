/*
 * MolProperties.cpp
 *
 *  Created on: Mar 5, 2015
 *      Author: dkoes
 */

#include "MolProperties.h"
#include <boost/filesystem.hpp>
#include <boost/assign.hpp>
#include <openbabel/descriptor.h>

using namespace OpenBabel;
using namespace boost;

void MolProperties::calculate(OpenBabel::OBMol& mol, unsigned long id)
{
	uniqueid = id;
	//various openbabel properties (what is calculated for obprop)
	vector<OBRing*> vr;
	vr = mol.GetSSSR();
	num_rings = vr.size();
	num_aromatics = 0;
	for (unsigned i = 0, n = vr.size(); i < n; i++)
	{
		if (vr[i]->IsAromatic())
			num_aromatics++;
	}
	OBDescriptor* desc;
	desc = OBDescriptor::FindType("logP");
	if (desc)
		logP = desc->Predict(&mol);
	desc = OBDescriptor::FindType("TPSA");
	if (desc)
		psa = desc->Predict(&mol);

	//HBA and HBD are set by pharmacophore detector
}

vector<const char*> MolProperties::fileNames = assign::list_of("prop_uniqueid")
		("prop_nrings")("prop_naromatics")("prop_logP")("prop_psa")("prop_hba")("prop_hbd");

//create files in spec
void MolProperties::createFiles(const boost::filesystem::path& dbpath,
		MolProperties::PropFiles& files)
{
	files.clear();

	//name files prop_<prop>
	for (unsigned i = 0, n = MolProperties::fileNames.size(); i < n; i++)
	{
		filesystem::path fpath = dbpath / fileNames[i];
		assert(!filesystem::exists(fpath));
		FILE *f = fopen(fpath.string().c_str(), "w+"); //create
		assert(f);
		files.push_back(f);
	}

}

void MolProperties::write(unsigned mid, const PropFiles& files)
{
	//there are gaps in the mids if molecules get too large
	//which certainly happens with our crazy names..

#define WRITE(NAME, I) fseek(files[I], sizeof(NAME)*mid, SEEK_SET); fwrite(&NAME, sizeof(NAME), 1, files[I]);
	//write to individual files that were opened with createFiles

	WRITE(uniqueid, 0);
	WRITE(num_rings, 1);
	WRITE(num_aromatics, 2);
	WRITE(logP, 3);
	WRITE(psa, 4);
	WRITE(hba, 5);
	WRITE(hbd, 6);
}


void MolProperties::setHB(const vector<PharmaPoint>& pts)
{
	hbd = hba = 0;

	for(unsigned i = 0, n = pts.size(); i < n; i++)
	{
		if(pts[i].pharma->name == "HydrogenDonor")
			hbd++;
		else if(pts[i].pharma->name == "HydrogenAcceptor")
			hba++;
	}
}
