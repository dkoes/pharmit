/*
 * MolProperties.h
 *
 *  Created on: Mar 5, 2015
 *      Author: dkoes
 */

#ifndef MOLPROPERTIES_H_
#define MOLPROPERTIES_H_

#include <vector>
#include <openbabel/mol.h>
#include "MMappedRegion.h"
#include "pharmarec.h"

using namespace std;

//various properties that we compute per molecule (not unique to conformers)
struct MolProperties
{
	unsigned long uniqueid; //user provided namespace for molecules
	//openbabel properties
	unsigned char num_rings;
	unsigned char num_aromatics;
	//NOTE: hydrogen bound counts are set in pharmacophore recognizer
	unsigned char hba;
	unsigned char hbd;
	float logP;
	float psa;

	enum PropIDs {UniqueID, NRings, NAromatics, LogP, PSA, HBA, HBD, None};

	//for querying, need to mmap a file for each of the above
	struct MolPropertyReader
	{
		MMappedRegion<unsigned long> uniqueid;
		MMappedRegion<unsigned char> num_rings;
		MMappedRegion<unsigned char> num_aromatics;
		MMappedRegion<float> logP;
		MMappedRegion<float> psa;
		MMappedRegion<unsigned char> hba;
		MMappedRegion<unsigned char> hbd;

		//access the approprate map at pos and return value casted to double
		double get(PropIDs kind, unsigned pos) const
		{
			switch(kind)
			{
			case UniqueID:
				return uniqueid[pos];
			case NRings:
				return num_rings[pos];
			case NAromatics:
				return num_aromatics[pos];
			case LogP:
				return logP[pos];
			case PSA:
				return psa[pos];
			case HBA:
				return hba[pos];
			case HBD:
				return hbd[pos];
			case None:
				return 0;
			}
			return 0;
		}
	};

	MolProperties(): uniqueid(0), num_rings(0), num_aromatics(0), hba(0), hbd(0), logP(0), psa(0) {}

	void calculate(OpenBabel::OBMol& mol, unsigned long id);

	//set hb counts
	void setHB(const vector<PharmaPoint>& pts);

	typedef vector<FILE*> PropFiles;
	void write(unsigned mid, const PropFiles& files); //write to individual files that were opened with createFiles
	static void createFiles(const boost::filesystem::path& dbpath, PropFiles& files);
	static void initializeReader(const boost::filesystem::path& dbpath, MolPropertyReader& reader);

	static vector<const char*> fileNames;
};

#endif /* MOLPROPERTIES_H_ */
