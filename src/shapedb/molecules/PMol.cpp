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
 * PMol.cpp
 *
 *  Created on: Aug 2, 2010
 *      Author: dkoes
 */

#include "PMol.h"

#include <openbabel/data.h>
#include <openbabel/mol.h>
#include "boost/foreach.hpp"
#include <boost/lexical_cast.hpp>
#include <iomanip>
#include <GraphMol/ROMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>

using namespace RDKit;

//copy data needed to write out pmol from obmol
void PMolCreator::copyFrom(OBMol& mol, bool deleteH)
{
	static OBIsotopeTable isotable;

	name = mol.GetTitle();
	//first construct atoms
	int typeindex[256]; //position in atoms vector of an atom type, indexed by atomic number
	int tmpatomindex[mol.NumAtoms()]; //position within the atom type vector
	memset(typeindex, -1, sizeof(typeindex));
	memset(tmpatomindex, -1, sizeof(tmpatomindex));


	for (OBAtomIterator aitr = mol.BeginAtoms(); aitr != mol.EndAtoms(); ++aitr)
	{
		OBAtom *atom = *aitr;
		unsigned int idx = atom->GetIdx();
		unsigned anum = atom->GetAtomicNum();
		assert(anum < 256);
		if(deleteH && anum == 1)
			continue;
		int pos = typeindex[anum];
		if (pos < 0) // no vector for this type yet
		{
			pos = typeindex[anum] = atoms.size();
			atoms.push_back(AtomGroup(anum));
		}

		tmpatomindex[idx] = atoms[pos].coords.size();
		atoms[pos].coords.push_back(FloatCoord(atom->x(), atom->y(), atom->z()));
	}
	//create mapping to actual atom index
	numAtoms = 0;
	BOOST_FOREACH(AtomGroup& ag, atoms)
	{	ag.startIndex = numAtoms;
		numAtoms += ag.coords.size();
	}

	int atomindex[mol.NumAtoms()]; //position within the atom type vector
	memset(atomindex, -1, sizeof(atomindex));
	for (OBAtomIterator aitr = mol.BeginAtoms(); aitr != mol.EndAtoms(); ++aitr)
	{
		OBAtom *atom = *aitr;
		unsigned int idx = atom->GetIdx();
		unsigned anum = atom->GetAtomicNum();
		if(deleteH && anum == 1)
			continue;
		unsigned ai = tmpatomindex[idx] + atoms[typeindex[anum]].startIndex;
		atomindex[idx] = ai;
		//also add properties
		int charge = atom->GetFormalCharge();
		if(charge != 0)
		{
			chg.push_back(Property(ai, charge));
		}

		unsigned isoi = atom->GetIsotope();
		if(isoi != 0)
		{
			int dmass = ::round(isotable.GetExactMass(isoi)-isotable.GetExactMass(0));
			if(dmass != 0)
			{
				iso.push_back(Property(ai, dmass));
			}
		}
	}

	//bonds are directly indexed by atom index
	for(unsigned i = 0; i < MAX_BONDS; i++)
	{
		bonds[i].resize(numAtoms);
	}

	nSrcs = 0;
	nDsts = 0;
	//construct bounds, always putting atom with highest degree first
	for (OBBondIterator bitr = mol.BeginBonds(); bitr != mol.EndBonds(); ++bitr)
	{
		OBBond *bond = *bitr;
		int btype = bond->GetBondOrder() -1;
		assert(btype >= 0 && btype < MAX_BONDS);
		OBAtom *a1 = bond->GetBeginAtom();
		OBAtom *a2 = bond->GetEndAtom();

		if(deleteH && (a1->GetAtomicNum() == 1 || a2->GetAtomicNum() == 1))
			continue;

		nDsts++; //basically nubmer of bonds
		if(a2->GetValence() > a1->GetValence())
		{
			//a2 goes first
			if(bonds[btype][atomindex[a2->GetIdx()]].size() == 0)
			{
				bndSize[btype] += 2*sizeof(unsigned char); //first time, size of src and size
				nSrcs++;
			}
			bonds[btype][atomindex[a2->GetIdx()]].push_back(atomindex[a1->GetIdx()]);
		}
		else
		{
			//a1 goes first
			if(bonds[btype][atomindex[a1->GetIdx()]].size() == 0)
			{
				bndSize[btype] += 2*sizeof(unsigned char); //first time, size of src and size
				nSrcs++;
			}
			bonds[btype][atomindex[a1->GetIdx()]].push_back(atomindex[a2->GetIdx()]);
		}
		bndSize[btype] += sizeof(unsigned char);
	}

}

//copy data needed to write out pmol from romol
void PMolCreator::copyFrom(ROMol& mol, bool deleteH)
{
	static OBIsotopeTable isotable;
	const Conformer& conf = mol.getConformer();

	mol.getProp("_Name",name);
	//first construct atoms
	int typeindex[256]; //position in atoms vector of an atom type, indexed by atomic number
	int tmpatomindex[mol.getNumAtoms()]; //position within the atom type vector
	memset(typeindex, -1, sizeof(typeindex));
	memset(tmpatomindex, -1, sizeof(tmpatomindex));


	for (ROMol::AtomIterator aitr = mol.beginAtoms(); aitr != mol.endAtoms(); ++aitr)
	{
		Atom *atom = *aitr;
		unsigned int idx = atom->getIdx();
		unsigned anum = atom->getAtomicNum();
		assert(anum < 256);
		if(deleteH && anum == 1)
			continue;
		int pos = typeindex[anum];
		if (pos < 0) // no vector for this type yet
		{
			pos = typeindex[anum] = atoms.size();
			atoms.push_back(AtomGroup(anum));
		}

		tmpatomindex[idx] = atoms[pos].coords.size();
		RDGeom::Point3D pt = conf.getAtomPos(atom->getIdx());
		atoms[pos].coords.push_back(FloatCoord(pt.x,pt.y,pt.z));
	}
	//create mapping to actual atom index
	numAtoms = 0;
	BOOST_FOREACH(AtomGroup& ag, atoms)
	{	ag.startIndex = numAtoms;
		numAtoms += ag.coords.size();
	}

	int atomindex[mol.getNumAtoms()]; //position within the atom type vector
	memset(atomindex, -1, sizeof(atomindex));
	for (ROMol::AtomIterator aitr = mol.beginAtoms(); aitr != mol.endAtoms(); ++aitr)
	{
		Atom *atom = *aitr;
		unsigned int idx = atom->getIdx();
		unsigned anum = atom->getAtomicNum();
		if(deleteH && anum == 1)
			continue;
		unsigned ai = tmpatomindex[idx] + atoms[typeindex[anum]].startIndex;
		atomindex[idx] = ai;
		//also add properties
		int charge = atom->getFormalCharge();
		if(charge != 0)
		{
			chg.push_back(Property(ai, charge));
		}

		int dmass = atom->getIsotope();
		if(dmass != 0)
		{
			iso.push_back(Property(ai, dmass));
		}
	}

	//bonds are directly indexed by atom index
	for(unsigned i = 0; i < MAX_BONDS; i++)
	{
		bonds[i].resize(numAtoms);
	}

	nSrcs = 0;
	nDsts = 0;
	//construct bounds, always putting atom with highest degree first
	for (ROMol::BondIterator bitr = mol.beginBonds(); bitr != mol.endBonds(); ++bitr)
	{
		Bond *bond = *bitr;
		int btype = bond->getBondTypeAsDouble() -1;
		assert(btype >= 0 && btype < MAX_BONDS);
		Atom *a1 = bond->getBeginAtom();
		Atom *a2 = bond->getEndAtom();

		if(deleteH && (a1->getAtomicNum() == 1 || a2->getAtomicNum() == 1))
			continue;

		nDsts++; //basically nubmer of bonds
		if(a2->getExplicitValence() > a1->getExplicitValence())
		{
			//a2 goes first
			if(bonds[btype][atomindex[a2->getIdx()]].size() == 0)
			{
				bndSize[btype] += 2*sizeof(unsigned char); //first time, size of src and size
				nSrcs++;
			}
			bonds[btype][atomindex[a2->getIdx()]].push_back(atomindex[a1->getIdx()]);
		}
		else
		{
			//a1 goes first
			if(bonds[btype][atomindex[a1->getIdx()]].size() == 0)
			{
				bndSize[btype] += 2*sizeof(unsigned char); //first time, size of src and size
				nSrcs++;
			}
			bonds[btype][atomindex[a1->getIdx()]].push_back(atomindex[a2->getIdx()]);
		}
		bndSize[btype] += sizeof(unsigned char);
	}

}

/* Writes a  very succinct representation of a small molecule.
 * Format:
 * total size : for later use, amount of memory needed to read in mol
 * #of atoms
 * #of atom types
 * #of iso properties
 * #of chg properties
 * # src atoms to bonds * MAX_BONDS
 * coordinates
 * [atom type | number of this type] - * #of atom type
 * [atom index | value ] - * # isos
 * [atom index | value ] - * # chgs
 * [srcIndex nDsts [dsts]] - * number of bonds for each bond type
 *
 * Returns false if cannot create molecule due to size constraints
 */
bool PMolCreator::writeBinary(ostream& out) const
{
	PMolHeader header;

	//check sizes
	if(nDsts >= PMOLHEADER_MAX || numAtoms >= PMOLHEADER_MAX) //assume others are < nAtoms
		return false;
	for (unsigned i = 0; i < MAX_BONDS; i++)
	{
		if(bndSize[i] >= PMOLHEADER_MAX)
			return false;
	}

	//calcualte size
	unsigned two = sizeof(unsigned char) * 2;
	unsigned short size = sizeof(header) + sizeof(FloatCoord) * numAtoms + two
			* atoms.size() + two * chg.size() + two * iso.size() + two * nSrcs
			+ sizeof(unsigned char) * nDsts + sizeof(char) * (name.size() + 1);
	unsigned check = 0;
	out.write((char*) &size, sizeof(unsigned short));
	header.nAtoms = numAtoms;
	header.nAtomTypes = atoms.size();
	header.nCHG = chg.size();
	header.nISO = iso.size();
	header.nBnds = nDsts;
	for (unsigned i = 0; i < MAX_BONDS; i++)
	{
		header.adjListSize[i] = bndSize[i]; //used for offsets
	}

	out.write((char*) &header, sizeof(header));
	check += sizeof(header);
	//coordinates
	BOOST_FOREACH(const AtomGroup& ag, atoms)
	{
		out.write((const char*)&ag.coords[0], sizeof(FloatCoord)*ag.coords.size());
		check += sizeof(FloatCoord)*ag.coords.size();
	}

	//number of each atom type
	BOOST_FOREACH(const AtomGroup& ag, atoms)
	{
		out.put(ag.atomic_number);
		out.put(ag.coords.size());
		check += 2;
	}

	//iso
	BOOST_FOREACH(const Property& p, iso)
	{
		out.put(p.atom);
		out.put(p.value);
		check += 2;
	}
	//chg
	BOOST_FOREACH(const Property& p, chg)
	{
		out.put(p.atom);
		out.put(p.value);
		check += 2;
	}

	//bonds
	for(unsigned i = 0; i < MAX_BONDS; i++)
	{
		for(unsigned j = 0, n = bonds[i].size(); j < n; j++)
		{
			if(bonds[i][j].size() > 0)
			{
				out.put(j); //atom index
				out.put(bonds[i][j].size()); //number of dsts
				check += 2;
				BOOST_FOREACH(unsigned char d, bonds[i][j])
				{
					out.put(d); //dsts
					check++;
				}
			}
		}
	}

	//name
	out.write(name.c_str(), name.size()+1); //include null
	check += name.size()+1;

	if(check != size)
	{
		cout << check << " " << size << "\n";
		cout << name << "\n";
		abort();
	}
	return true;
}


void* PMolReaderMalloc::allocate(unsigned size)
{
	return malloc(size);
}

void* PMolReaderSingleAlloc::allocate(unsigned sz)
{
	if(sz >= bsize)
	{
		bsize = 2*sz;
		free(buffer);
		buffer = malloc(bsize);
	}
	return buffer;
}
//allocate a pmol from data and return it
PMol* PMolReader::readPMol(const char *data)
{
	unsigned short size;
	memcpy(&size, data, sizeof(size));
	void *mem = allocate(sizeof(PMol) + size - sizeof(PMolHeader));

	PMol *ret = new (mem) PMol();
	memcpy(&ret->header, data + sizeof(size), size);
	unsigned off = ret->setup();
	assert(size == off+sizeof(PMolHeader));
	return ret;
}

//allocate a pmol from f and return it
PMol* PMolReader::readPMol(FILE *f)
{
	unsigned short size;
	int check = fread(&size, sizeof(size), 1, f);
	assert(check == 1);
	void *mem = allocate(sizeof(PMol) + size - sizeof(PMolHeader));

	PMol *ret = new (mem) PMol();
	check = fread(&ret->header, size, 1, f);
	assert(check == 1);
	unsigned off = ret->setup();
	assert(size == off+sizeof(PMolHeader));
	return ret;
}

//data is already initialized, just setup pointers
unsigned PMol::setup()
{
	unsigned offset = sizeof(FloatCoord) * header.nAtoms;
	atomtypes = (AtomTypeCnts*) &buffer[offset];

	offset += header.nAtomTypes * sizeof(AtomTypeCnts);
	iso = (Property*) &buffer[offset];
	offset += header.nISO * sizeof(Property);
	chg = (Property*) &buffer[offset];

	offset += header.nCHG * sizeof(Property);

	for (unsigned i = 0; i < MAX_BONDS; i++)
	{
		adjlists[i] = &buffer[offset];
		offset += header.adjListSize[i];
	}
	name = &buffer[offset];
	offset += strlen(name);
	return offset + 1; //null terminator
}

//return molecular weight, adjusts average weight by iso,
//but does not use weights specific to isomeric atoms
double PMol::getMolWeight() const
{
	double ret = 0;
	for (unsigned i = 0, n = header.nAtomTypes; i < n; i++)
	{
		ret += etab.GetMass(atomtypes[i].atomic_number)
				* atomtypes[i].cnt;
	}

	for (unsigned i = 0, n = header.nISO; i < n; i++)
	{
		ret += iso[i].value;
	}

	return ret;
}


//write sdf with associated meta data
void PMol::writeSDF(ostream& out, const vector<ASDDataItem>& sddata,
		const RMSDResult& rms)
{
	vector<FloatCoord> transformed;
	FloatCoord *coords = header.coords;
	if (rms.value() > 0)
	{
		//transform a copy of the points
		transformed.insert(transformed.begin(), coords, coords + header.nAtoms);
		coords = &transformed[0];
		rms.reorient(transformed.size(), (float*)coords);
	}

	//line 1 - name
	out << name << "\n";
	//line 2 - blank
	out << "\n";
	//line 3 - blank
	out << "\n";
	//connection table
	//number of atoms, bonds,  etc, fixed width

	const int BUFFSIZE = 256;
	char buff[BUFFSIZE];

	snprintf(buff, BUFFSIZE, "%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d V2000\n",
			header.nAtoms, header.nBnds, 0, 0, 0, 0, 0, 0, 0, 0, 999);
	out << buff;

	//atom block
	//internal atom numbers are indexed from zero, sdf number from one
	unsigned curType = 0;
	unsigned typecnt = 0;
	for (unsigned ai = 0; ai < header.nAtoms; ai++)
	{
		snprintf(buff, BUFFSIZE, "%10.4f%10.4f%10.4f %-3s%2d%3d%3d%3d%3d\n",
				coords[ai].x, coords[ai].y, coords[ai].z, etab.GetSymbol(
						atomtypes[curType].atomic_number), 0, 0, 0, 0, 0);
		out << buff;
		typecnt++;
		if (typecnt >= atomtypes[curType].cnt)
		{
			//go to next type
			curType++;
			typecnt = 0;
		}
	}

	//bond list
	for (unsigned i = 0; i < MAX_BONDS; i++)
	{
		unsigned curpos = 0;
		unsigned maxpos = header.adjListSize[i];
		while (curpos < maxpos)
		{
			AdjList *list = (AdjList*) (adjlists[i] + curpos);
			for (unsigned k = 0; k < list->nDsts; k++)
			{
				snprintf(buff, BUFFSIZE, "%3d%3d%3d%3d%3d%3d\n",
						list->src + 1, list->dsts[k] + 1, i + 1, 0, 0, 0);
				out << buff;
			}
			curpos += sizeof(AdjList) + sizeof(unsigned char) * list->nDsts;
		}
		assert(curpos == maxpos);
	}

	//chg properties
	if (header.nCHG > 0)
	{
		//output max 8 per line
		int counter = 0;
		for (unsigned i = 0, n = header.nCHG; i < n; i++, counter++)
		{
			if (counter % 8 == 0)
			{
				if (counter > 0)
					out << "\n";
				//number of entries
				out << "M  CHG" << setw(3) << min(n - counter, 8U);
			}
			out << setw(4) << (int)chg[i].atom+1 << setw(4) << (int)chg[i].value;
		}
		out << "\n";
	}
	//iso props
	if (header.nISO > 0)
	{
		//output max 8 per line
		int counter = 0;
		for (unsigned i = 0, n = header.nISO; i < n; i++, counter++)
		{
			if (counter % 8 == 0)
			{
				if (counter > 0)
					out << "\n";
				//number of entries
				out << "M  ISO" << setw(3) << min(n - counter, 8U);
			}
			out << setw(4) << (int)iso[i].atom+1 << setw(4) << (int)iso[i].value;
		}
		out << "\n";
	}
	out << "M  END\n";

	//output sd data
	BOOST_FOREACH(const ASDDataItem& data, sddata)
	{	out << "> <" << data.tag << ">\n";
	out << boost::lexical_cast<string>(data.value) << "\n\n";
	}
	out << "$$$$" << endl;
}


