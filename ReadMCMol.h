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
 * ReadMCMol.h
 *
 *  Created on: Feb 2, 2011
 *      Author: dkoes
 */

#ifndef READMCMOL_H_
#define READMCMOL_H_
#include <openbabel/mol.h>

using namespace boost;
using namespace std;
using namespace OpenBabel;

//openbabel currently doesn't have native support of multiple conformer input
//this just assumes sequential mols with the same name are conformers
//the read-ahead mol is cached
class ReadMCMol
{
	istream& infile;
	unsigned stride;
	unsigned offset;
	unsigned reduceConfs;
	OBConversion conv;
	unsigned molcnt;

	struct MInfo
	{
		string title;
		OBMol mol;
		string data;

		bool load(OBConversion& conv)
		{
			title.clear();
			mol.Clear();
			data.clear();
			//openbabel is SLOW due to (unneccessary I think) perception routines that
			//get run when construction an obmol from an sdf, so optimize for it separately
			if(strcmp(conv.GetInFormat()->GetID(), "sdf") == 0)
			{
				//parse out a mol ourselves
				istream* in = conv.GetInStream();
				getline(*in, title);
				stringstream d;
				d << title << "\n";
				char c;
				unsigned dollarcnt = 0;
				while(in->get(c))
				{
					d.put(c);
					if(c == '$') dollarcnt++;
					else if(dollarcnt == 4) //readnewline
						break;
					else
						dollarcnt = 0;
				}
				if(dollarcnt != 4)
					return false;
				data = d.str();
				return *in;
			}
			else
			{
				if(!conv.Read(&mol))
					return false;
				title = mol.GetTitle();
			}
			return true;
		}

		OBMol& getMol()
		{
			if(data.length() > 0)
			{
				OBConversion strconv;
				strconv.SetInFormat("sdf");
				strconv.ReadString(&mol, data);
				data.clear();
			}
			return mol;
		}

		bool isValid()
		{
			return title.length() > 0;
		}

	};

	MInfo next;

public:
	ReadMCMol(istream& in, OBFormat* f, unsigned st, unsigned o, unsigned reduce) :
		infile(in), stride(st), offset(o), reduceConfs(reduce > 0 ? reduce : UINT_MAX), molcnt(0)
	{
		conv.SetInFormat(f);
		conv.SetOutFormat("sdf");
		//must do a read to setup conv
		conv.Read(&next.mol, &infile);
		next.title = next.mol.GetTitle();
	}

	//put multi-conformer molecule in mol
	//include weight to ensure value used for filtering matches stored value
	bool read(OBMol& mol)
	{
		mol.Clear();
		if (!next.isValid())
			return false; //nothing else
		while(next.isValid())
		{
			if((molcnt % stride) == offset)
			{
				molcnt++;
				unsigned confcnt = 1;
				//read the rest of conformers for next
				mol = next.getMol();
				unsigned N = mol.NumAtoms();
				string curtitle = next.title;

				while (next.load(conv))
				{
					if (next.title == curtitle)
					{
						if (confcnt < reduceConfs)
						{
							if (next.getMol().NumAtoms() != N)
							{
								cerr
										<< "Invalid Input. Sequential molecules with the same name must be conformers.\n";
								continue; //don't fail, just ignore
							}

							double *confdata = new double[N * 3];
							memcpy(confdata, next.getMol().GetConformer(0),
									sizeof(double) * N * 3);
							mol.AddConformer(confdata); //mol is now managing the memory
						}
						confcnt++;
					}
					else //all done
					{
						break;
					}
				}
				return true;
			}
			else
			{
				//read next mc mol
				molcnt++;
				string curtitle = next.title;
				while(next.load(conv))
				{
					if(next.title != curtitle)
						break;
				}
			}
		}

		return false;
	}

	unsigned molsRead() const { return molcnt; }
};

#endif /* READMCMOL_H_ */
