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
 * ReadMCMol.h
 *
 *  Created on: Feb 2, 2011
 *      Author: dkoes
 */

#ifndef PHARMITSERVER_READMCMOL_H_
#define PHARMITSERVER_READMCMOL_H_
#include <openbabel/mol.h>
#include <openbabel/generic.h>

typedef OpenBabel::OBPairTemplate< vector<string> > OBVecData;

//openbabel currently doesn't have native support of multiple conformer input
//this just assumes sequential mols with the same name are conformers
//the read-ahead mol is cached
class ReadMCMol
{
	istream& infile;
	unsigned stride;
	unsigned offset;
	unsigned reduceConfs;
	OpenBabel::OBConversion conv;
	unsigned molcnt;
	vector<string> pharmacophores; //accumulate pharmacophore data if available

	struct MInfo
	{
		string title;
		OpenBabel::OBMol mol;
		string data;
		string pharmacophore;
		bool valid;

		MInfo(): valid(false) {}

		bool load(OpenBabel::OBConversion& conv)
		{
			title.clear();
			mol.Clear();
			data.clear();
			valid = false;
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
				valid = bool(*in);
				return bool(*in);
			}
			else
			{
				if(!conv.Read(&mol))
					return false;
				title = mol.GetTitle();
				pharmacophore.clear();
				if(mol.HasData("pharmacophore"))
				{
					pharmacophore = mol.GetData("pharmacophore")->GetValue();
				}
			}
			valid = true;
			return true;
		}

		OpenBabel::OBMol& getMol()
		{
			if(data.length() > 0)
			{
				OpenBabel::OBConversion strconv;
				strconv.SetInFormat("sdf");
				strconv.ReadString(&mol, data);
				data.clear();
				pharmacophore.clear();

				if(mol.HasData("pharmacophore"))
				{
					pharmacophore = mol.GetData("pharmacophore")->GetValue();
				}
			}

			return mol;
		}

		bool isValid()
		{
			return valid;
		}

	};

	MInfo next;

public:

	ReadMCMol(istream& in, OpenBabel::OBFormat* f, unsigned st, unsigned o, unsigned reduce) :
		infile(in), stride(st), offset(o), reduceConfs(reduce > 0 ? reduce : UINT_MAX), molcnt(0)
	{
		conv.SetInFormat(f);
		conv.SetOutFormat("sdf");
		conv.SetInStream(&infile);
		next.load(conv);
	}

	//put multi-conformer molecule in mol
	//include weight to ensure value used for filtering matches stored value
	bool read(OpenBabel::OBMol& mol)
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

				pharmacophores.clear();
				if(next.pharmacophore.size())
					pharmacophores.push_back(next.pharmacophore);

				while (next.load(conv))
				{
					if (next.title == curtitle)
					{
						if (confcnt < reduceConfs)
						{
							if (next.getMol().NumAtoms() != N)
							{
								cerr
										<< "Warning: Invalid Input. Sequential molecules with the same name should be conformers.\n";
								break; //don't fail, just ignore
							}

							double *confdata = new double[N * 3];
							memcpy(confdata, next.getMol().GetConformer(0),
									sizeof(double) * N * 3);
							mol.AddConformer(confdata); //mol is now managing the memory

							if(pharmacophores.size() > 0)
							{
								pharmacophores.push_back(next.pharmacophore);
							}
						}
						confcnt++;
					}
					else //all done
					{
						break;
					}
				}

				if(pharmacophores.size() > 0)
				{
					OBVecData* sddata = new OBVecData();
					sddata->SetAttribute("pharmacophores"); //note the s
					sddata->SetValue(pharmacophores);
					mol.SetData(sddata);
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

#endif /* PHARMITSERVER_READMCMOL_H_ */
