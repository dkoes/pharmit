/*
 * ResultMolecules.h
 *
 *  Created on: Mar 20, 2012
 *      Author: dkoes
 *
 *  Container class for reading and storing molecules from the molecular data.
 */

#ifndef RESULTMOLECULES_H_
#define RESULTMOLECULES_H_

#include <iostream>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/unordered_map.hpp>
#include "PMol.h"
#include "Results.h"
using namespace std;
using namespace boost;

class ResultMolecules: public Results
{
	PMolReaderMalloc reader;
	vector<PMol*> mols;
	vector<double> scores;

	bool uniqueify; //true if should only grab one conformer
	unordered_map<string, unsigned > seen; //map titles to position for uniquification
public:
	ResultMolecules(bool uniq=false): uniqueify(uniq)
	{

	}
	virtual ~ResultMolecules()
	{
		clear();
	}

	void clear()
	{
		for(unsigned i = 0, n = mols.size(); i < n; i++)
		{
			if(mols[i]) free(mols[i]);
		}
		mols.clear();
	}

	//add a molecule to the result set; written in pmolf format at data
	//support a single sddata entry for the goodness of the result
	void add(const char *data, double score)
	{
		PMol *mol = (PMol*)reader.readPMol(data);

		if(uniqueify)
		{
			string title = mol->getTitle();
			if(seen.count(title) == 0)
			{
				seen[title] = mols.size();
				mols.push_back(mol);
				scores.push_back(score);
			}
			else
			{
				//update to keep lowest score
				unsigned pos = seen[title];
				if(score < scores[pos])
				{
					free(mols[pos]);
					mols[pos] = mol;
					scores[pos] = score;
				}
				else //delete this one
				{
					free(mol);
				}
			}
		}
		else
		{
			mols.push_back(mol);
			scores.push_back(score);
		}
	}

	void reserve(unsigned n)
	{
		mols.reserve(n);
	}

	//number of molecules in set
	unsigned size() const
	{
		return mols.size();
	}

	//write ith mol
	void writeSDF(ostream& out, unsigned i) const
	{
		assert(i < mols.size());
		vector<ASDDataItem> data;
		data.push_back(ASDDataItem("score",boost::lexical_cast<string>(scores[i])));
		mols[i]->writeSDF(out, data);
	}

	void writeTitleScore(ostream& out, unsigned i) const
	{
		assert(i < mols.size());
		out << mols[i]->getTitle() << "\t" << scores[i] << "\n";
	}

	//write structures to outfile (if defined) and output
	//title/scores to stdout if print is true
	void writeOutput(const std::string outfile, bool print) const
	{
		std::cout.precision(12);
		if(print)
		{
			for (unsigned i = 0, n = size(); i < n; i++)
				writeTitleScore(std::cout, i);
		}

		if(outfile.size() > 0)
		{
			ofstream out(outfile.c_str());
			for (unsigned i = 0, n = size(); i < n; i++)
				writeSDF(out, i);
		}
	}
};

#endif /* RESULTMOLECULES_H_ */
