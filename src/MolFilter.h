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
 * MolFilter.h
 *
 *  Created on: Jan 28, 2011
 *      Author: dkoes
 */

#ifndef PHARMITSERVER_MOLFILTER_H_
#define PHARMITSERVER_MOLFILTER_H_

#include <openbabel/mol.h>
using namespace std;

//base class for specifying to creator to skip certain mols
class MolFilter {
public:
	virtual bool skip(OpenBabel::OBMol& mol) = 0;

};

class WeightRangeFilter: public  MolFilter
{
	double min;
	double max;
public:
	WeightRangeFilter():min(0),max(HUGE_VAL) {}

	WeightRangeFilter(double min_, double max_): min(min_), max(max_)
	{

	}

	bool badWeight(double weight)
	{
		if(weight < min)
			return true;
		if(weight >= max)
			return true;
		return false;
	}
	//include min, not max
	bool skip(OpenBabel::OBMol& mol)
	{
		return badWeight(mol.GetMolWt());
	}

};


#endif /* PHARMITSERVER_MOLFILTER_H_ */
