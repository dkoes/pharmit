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
 * QueryTripletFingerprint.h
 *
 *  Created on: Nov 2, 2010
 *      Author: dkoes
 *
 *      Checks to see if a fingerprint potentially matches the query
 */

#ifndef PHARMITSERVER_QUERYTRIPLETFINGERPRINT_H_
#define PHARMITSERVER_QUERYTRIPLETFINGERPRINT_H_

#include "TripletFingerprint.h"
#include "SimpleFingers.h"
#include "BitSetTree.h"
class PharmerQuery;

class QueryTripletFingerprint
{
	//multi resolution search, first check bigqfingers (first DISTPACE), then if there is a match,
	//check corresponding smallqfingers (second DISTPACE)
	vector< vector<TripletFingerprint> > bigqfingers; //all possible fingerprints of each query point not part of this triplet
	vector< vector< vector<TripletFingerprint> > > smallqfingers;
public:
	QueryTripletFingerprint() { assert(TripletFingerprint::NUMDISTSPACES <= 2);}

	 ~QueryTripletFingerprint() {
	 }

	 void set(unsigned i, unsigned j, unsigned k, const PharmerQuery& query);

	 bool isValid(const TripletFingerprint& f) const;

};

#endif /* PHARMITSERVER_QUERYTRIPLETFINGERPRINT_H_ */
