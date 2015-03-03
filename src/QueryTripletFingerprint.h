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
