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
 * pharminfo.h
 *
 *  Created on: August 11, 2015
 *      Author: dkoes
 *
 *  This file contains routines for reading and writing slightly condensed pharmacophore
 *  features for a whole molecule.  It is used by shape search to filter by pharmacophores.
 */
#include "pharmarec.h"
#include <cstdio>
#include <vector>

//write points to file and return the starting offset
unsigned long writePharmacophoreInfo(FILE *f, const vector<PharmaPoint>& points, const Pharmas& pharmas);

//returns true if each query point matches at least one member of the pharmacophore
//pointed to by pharmacophore - avoid copying into memory
bool pharmacophoreMatchesQuery(const char *pharmacophore, const vector<PharmaPoint>& querypoints, const Pharmas& pharmas);
