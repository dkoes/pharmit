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
 * dbloader.h
 *
 * Helper routines for quickly loading databases.
 *  Created on: Mar 11, 2015
 *      Author: dkoes
 */

#ifndef DBLOADER_H_
#define DBLOADER_H_

#include "PharmerQuery.h"
#include <boost/filesystem.hpp>
#include <boost/unordered_map.hpp>


//fill in databases present in dbpaths
void loadDatabases(vector<boost::filesystem::path>& dbpaths, StripedSearchers& databases);

//populated databases using prefixes
void loadFromPrefixes(vector<boost::filesystem::path>& prefixes, boost::unordered_map<string, StripedSearchers >& databases);

void loadNewFromPrefixes(vector<boost::filesystem::path>& prefixes,
		boost::unordered_map<string, StripedSearchers >& databases,
		const boost::unordered_map<string, StripedSearchers >& olddatabases, bool deactivate=false);
#endif /* DBLOADER_H_ */
