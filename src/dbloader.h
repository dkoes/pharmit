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
		const boost::unordered_map<string, StripedSearchers >& olddatabases);
#endif /* DBLOADER_H_ */
