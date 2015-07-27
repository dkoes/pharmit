/*
 * GSSTypes.h
 *
 *  Created on: Oct 18, 2011
 *      Author: dkoes
 *
 *      Simple types used by GSSTrees
 */

#ifndef GSSTYPES_H_
#define GSSTYPES_H_

#define BOOST_FILESYSTEM_VERSION 3

#include <boost/filesystem.hpp>
#include <vector>
#include <cassert>
using namespace std;

typedef unsigned long file_index;

struct result_info
{
	file_index pos; //position of result
	double val; //some measure of goodness of result
	result_info(): pos(0), val(0) {}
	result_info(file_index p, double v): pos(p), val(v) {}

	//sort by position for better access
	bool operator<(const result_info& rhs) const
	{
		return pos < rhs.pos;
	}
};

class DataViewer;
struct Cluster;



#endif /* GSSTYPES_H_ */
