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
