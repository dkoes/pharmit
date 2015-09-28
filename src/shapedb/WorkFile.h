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
 * WorkFile.h
 *
 *  Created on: Oct 18, 2011
 *      Author: dkoes
 *
 *      Simple wrapper for a file that is written to and then mapped.
 */

#ifndef WORKFILE_H_
#define WORKFILE_H_
#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/filesystem.hpp>

using namespace std;

//store infor for files that are being created and will be memory mapped
struct WorkFile
{
	//none of these has a real copy constructor
	ofstream *file;
	boost::interprocess::mapped_region *map;
	boost::interprocess::file_mapping *mapping;

	WorkFile(): file(NULL), map(NULL), mapping(NULL) {}
	WorkFile(const char *name);
	~WorkFile();

	void switchToMap();
	void set(const char *name);
	void clear();
	void remove(); //removes file

	friend void swap(WorkFile& a, WorkFile& b)
	{
		swap(a.file,b.file);
		swap(a.map, b.map);
		swap(a.mapping, b.mapping);
	}
};



#endif /* WORKFILE_H_ */
