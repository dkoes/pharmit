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
