/*
 * WorkFile2.cpp
 *
 *  Created on: Oct 18, 2011
 *      Author: dkoes
 */

#include "WorkFile.h"
#include <fstream>
#include <iostream>

using namespace boost::interprocess;

WorkFile::WorkFile(const char *name): map(NULL)
{
	file = new ofstream(name);
	assert(*file);
	mapping = new file_mapping(name, read_only);
}

WorkFile::~WorkFile()
{
	//these must be manually cleared
}

void WorkFile::set(const char *name)
{
	clear();
	file = new ofstream(name);
	mapping = new file_mapping(name, read_only);
}

void WorkFile::switchToMap()
{
	if(map == NULL)
	{
		file->close();
		map = new mapped_region(*mapping, read_only);
	}
}

//deallocate and reset
void WorkFile::clear()
{
	if(file) delete file;
	if(mapping) delete mapping;
	if(map) delete map;
	file = NULL;
	mapping = NULL;
	map = NULL;
}

void WorkFile::remove()
{
	if(file == NULL)
		return;
	file->close();
	boost::filesystem::remove(mapping->get_name());

	clear();
}
