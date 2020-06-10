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
 * MMappedRegion.h
 *
 *  Created on: Aug 2, 2010
 *      Author: dkoes
 */

#ifndef PHARMITSERVER_MMAPPEDREGION_H_
#define PHARMITSERVER_MMAPPEDREGION_H_
#include <sys/types.h>
#include <sys/mman.h>
#include <iostream>
#include <cstdio>
#include <fcntl.h>
#include <unistd.h>
#include <boost/filesystem.hpp>
#include <string>

using namespace std;

//a memory mapped file viewed as an array of T
template<class T>
class MMappedRegion
{
	T *data;
	unsigned long long size;

	//sorry, can't copy these
	MMappedRegion(const MMappedRegion& rhs)
	{
		abort();
	}
public:
	MMappedRegion() :
		data(NULL), size(0)
	{
	}

	MMappedRegion(const string& fname, bool readOnly)
	{
		map(fname, readOnly);
	}

	~MMappedRegion()
	{
		//free and sync back to disk
		munmap(data, size);
	}

	//creating mapping to filename fname
	//readOnly should be set if only for reading
	void map(const string& fname, bool readOnly, bool sequential=true, bool populate=false, bool readonce=false)
	{
		using namespace boost;
		if (data != NULL)
			munmap(data, size);
		unsigned flags = readOnly ? O_RDONLY : O_RDWR;
		int fd = open(fname.c_str(), flags);
		if(fd < 0)
		{
			perror(fname.c_str());
			abort();
		};
#ifdef POSIX_FADV_SEQUENTIAL
		if(sequential)
			posix_fadvise(fd, 0,0,POSIX_FADV_SEQUENTIAL);
#endif
#ifdef POSIX_FADV_NOREUSE
		if(readonce)
			posix_fadvise(fd, 0, 0, POSIX_FADV_NOREUSE);
#endif
		flags = readOnly ? PROT_READ : (PROT_READ | PROT_WRITE);
		if (filesystem::file_size(fname) > 0)
		{
			size = filesystem::file_size(fname);
			data = (T*) mmap(NULL, size, flags, (readOnly ? MAP_PRIVATE|MAP_NORESERVE
					: MAP_SHARED) |(populate ? MAP_POPULATE : 0), 
					fd, 0);
			close(fd); //once mmapped, don't need to keep open             
		}
		else
		{
			size = 0;
			data = NULL;
		}

		if ((long) data == -1)
		{
			perror("mmap  failed: ");
			abort();
		}

	}
	T& operator[](unsigned long long i)
	{
		return data[i];
	}
	const T& operator[](unsigned long long i) const
	{
		return data[i];
	}
	unsigned long long length() const
	{
		assert(size % sizeof(T) == 0);
		return size / sizeof(T);
	}
	T* begin()
	{
		return data;
	}
	T* end()
	{
		return data + length();
	}
	void sync() //flush back to disk
	{
		msync(data, size, MS_SYNC);
	}
};

#endif /* PHARMITSERVER_MMAPPEDREGION_H_ */
