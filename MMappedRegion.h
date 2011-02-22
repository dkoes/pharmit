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
 * MMappedRegion.h
 *
 *  Created on: Aug 2, 2010
 *      Author: dkoes
 */

#ifndef MMAPPEDREGION_H_
#define MMAPPEDREGION_H_
#include <sys/types.h>
#include <sys/mman.h>
#include <iostream>
#include <cstdio>
#include <fcntl.h>
#include <unistd.h>

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

	MMappedRegion(string fname, bool readOnly)
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
	void map(string fname, bool readOnly, bool sequential, bool populate=false, bool readonce=false)
	{
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
			data = (T*) mmap(NULL, size, flags, (readOnly ? MAP_PRIVATE
					: MAP_SHARED) |(populate ? MAP_POPULATE : 0) , fd, 0);
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

#endif /* MMAPPEDREGION_H_ */
