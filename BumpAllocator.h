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
 * BumpAllocator.h
 *
 *  Created on: Feb 23, 2010
 *      Author: dkoes
 *
 *      A thread-safe bump allocator templatized on the chunk size.
 */

#ifndef BUMPALLOCATOR_H_
#define BUMPALLOCATOR_H_

#include <vector>
#include <boost/foreach.hpp>
#include <iostream>
#include <cstdio>

using namespace std;

template <unsigned ChunkSize, bool ThreadSafe=true>
class BumpAllocator
{
	vector<void *> chunks;
	unsigned lock;
	char * chunk;
	unsigned offset;

public:
	BumpAllocator(): lock(0), offset(0)
	{
		chunks.reserve(1024); //avoid unnecessary resizing
		chunk = (char*)malloc(ChunkSize);
		chunks.push_back(chunk);
	}

	~BumpAllocator()
	{
		clear();
	}

	//free all mem
	void clear()
	{
		BOOST_FOREACH(void *ptr, chunks)
		{
			free(ptr);
		}
		chunks.clear();
		offset = ChunkSize; //force new allocation
	}

	unsigned numChunks() const { return chunks.size(); }

	//allocate size bytes and return
	//thread-safe
	//size must be less than ChunkSize
	void* alloc(unsigned size)
	{
		//spin lock - this function should be very fast
		if(ThreadSafe)
		{ //if we know there are no other threads, don't incur synchronization penalty
			while(__sync_lock_test_and_set(&lock, 1) != 0)
				;
		}
		if(offset+size >= ChunkSize)
		{
			//need new chunk
			chunk = (char*)malloc(ChunkSize);

			chunks.push_back(chunk);
			offset = 0;
		}
		void *ptr = chunk+offset;
		offset += size;

		//release lock
		if(ThreadSafe) {
			__sync_lock_release(&lock, 0);
		}
		return ptr;
	}
};

#endif /* BUMPALLOCATOR_H_ */
