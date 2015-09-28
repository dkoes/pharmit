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
 * BumpAllocator.h
 *
 *  Created on: Feb 23, 2010
 *      Author: dkoes
 *
 *      A thread-safe bump allocator templatized on the chunk size.
 */

#ifndef PHARMITSERVER_BUMPALLOCATOR_H_
#define PHARMITSERVER_BUMPALLOCATOR_H_

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

#endif /* PHARMITSERVER_BUMPALLOCATOR_H_ */
