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
 * MTQueue.h
 *
 *  Created on: Jan 20, 2010
 *      Author: dkoes
 *
 *      A very simple multi-threaded queue. All operations are atomic.
 *      Optionally, can support having a limited size where push operations are
 *      blocked until data is consumed.
 *      pop operations block on an empty queue until all data is produced.
 *      Data is copied into the queue.
 */

#ifndef PHARMITSERVER_MTQUEUE_H_
#define PHARMITSERVER_MTQUEUE_H_

#include <boost/thread.hpp>
#include <vector>
#include <queue>
#include <functional>
#include <boost/unordered_map.hpp>
#include "SpinLock.h"

using namespace std;

template <class T>
class MTQueue
{
	deque<T> data;

	SpinMutex mutex;
	volatile unsigned registeredProducers;
	volatile unsigned num; //equals data.size, but volatile for MT
	unsigned maxSize;
public:
	MTQueue(): registeredProducers(0), num(0), maxSize(0)
	{
	}

	//limit quue to at most max entries -- this is approximate though to
	//limit synchronization penalties
	MTQueue(unsigned max): registeredProducers(0), num(0), maxSize(max)
	{

	}

	virtual ~MTQueue() {}

	//register the existence of a producer into this queue
	//only use if multithreaded
	void addProducer()
	{
		__sync_fetch_and_add(&registeredProducers, 1);
	}

	//remove  a producer, indicating no more data
	//only use if multithreaded
	void removeProducer()
	{
		assert(registeredProducers > 0);
		__sync_fetch_and_sub(&registeredProducers, 1);
	}

	unsigned numProducers() { return registeredProducers; }

	//set max allowed elements in Q
	void setMax(unsigned m)
	{
		maxSize = m;
	}
	//return number of elements in queue
	unsigned size()
	{
		return num;
	}

	//push val onto q
	void push(const T& val)
	{
		SpinLock lock(mutex, false);

		if (maxSize > 0) //bounded
		{
			while (true)
			{
				while (num >= maxSize)
				{
					boost::this_thread::yield();
				}
				lock.acquire();
				if (num < maxSize)
					break; //good state
				lock.release(); //keep waiting
			}
		}
		else
		{
			lock.acquire();
		}

		data.push_back(val);
		num++;
	}

	//return true if queue is not empty
	//and fill out val
	//only return false if there is no data available AND there is no chance
	//of there being data available
	bool pop(T& val)
	{
		SpinLock lock(mutex, false);
		//wait for data
		while (true)
		{
			while (num == 0)
			{
				if (registeredProducers == 0 && num == 0) //nothing is producing
					return false;
				boost::this_thread::yield();
			}
			lock.acquire();
			if (num != 0)
				break; //something to pop
			lock.release(); //missed our chance, go aroudn again
		}

		assert(mutex.val == SpinMutexLocked);
		assert(num > 0);
		assert(data.size() > 0);
		val = data.front();
		data.pop_front();
		num--;
		return true;
	}

	//return true if nothing is in the queue nor will there be anything in the queue (unless producer added)
	bool finished()
	{
		if (registeredProducers == 0 && num == 0) //nothing is producing
			return true;
		return false;
	}

	//return true if queue is not empty
	//and fill out val
	//only return false if there is no data available AND there is no chance
	//of there being data available
	//will not block waiting for data
	bool popAll(vector<T>& val)
	{
		SpinLock lock(mutex, false);
		if (num == 0)
		{
			if (finished()) //nothing is producing
				return false;
			return true;
		}

		lock.acquire();

		assert(mutex.val == SpinMutexLocked);
		assert(data.size() > 0);

		val.insert(val.end(), data.begin(), data.end());
		data.clear();
		num = 0;
		return !finished();
	}

};

#endif /* PHARMITSERVER_MTQUEUE_H_ */
