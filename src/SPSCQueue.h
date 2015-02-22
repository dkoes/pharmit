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
 * SPSCQueue.h
 *
 *  Created on: Feb 11, 2010
 *      Author: dkoes
 *
 *      A single producer/single consumer thread-safe non-blocking queue
 */

#ifndef SPSCQUEUE_H_
#define SPSCQUEUE_H_

using namespace std;
template < typename T, unsigned N>
class SPSCQueue
{
	T* buffer[N];
	//buffering is to avoid false sharing
    volatile unsigned     pread;
    unsigned              padding[15];
    volatile unsigned     pwrite;
    unsigned              padding2[15];
    bool done;
public:
	SPSCQueue(): pread(0), pwrite(0), done(false)
	{
		memset(buffer, 0, sizeof(buffer));
	}
	virtual ~SPSCQueue() {}

	//removes all data, calling delete on the pointers and resets done
	//NOT THREAD SAFE
	void clearAndDestroy()
	{
		for(unsigned i = 0; i < N; i++)
		{
			if(buffer[i])
			{
				delete buffer[i];
				buffer[i] = NULL;
			}
		}
		done = false;
		pread = pwrite = 0;
	}

	//return true if succesful (queue not blocked, val not null)
	bool push(T* val)
	{
        if (val == NULL) return false;

        if (buffer[pwrite] == NULL) {
        	//memory barrier necessary here if stores become visible out of order
        	__sync_synchronize();
            buffer[pwrite] = val;
            pwrite = (pwrite+1)%N;

            return true;
        }
        return false;
	}

	//return true if successful (data currently available)
	bool pop(T*& val)
	{
		if(buffer[pread] == NULL)
			return false;
		val = buffer[pread];
        buffer[pread] = NULL;
        pread = (pread+1)%N;

        return true;
	}

	//empty and will continue to be empty
	bool empty()
	{
		return buffer[pread] == NULL && done;
	}

	//producer has no more data
	void setDone() { done =true; }
};

//stores data, not pointers
template < typename T, unsigned N>
class SPSCDataQueue
{
	struct Node {
		T val;
		bool isEmpty;

		Node(): isEmpty(true) {}
	};

	Node buffer[N];
	//buffering is to avoid false sharing
    volatile unsigned     pread;
    unsigned              padding[15];
    volatile unsigned     pwrite;
    unsigned              padding2[15];
    bool done;
public:
    SPSCDataQueue(): pread(0), pwrite(0), done(false)
	{
	}
	virtual ~SPSCDataQueue() {}

	//rblocks on full queue
	void push(const T& val)
	{
		while(true)
		{
			if (buffer[pwrite].isEmpty)
			{
				//memory barrier necessary here if stores become visible out of order
				__sync_synchronize();
				buffer[pwrite].val = val;
				buffer[pwrite].isEmpty = false;
				if (pwrite >= N - 1)
					pwrite = 0;
				else
					pwrite++;
				return;
			}
			this_thread::yield();
			this_thread::interruption_point();
		}
	}

	//return true if successful (data currently available)
	bool pop(T& val)
	{
		if(buffer[pread].isEmpty)
			return false;
		swap(val,buffer[pread].val);
        buffer[pread].isEmpty = true;;
        if(pread >= N-1)
        	pread = 0;
        else
        	pread++;
        return true;
	}

	T* top()
	{
		if(buffer[pread].isEmpty)
			return NULL;
		return &buffer[pread].val;
	}
	//empty and will continue to be empty
	bool empty()
	{
		return buffer[pread].isEmpty && done;
	}

	//producer has no more data
	void setDone() { done =true; }
};

#endif /* SPSCQUEUE_H_ */
