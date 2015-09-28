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
 * SPSCQueue.h
 *
 *  Created on: Feb 11, 2010
 *      Author: dkoes
 *
 *      A single producer/single consumer thread-safe non-blocking queue
 */

#ifndef PHARMITSERVER_SPSCQUEUE_H_
#define PHARMITSERVER_SPSCQUEUE_H_

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

#endif /* PHARMITSERVER_SPSCQUEUE_H_ */
