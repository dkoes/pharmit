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
 * SpinLock.h
 *
 *  Created on: Feb 10, 2010
 *      Author: dkoes
 *
 *      A very simple spinlock.  Separtes lock var (mutex) from locker so destructor can release lock.
 */

#ifndef PHARMITSERVER_SPINLOCK_H_
#define PHARMITSERVER_SPINLOCK_H_
#include <cassert>

enum SpinMutexEnum {SpinMutexUnlocked, SpinMutexLocked};

//make class for automatic initialization
struct SpinMutex
{
	struct {
		long prevent_false_sharing[7];
		volatile SpinMutexEnum val;
		long prevent_false_sharing_after[8];
	};
	SpinMutex(): val(SpinMutexUnlocked) {}
	SpinMutex(SpinMutexEnum v): val(v) {}
};

class SpinLock
{
	SpinMutex& mutex;
	bool holdsLock; //true if this lock has the mutex
public:

	static void acquire(SpinMutex& m)
	{
		while(!__sync_bool_compare_and_swap (&m.val, SpinMutexUnlocked, SpinMutexLocked))
			; //spin if was already locked
		assert(m.val == SpinMutexLocked);
	}

	static void release(SpinMutex& m)
	{
		assert(m.val == SpinMutexLocked);
		m.val = SpinMutexUnlocked;
	}

	void acquire() { acquire(mutex); holdsLock = true; }
	void release() { if(holdsLock) release(mutex); holdsLock = false;}

	//acquire mutex
	SpinLock(SpinMutex& m, bool lock=true): mutex(m)
	{
		if(lock)
			acquire();
		holdsLock = lock;
	}

	virtual ~SpinLock()
	{
		//always release if we're holding the lock
		if(holdsLock)
			release();
	}
};

#endif /* PHARMITSERVER_SPINLOCK_H_ */
