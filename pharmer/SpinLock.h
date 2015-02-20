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
 * SpinLock.h
 *
 *  Created on: Feb 10, 2010
 *      Author: dkoes
 *
 *      A very simple spinlock.  Separtes lock var (mutex) from locker so destructor can release lock.
 */

#ifndef SPINLOCK_H_
#define SPINLOCK_H_
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

#endif /* SPINLOCK_H_ */
