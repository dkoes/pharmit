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
 * ThreadCounter.h
 *
 *  Created on: Aug 2, 2010
 *      Author: dkoes
 */

#ifndef PHARMITSERVER_THREADCOUNTER_H_
#define PHARMITSERVER_THREADCOUNTER_H_

//thread safe counter for running threads (upto hardware supported)
class ThreadCounter
{
	unsigned maxThreads;
	unsigned count;
public:

	//override hardware
	ThreadCounter(unsigned n): maxThreads(n), count(0) {}

	//return true if there's an available thread
	//it is assumed the caller will then spawn a thread and call reduceThreadCount
	bool threadAvailable()
	{
		if(count < maxThreads)
		{
			unsigned cur = __sync_fetch_and_add(&count, 1);
			if(cur < maxThreads)
				return true;
			//someone else beat us to it
			__sync_fetch_and_sub(&count, 1);
			return false;
		}
		else
			return false;
	}

	void reduceThreadCount()
	{
		__sync_fetch_and_sub(&count, 1);
	}
};

#endif /* PHARMITSERVER_THREADCOUNTER_H_ */
