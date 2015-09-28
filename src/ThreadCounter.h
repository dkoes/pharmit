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
