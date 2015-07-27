/*
 * Timer.h
 *
 *  Created on: Jan 19, 2010
 *      Author: dkoes
 *
 *      My own class for  timing things.  Modeled after boost timer, but
 *      actually measures wall clock time.
 */

#ifndef TIMER_H_
#define TIMER_H_

#include <sys/times.h>
#include <sys/time.h>
#include <cmath>
using namespace std;

class Timer
{
	struct tms start;
	struct timeval startwall;
	const unsigned long ticks;
	//total all the times fields
	clock_t total(const struct tms& x) const
	{
		return x.tms_cstime+x.tms_cutime+x.tms_stime+x.tms_utime;
	}

	//comptue diff in seconds
	double timevalDifF(const struct timeval& start, const struct timeval& end) const
	{
		double s = (double(start.tv_sec) *1000000 + start.tv_usec)/1000000.0;
		double e = (double(end.tv_sec) *1000000 + end.tv_usec)/1000000.0;
		return e - s;
	}

public:
	Timer(): ticks(sysconf(_SC_CLK_TCK))
	{
		times(&start);
		gettimeofday(&startwall, NULL);
	}

	void restart()
	{
		times(&start);
		gettimeofday(&startwall, NULL);
	}

	double elapsed() const // return elapsed wall clock time in seconds
	{
		struct timeval now;
		gettimeofday(&now, NULL);
		//round result for pretty printing..
		return timevalDifF(startwall, now);
	}


	double elapsedProcess() const // return time spent in this process (or children)
	{
		struct tms now;
		times(&now);

		return ((double)(total(now)-total(start))) / ticks;
	}


	double elapsedUser() const //return elapsed user time in seconds for this process
	{
		struct tms now;
		times(&now);

		return (double(now.tms_utime-start.tms_utime))/ticks;
	}

	double elapsedSystem() const //return elapsed system time in seconds for this process
	{
		struct tms now;
		times(&now);

		return (double(now.tms_stime-start.tms_stime))/ticks;
	}
};


#endif /* TIMER_H_ */
