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
 * Results.h
 *
 *  Created on: Jun 3, 2013
 *      Author: dkoes
 *
 *  Base class for results containers. Default behavior is to do nothing at all.
 */

#ifndef RESULTS_H_
#define RESULTS_H_

#include <vector>
#include <string>

class Results
{
public:
	Results()
	{
	}
	virtual ~Results()
	{
	}

	virtual void clear()
	{
	}

	virtual void add(const char *data, double score) = 0;

	virtual void reserve(unsigned n)
	{
	}
	virtual unsigned size() const
	{
		return 0;
	}

	virtual bool stopEarly() const { return false; }
};

//for ojects that just store a string identifier (null terminated)
class StringResults: public Results
{
	std::vector<std::string> strs;
	std::vector<double> scores;
	public:
	StringResults()
	{
	}
	virtual ~StringResults()
	{
	}

	virtual void clear()
	{
		strs.clear();
	}
	virtual void add(const char *data, double score)
	{
		strs.push_back(data); //better be null terminated
		scores.push_back(score);
	}
	virtual void reserve(unsigned n)
	{
		strs.reserve(n);
		scores.reserve(n);
	}
	virtual unsigned size() const
	{
		return strs.size();
	}

	const std::string& getString(unsigned i) const { return strs[i]; }
	double getScore(unsigned i) const { return scores[i]; }
};

#endif /* RESULTS_H_ */
