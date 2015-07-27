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

	const string& getString(unsigned i) const { return strs[i]; }
	double getScore(unsigned i) const { return scores[i]; }
};

#endif /* RESULTS_H_ */
