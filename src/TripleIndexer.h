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
 * TripleIndexer.h
 *
 *  Created on: Aug 3, 2010
 *      Author: dkoes
 *
 *  Class for indexing into a multiset of triples.
 */

#ifndef PHARMITSERVER_TRIPLEINDEXER_H_
#define PHARMITSERVER_TRIPLEINDEXER_H_

#include <boost/multi_array.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <set>

using namespace std;

namespace boost {
static size_t hash_value(const array<unsigned,3>& x)
{
	size_t seed = 0;
	hash_combine(seed, x[0]);
	hash_combine(seed, x[1]);
	hash_combine(seed, x[2]);
	return seed;
}
}

class TripleIndexer
{
	vector< boost::array<unsigned,3> > reverse;
	boost::multi_array<unsigned, 3> lookup; //n*n*n matrix with indices for any ordering
	unsigned sz;
public:

	TripleIndexer() {}

	TripleIndexer(unsigned n): lookup(boost::extents[n][n][n])
	{
		set(n);
	}

	void set(unsigned n)
	{
		using namespace boost;
		lookup.resize(extents[n][n][n]);
		//there are (n+2)!/(3!) sets of triples
		//for fast lookup, we create a table of all n^3 orderings
		//that returns the unique index
		unsigned cnt = math::binomial_coefficient<double>(n+2, 3U);
		reverse.resize(cnt);

		boost::unordered_map< boost::array<unsigned,3>, unsigned > counter; //somewhat inefficient
		for(unsigned i = 0; i < n; i++)
		{
			for(unsigned j = 0; j < n; j++)
			{
				for(unsigned k = 0; k < n; k++)
				{
					boost::array<unsigned,3> ijk = assign::list_of(i)(j)(k);
					sort(ijk.begin(), ijk.end());

					boost::unordered_map< boost::array<unsigned,3>, unsigned >::iterator pos = counter.find(ijk);
					if(pos != counter.end()) //already has index
					{
						lookup[i][j][k] = pos->second;
					}
					else
					{
						unsigned index = counter.size();
						lookup[i][j][k] = index;
						counter[ijk] = index;
						assert(index < reverse.size());
						reverse[index] = ijk;
					}
				}
			}
		}

		//sanity check
		assert(cnt == counter.size());
		sz = cnt;
	}

	unsigned operator()(unsigned i, unsigned j, unsigned k) const
	{
		return lookup[i][j][k];
	}

	void getIJK(unsigned pos, unsigned& i, unsigned& j, unsigned& k) const
	{
		assert(pos < reverse.size());
		i = reverse[pos][0];
		j = reverse[pos][1];
		k = reverse[pos][2];
	}
	unsigned size() const { return sz; }

	void dump(FILE *out) const
	{
		unsigned i,j,k;
		for(unsigned p = 0, n = reverse.size(); p < n; p++)
		{
			getIJK(p, i,j,k);
			fprintf(out, "%d %d %d %d\n", i, j, k, p);
		}

	}
};

#endif /* PHARMITSERVER_TRIPLEINDEXER_H_ */
