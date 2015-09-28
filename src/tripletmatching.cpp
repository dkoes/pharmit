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
 * tripletmatching.cpp
 *
 *  Created on: Aug 6, 2010
 *      Author: dkoes
 *
 *      Implementations of triplet matching.
 */
#include "tripletmatching.h"

unsigned long const TripletMatchHash::prime_list[] =
{ 53ul, 97ul, 193ul, 389ul, 769ul, 1543ul, 3079ul, 6151ul, 12289ul, 24593ul,
		49157ul, 98317ul, 196613ul, 393241ul, 786433ul, 1572869ul, 3145739ul,
		6291469ul, 12582917ul, 25165843ul, 50331653ul, 100663319ul,
		201326611ul, 402653189ul, 805306457ul, 1610612741ul, 3221225473ul,
		4294967291ul };
const unsigned TripletMatchHash::endPrimeIndex =
		sizeof(TripletMatchHash::prime_list) / sizeof(unsigned long);

TripletMatchAllocator::TripletMatchAllocator(unsigned qsz) :
	qsize(qsz), PMsize(8 * (((sizeof(TripletMatch) + qsz
			* sizeof(TripletMatchInfoArray)) - 1) / 8 + 1)) //align to 8 byte boudnary
{

}

TripletMatch* TripletMatchAllocator::newTripletMatch(unsigned mid,
		const ThreePointData& tdata)
{
	void *ptr = allocator.alloc(PMsize);
	return new (ptr) TripletMatch(mid, tdata, qsize);
}

unsigned TripletMatchInfoArray::allocSize(unsigned numEl)
{
	return sizeof(TripletMatchInfoArray) + (numEl-TMINFO_ARRAY_STARTSIZE)*sizeof(TripletMatchInfo);
}

TripletMatchInfoArray* TripletMatchAllocator::newTripletMatchInfoArray(unsigned numEl)
{
	TripletMatchInfoArray *ptr = (TripletMatchInfoArray*)allocator.alloc(TripletMatchInfoArray::allocSize(numEl));
	//avoid constructor to avoid initializing all the array contents
	ptr->next = NULL;
	ptr->num = 0;
	return ptr;
}

//double the size of the table
void TripletMatchHash::grow()
{
	TripletMatch **oldtable = table;
	unsigned long oldtable_size = table_size;
	prime_index+=2;
	if(prime_index >= endPrimeIndex)
		table_size = table_size*2+1; //off the end
	else
		table_size = prime_list[prime_index];

	table = (TripletMatch**)malloc(table_size*sizeof(TripletMatch*));
	memset(table, 0, sizeof(TripletMatch*)*table_size);

	//rehash the old table
	for(unsigned long i = 0; i < oldtable_size; i++)
	{
		if(oldtable[i] != NULL)
		{
			TripletMatch *m = oldtable[i];
			unsigned long pos = getPos(m->id);
			table[pos] = m;
		}
	}

	free(oldtable);
}
