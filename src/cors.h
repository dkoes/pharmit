/*
Pharmit
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
 * cors.h
 *
 * Classes for maintaining correspondences between point sets.
 * A correspondence just associates one set of indices with another;
 * in this case the indices of the query with the indices of the mol.
 *
 *  Created on: Aug 8, 2010
 *      Author: dkoes
 */

#ifndef PHARMITSERVER_CORS_H_
#define PHARMITSERVER_CORS_H_

#include "pharmerdb.h"
#include "tripletmatching.h"
#include "RMSD.h"

typedef  unsigned char matchType;
//a correspondence that stores data inline making it suitable for bump allocation
class CorrespondenceInline
{
	unsigned char size;
	matchType match[];

	friend class Correspondence;
public:
	//must be allocated propertly so there is space in match
	CorrespondenceInline(unsigned n) :
		size(n)
	{
		clear();
	}

	CorrespondenceInline(const CorrespondenceInline& rhs)
	{
		size = rhs.size; //assume mem allocated
		memcpy(match, rhs.match, size * sizeof(matchType));
	}

	//return index of  point matching the ith query point
	int matchIndex(int i) const
	{
		if(match[i] == UCHAR_MAX) return -1;
		return match[i];
	}

	void print(ostream& out) const
	{
		for (int i = 0; i < size; i++)
		{
			out << match[i];
			if(i < size-1) out << " ";
		}
	}

	//set the  index s of the ith model point
	void setIndex(int i, int s)
	{
		match[i] = s;
	}

	void clear()
	{
		memset(match, -1, sizeof(matchType) * size);
	}

	//write out as a json array
	void writeJSON(ostream& out)
	{
		out << "[";
		for (int i = 0; i < size; i++)
		{
			out << match[i];
			if(i < size-1) out << ",";
		}
		out << "]";
	}

};


//a single correspondence result, must be allocated with a cor allocator
//to properly reserve mem for the inline correspondence
struct CorrespondenceResult
{
	unsigned long location; //incorporates db
	RMSDResult rmsd; //weighted by radii
	float val; //unweighted rmsd val of weighted orientation
	float weight;
	unsigned molid;// also incorporates db - used for hit reduction
	unsigned dbid; //molid for db
	unsigned short nRBnds;
	CorrespondenceInline cor;

	CorrespondenceResult(unsigned numPoints) :
		location(0), val(0), weight(0), molid(0), dbid(0), nRBnds(0), cor(numPoints)
	{
	}

	void reinitialize(const TripletMatch& tm, unsigned db, unsigned numdb)
	{
		location = tm.id*numdb+db;
		weight = tm.weight;
		nRBnds = tm.nRBnds;
		molid = tm.molid*numdb+db;
		dbid = tm.molid;
		rmsd.clear();
		cor.clear();
	}

};

//create CorrespondenceInlines
class CorAllocator
{
	unsigned qsize; //number of points
	unsigned Csize; //size of a CorrespondenceInline
	unsigned CRsize; //size of CorrespondenceResult
	BumpAllocator<1024 * 1024> allocator; //1MB chunks

public:
	CorAllocator() :
		qsize(0), Csize(0), CRsize(0)
	{
	}

	CorAllocator(unsigned qsz) :
		qsize(qsz), Csize(sizeof(CorrespondenceInline) + sizeof(matchType) * qsz),
				CRsize(sizeof(CorrespondenceResult) + sizeof(matchType) * qsz)
	{
	}

	void setSize(unsigned qsz)
	{
		qsize = qsz;
		Csize = sizeof(CorrespondenceInline) + sizeof(matchType) * qsz;
		CRsize = sizeof(CorrespondenceResult) + sizeof(matchType) * qsz;
	}

	CorrespondenceInline* newCor()
	{
		void *ptr = allocator.alloc(Csize);
		return new (ptr) CorrespondenceInline(qsize);
	}

	CorrespondenceResult* newCorResult(const CorrespondenceResult& c)
	{
		void *ptr = allocator.alloc(CRsize);
		return new (ptr) CorrespondenceResult(c);
	}

	CorrespondenceResult* newCorResult()
	{
		void *ptr = allocator.alloc(CRsize);
		return new (ptr) CorrespondenceResult(qsize);
	}

	void clear()
	{
		allocator.clear();
	}

	unsigned numChunks() const { return allocator.numChunks(); }
};



#endif /* PHARMITSERVER_CORS_H_ */
