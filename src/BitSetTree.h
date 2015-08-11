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
 * BitSetTree.h
 *
 *  A bitset tree.  A balanced tree that organizes bitsets, in this case
 *  triplet fingerprints.  Right now this is expermental and assumes
 *  everything can fit in memory.
 *  Created on: Jan 4, 2011
 *      Author: dkoes
 */

#ifndef PHARMITSERVER_BITSETTREE_H_
#define PHARMITSERVER_BITSETTREE_H_

#include "TripletFingerprint.h"
#include "QueryTripletFingerprint.h"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>
#include "Timer.h"
#include <boost/foreach.hpp>

class BitSetTree
{
	struct Node {
		TripletFingerprint ored;
		TripletFingerprint anded;
		unsigned left;
		unsigned right;
		unsigned short dist;
		Node (): left(0), right(0), dist(0) {}

		Node(const TripletFingerprint& v): ored(v), anded(v), left(0), right(0), dist(0) {}
		Node(unsigned l, unsigned r, const Node& lc, const Node& rc):
			ored(lc.ored | rc.ored), anded(lc.anded & rc.anded), left(l), right(r), dist(max(lc.dist, rc.dist)+1) {}
	};

	vector <Node> tree;
	unsigned root;
	unsigned levels;
	vector <unsigned> nonzeroRoots;
	vector<TripletFingerprint> fingerHeap;
	unsigned heapSkip;

	//distance - size of union - want to minimize
	int distanceGBU(const Node& a, const Node& b)
	{
		TripletFingerprint u = a.ored | b.ored;
		TripletFingerprint i = a.anded & b.anded;
		u ^= i;

		return u.bitcnt();
	}

	//create an upper level for the nodes betweeen start and end
	void joinGreedyBU(unsigned start, unsigned end)
	{
		unsigned len = end-start;
		vector<bool> merged(len, false);

		//n^2
		for(unsigned i = 0; i < len; i++)
		{
			if(!merged[i])
			{
				//find most similar
				unsigned bestpos = 0;
				int min = INT_MAX;

				for(unsigned j = i+1; j < len; j++)
				{
					if(!merged[j])
					{
						int d = distanceGBU(tree[start+i], tree[start+j]);
						if(d < min)
						{
							min = d;
							bestpos = j;
						}
					}
				}

				if(bestpos == 0) //didn't find one - only one left
				{
					tree.push_back(Node(start+i, 0, tree[start+i], tree[start+i]));
				}
				else
				{
					tree.push_back(Node(start+i,start+bestpos,tree[start+i],tree[start+bestpos]));
					merged[bestpos] = true;
				}
				merged[i] = true;
			}
		}
	}

	//this is the algorithm described in the paper
	//randomly (first) select a finger, find the closest and combine
	void constructGreedyBU(const vector<TripletFingerprint>& fingers)
	{
		Timer t;
		if(fingers.size() == 0) return;
	//	cout << "Building bitset tree\n";
		tree.reserve(fingers.size()*2+1);
		tree.resize(1); //blank node at index zero
		for(unsigned i = 0, n = fingers.size(); i < n; i++)
		{
			tree.push_back(Node(fingers[i]));
		}

		unsigned curstart = 1;
		unsigned curend = tree.size();
		levels = 1;

		while(curend-curstart > 1)
		{
			//cout << " Level " << levels << " (" << curend-curstart << ")\n";
			joinGreedyBU(curstart, curend);
			curstart = curend;
			curend = tree.size();
			levels++;
		}

		root = tree.size()-1;

		//cout << "Built bitset tree " << t.elapsed() << "\n";
	}

	void cntbits(const vector<TripletFingerprint>& fingers, vector<unsigned>& bitcnts)
	{
		bitcnts.assign(256, 0);

		for(unsigned f = 0, n = fingers.size(); f < n; f++)
		{
			for(unsigned i = 0; i < 256; i++)
			{
				if(fingers[f].getBit(i))
					bitcnts[i]++;
			}
		}

	}


	//divide based on bits
	void split(const vector<TripletFingerprint>& fingers, vector<TripletFingerprint>& left, vector<TripletFingerprint>& right)
	{
		vector<unsigned> bitcnts;

		left = fingers;
		right.clear();
		unsigned half = (1+fingers.size()) / 2;

		//peel off a minimally informative bit at a time in attempt to
		//maximize the number of distinguishing bits at each level
		while (true)
		{
			cntbits(left, bitcnts);
			unsigned maxset = 0;
			unsigned maxsetpos = 0;
			unsigned maxunset = 0;
			unsigned maxunsetpos = 0;
			for (unsigned i = 0, n = bitcnts.size(); i < n; i++)
			{
				if (bitcnts[i] > maxset && bitcnts[i] != left.size())
				{
					maxset = bitcnts[i];
					maxsetpos = i;
				}
				if (left.size() - bitcnts[i] > maxunset && bitcnts[i] != 0)
				{
					maxunset = left.size() - bitcnts[i];
					maxunsetpos = i;
				}
			}
//cout << " left size " << left.size() <<" maxunset " << maxunset << " pos " << maxunsetpos << "\n";
			if(false && maxset >= maxunset)
			{
				//all set bits stay in left
				vector<TripletFingerprint> newleft;
				for(unsigned i = 0, n = left.size(); i < n; i++)
				{
					if(left[i].getBit(maxsetpos))
						newleft.push_back(left[i]);
					else
					{
						//only do a partial split if necessary to do a binary division
						if(right.size() >= half)
							newleft.push_back(left[i]);
						else
							right.push_back(left[i]);
					}
				}
				swap(left, newleft);
			}
			else
			{
				//all unset bits stay in left
				vector<TripletFingerprint> newleft;
				for(unsigned i = 0, n = left.size(); i < n; i++)
				{
					if(maxunset > 0 && !left[i].getBit(maxunsetpos))
						newleft.push_back(left[i]);
					else
					{
						if(right.size() >= half)
							newleft.push_back(left[i]);
						else
							right.push_back(left[i]);
					}
				}
				swap(left, newleft);
			}

			if(left.size() == half || right.size() == half)
				break;
		}
	}

	//create a node for fingers, split into children and recurse
	unsigned createPeelNode(const vector<TripletFingerprint>& fingers, unsigned l)
	{
		if(fingers.size() == 0)
			return 0;

		if(l > levels)
			levels = l;
		unsigned ret = tree.size();
		tree.push_back(Node());

		tree[ret].anded.setAll();
		for(unsigned i = 0, n = fingers.size(); i < n; i++)
		{
			tree[ret].anded &= fingers[i];
			tree[ret].ored |= fingers[i];
		}

		if(l == 10) cout << "LevelSize " << l << " " << tree[ret].ored.bitcnt() << "\n";

		if(tree[ret].anded == tree[ret].ored)
		{
			tree[ret].left = 0;
			tree[ret].right = 0;
		}
		else if(fingers.size() == 1)
		{
			tree[ret].left = 0;
			tree[ret].right = 0;
		}
		else
		{
			vector<TripletFingerprint> left;
			vector<TripletFingerprint> right;

			split(fingers, left, right);

			tree[ret].left = createPeelNode(left, l+1);
			left.clear();
			tree[ret].right = createPeelNode(right, l+1);
		}
		return ret;
	}

	//peel off most discrimative bits
	void construectPeelTD(const vector<TripletFingerprint>& fingers)
	{
		Timer t;
		cout << "Building bitset tree\n";
		tree.reserve(fingers.size()*2+1);
		tree.resize(1); //blank node at index zero
		levels = 0;
		root = createPeelNode(fingers, 1);
		cout << "Built bitset tree " << t.elapsed() << "\n";

	}

	typedef Eigen::Matrix<float, 256, 256> M;
	typedef Eigen::Matrix<float, 256, 1> V;

	struct FInfo
	{
		double val;
		unsigned pos;

		FInfo(): val(0), pos(0) {}

		bool operator<(const FInfo& rhs) const
		{
			return val < rhs.val;
		}
	};

	void setBits(const TripletFingerprint& f, V& bits)
	{
		for(unsigned b = 0; b < 256; b++)
		{
			bits[b] = f.getBit(b);
		}
	}

	//compute the first principal component, and divide evenly
	void splitpca(const vector<TripletFingerprint>& fingers, vector<TripletFingerprint>& left, vector<TripletFingerprint>& right)
	{
		M c = M::Zero();
		V bcnt = V::Zero();
		unsigned N = fingers.size();

		for(unsigned i = 0; i < N; i++)
		{
			const TripletFingerprint& f = fingers[i];
			V bits;
			setBits(f, bits);

			bcnt += bits;
			c += bits * bits.transpose();
		}

		V mean = bcnt / (double)N;
		M mc = bcnt * mean.transpose();
		M C = c - mc - mc.transpose() + (N * (mean * mean.transpose()));
		C /= (N-1);

		Eigen::SelfAdjointEigenSolver<M> eigensolver(C);
		//we'll assume the largest is the last
		V evec = eigensolver.eigenvectors().col(eigensolver.eigenvectors().cols()-1);

		vector<FInfo> finfos(N);
		for(unsigned i = 0; i < N; i++)
		{
			finfos[i].pos = i;
			V bits;
			setBits(fingers[i], bits);
			finfos[i].val = evec.transpose() * bits;
		}

		sort(finfos.begin(), finfos.end());

		left.clear(); left.reserve(N/2+1);
		right.clear(); right.reserve(N/2+1);

		unsigned i;
		for(i = 0; i < N/2; i++)
		{
			left.push_back(fingers[finfos[i].pos]);
		}
		for( ; i < N; i++)
		{
			right.push_back(fingers[finfos[i].pos]);
		}
	}

	//create a node for fingers, split into children and recurse (using PCA!)
	unsigned createPCANode(const vector<TripletFingerprint>& fingers, unsigned l)
	{
		if(fingers.size() == 0)
			return 0;

		if(l > levels)
			levels = l;
		unsigned ret = tree.size();
		tree.push_back(Node());

		tree[ret].anded.setAll();
		for(unsigned i = 0, n = fingers.size(); i < n; i++)
		{
			tree[ret].anded &= fingers[i];
			tree[ret].ored |= fingers[i];
		}

		if(l == 10) cout << "LevelSize " << l << " " << tree[ret].ored.bitcnt() << "\n";

		if(tree[ret].anded == tree[ret].ored)
		{
			tree[ret].left = 0;
			tree[ret].right = 0;
		}
		else if(fingers.size() == 1)
		{
			tree[ret].left = 0;
			tree[ret].right = 0;
		}
		else
		{
			vector<TripletFingerprint> left;
			vector<TripletFingerprint> right;

			splitpca(fingers, left, right);

			tree[ret].left = createPCANode(left, l+1);
			left.clear();
			tree[ret].right = createPCANode(right, l+1);
		}
		return ret;
	}

	void constructPCA(const vector<TripletFingerprint>& fingers)
	{
		Timer t;
		cout << "Building bitset tree\n";
		tree.reserve(fingers.size()*2+1);
		tree.resize(1); //blank node at index zero
		levels = 0;
		root = createPCANode(fingers, 1);
		cout << "Built bitset tree " << t.elapsed() << "\n";
	}

	void computeBitCnts(vector<vector<int> >& levelbitcnts, const Node& N, unsigned level)
	{
		if(levelbitcnts.size() <= level)
			levelbitcnts.resize(level+1);

		int b = N.ored.bitcnt();
		if(N.left == 0 && N.right == 0)
			b = -b;
		levelbitcnts[level].push_back(b);

		if(N.left > 0)
			computeBitCnts(levelbitcnts, tree[N.left], level+1);
		if(N.right > 0)
			computeBitCnts(levelbitcnts, tree[N.right], level+1);
	}

	bool hasContainedIn(const Node& N, const TripletFingerprint& bset) const
	{
		if(bset.contains(N.anded)) //possibly true
		{
			if(N.left == 0 && N.right == 0)
				return true; //a leaf
			if(N.left > 0 && hasContainedIn(tree[N.left], bset))
				return true;
			if(N.right > 0 && hasContainedIn(tree[N.right], bset))
				return true;
			return false;
		}
		return false;
	}

	//first the first nodes with non-zero andeds
	void addNonZeroRoots(unsigned n)
	{
		if(n == 0) return;
		if(tree[n].anded.isZero())
		{
			addNonZeroRoots(tree[n].left);
			addNonZeroRoots(tree[n].right);
		}
		else
		{
			nonzeroRoots.push_back(n);
		}
	}

	void createMiniHeap(unsigned n, vector<TripletFingerprint>& heap, unsigned pos, unsigned max)
	{
		if(pos >= max)
		{
			assert(n == 0);
			return;
		}

		if(n == 0)
		{
			heap[pos].setAll();
			createMiniHeap(0, heap, 2*pos, max);
			createMiniHeap(0, heap, 2*pos+1, max);
		}
		else
		{
			heap[pos] = tree[n].anded;
			createMiniHeap(tree[n].left, heap, 2*pos, max);
			createMiniHeap(tree[n].right, heap, 2*pos+1, max);
		}
	}
	//create a set of fixed sized heaps off of each nonzero internal node
	void createFingerHeap_r(unsigned n, unsigned depth)
	{
		if(n == 0) return;
		if(tree[n].dist > depth)
		{
			createFingerHeap_r(tree[n].left, depth);
			createFingerHeap_r(tree[n].right, depth);
		}
		else if(tree[n].dist == depth)
		{
			vector<TripletFingerprint> miniheap(heapSkip);
			createMiniHeap(n, miniheap, 1, heapSkip);
			fingerHeap.insert(fingerHeap.end(), miniheap.begin(), miniheap.end());
		}
		else
			abort();
	}

	//make sure dists are right
	void computeDists(unsigned n)
	{
		if(n == 0) return;
		 if(tree[n].left == 0 && tree[n].right == 0)
			 tree[n].dist = 0;
		 else if(tree[n].left == 0)
		 {
			 computeDists(tree[n].right);
			 tree[n].dist = tree[tree[n].right].dist+1;
		 }
		 else if(tree[n].right == 0)
		 {
			 computeDists(tree[n].left);
			 tree[n].dist = tree[tree[n].left].dist+1;
		 }
		 else
		 {
			 computeDists(tree[n].left);
			 computeDists(tree[n].right);
			 tree[n].dist = max(tree[tree[n].left].dist, tree[tree[n].right].dist)+1;
		 }
	}

	//check miniheap (not really a heap, I know) at possition off
	bool containedInMiniHeap(unsigned start, unsigned off, const TripletFingerprint& bset) const
	{
		if(off >= heapSkip)
			return true;
		if(bset.contains(fingerHeap[start+off]))
		{
			if(containedInMiniHeap(start, 2*off, bset))
				return true;
			if(containedInMiniHeap(start, 2*off+1, bset))
				return true;
		}
		return false;
	}

public:

	enum ConstructAlg {GreedyBottomUp, PeelTopDown, PCASplit};

	BitSetTree(): root(0), levels(0), heapSkip(0)
	{

	}

	void construct(const vector<TripletFingerprint>& fingers, ConstructAlg alg)
	{
		tree.clear();
		root = 0;
		levels = 0;
		if(alg == GreedyBottomUp)
			constructGreedyBU(fingers);
		else if(alg == PeelTopDown)
			construectPeelTD(fingers);
		else if(alg == PCASplit)
			constructPCA(fingers);

		computeDists(root);
		addNonZeroRoots(root);
		if(nonzeroRoots.size() > 0)
		{
			unsigned d = tree[nonzeroRoots[0]].dist;
			heapSkip =  1<<(d+1);
			createFingerHeap_r(root,d);
		}
		//printInfo();
	}

	BitSetTree(const vector<TripletFingerprint>& fingers, ConstructAlg alg): root(0), levels(0)
	{
		construct(fingers, alg);
		printInfo();
	}
	~BitSetTree() {}

	void printInfo()
	{
		if(root == 0) return;
		//first count the number of identical nodes that aren't leaves
		unsigned identicalCnt = 0;
		unsigned leafCnt = 0;
		for(unsigned i = 1, n = tree.size(); i < n; i++)
		{
			if(tree[i].left == 0 && tree[i].right == 0)
				leafCnt++;
			if(tree[i].left != 0 && tree[i].anded == tree[i].ored)
			{
				identicalCnt++;
			}

		}

		cout << "Leaves: " << leafCnt << "\n";
		cout << "Identical: " << identicalCnt << "\n";

		vector<vector<int> > levelbitcnts;

		computeBitCnts(levelbitcnts, tree[root], 0);

		BOOST_FOREACH(vector<int>& vec, levelbitcnts)
		{
			BOOST_FOREACH(int b, vec)
			{
				cout << b << " ";
			}
			cout << "\n";
		}
	}

	struct CheckCnt
	{
		unsigned cnt;
		CheckCnt(): cnt(0) {}
		void incr() {cnt++;}
		~CheckCnt() { cout << "CHCKCNT " << cnt << "\n";}
	};
	//return true if there exists a finger stored in the bitset tree that is contained in bset
	bool hasContainedIn(const TripletFingerprint& bset) const
	{
		for(unsigned i = 0, n = fingerHeap.size(); i < n; i += heapSkip)
		{
			if(containedInMiniHeap(i, 1, bset))
				return true;
			/* avoiding recursion.. no any faster
			bool keepgoing = true;
			for(unsigned off = 1, half = heapSkip/2; off < half; off++)
			{
				if((off & (off-1)) == 0) //powr of 2
				{
					if (!keepgoing)
						break;
					keepgoing = false;
				}
				if (bset.contains(fingerHeap[i + off]))
					keepgoing = true;
			}
			if (keepgoing)
			{
				for (unsigned off = heapSkip / 2, stop = heapSkip; off < stop; off++)
				{
					//leaves, an positive stops
					if (bset.contains(fingerHeap[i + off]))
						return true;
				}
			}
			*/
		}
		return false;
	}
};

#endif /* PHARMITSERVER_BITSETTREE_H_ */
