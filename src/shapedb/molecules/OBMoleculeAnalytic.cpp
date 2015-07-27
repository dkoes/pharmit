/*
 * OBMoleculeAnalytic.cpp
 *
 *  Created on: Mar 14, 2012
 *      Author: dkoes
 */

#include "OBMoleculeAnalytic.h"
#include <boost/multi_array.hpp>
#include "CommandLine2/CommandLine.h"

extern cl::opt<bool> Verbose;
using namespace boost;
//maintain a bounding box
struct BoundingBox
{
	double minx;
	double maxx;
	double miny;
	double maxy;
	double minz;
	double maxz;

	BoundingBox() :
			minx(HUGE_VAL), maxx(-HUGE_VAL), miny(HUGE_VAL), maxy(-HUGE_VAL), minz(
					HUGE_VAL), maxz(-HUGE_VAL)
	{
	}

	void update(double x, double y, double z)
	{
		if (x < minx)
			minx = x;
		if (x > maxx)
			maxx = x;
		if (y < miny)
			miny = y;
		if (y > maxy)
			maxy = y;
		if (z < minz)
			minz = z;
		if (z > maxz)
			maxz = z;
	}

	//add amount to all dimensions
	void extend(double amount)
	{
		minx -= amount;
		maxx += amount;
		miny -= amount;
		maxy += amount;
		minz -= amount;
		maxz += amount;
	}

	//x,y,z within box
	bool contains(double x, double y, double z) const
	{
		if (x >= minx && x <= maxx && y >= miny && y <= maxy && z >= minz
				&& z <= maxz)
		{
			return true;
		}
		return false;
	}
};

struct AtomPoint
{
	double x;
	double y;
	double z;
	unsigned interactingCnt;

	AtomPoint() :
			x(0), y(0), z(0), interactingCnt(0)
	{
	}
	AtomPoint(double _x, double _y, double _z) :
			x(_x), y(_y), z(_z), interactingCnt(0)
	{
	}

	double distSq(double X, double Y, double Z) const
	{
		return (x - X) * (x - X) + (y - Y) * (y - Y) + (z - Z) * (z - Z);
	}

};

typedef pair<unsigned, vector<unsigned> > ClusterAdjList;

static bool neighSorter(const ClusterAdjList& lhs, const ClusterAdjList& rhs)
{
	return lhs.second.size() < rhs.second.size();
}

//checks that adding c to cluster won't violate maxClusterDist
static bool goodForCluster(const vector<unsigned>& cluster, unsigned c, const multi_array<double, 2>& distances, double maxdistSq)
{
	for(unsigned i = 0, n = cluster.size(); i < n; i++)
	{
		unsigned index = cluster[i];
		if(distances[index][c] > maxdistSq)
		{
			return false;
		}
	}
	return true;
}

//create clusters of points; clusters contains indices into points
static void clusterPoints(const vector<AtomPoint>& points,
		double maxClusterDist, vector<vector<unsigned> >& clusters)
{
	unsigned N = points.size();
	multi_array<double, 2> distances(extents[N][N]);

	//compute distances array and count close neighbors
	vector<ClusterAdjList> closeNeighbors(N);
	double maxdistSq = maxClusterDist * maxClusterDist;

	for (unsigned i = 0; i < N; i++)
	{
		closeNeighbors[i].first = i;
		for (unsigned j = 0; j < N; j++)
		{
			distances[i][j] = points[i].distSq(points[j].x, points[j].y,
					points[j].z);

			if (i != j && distances[i][j] <= maxdistSq)
			{
				closeNeighbors[i].second.push_back(j);
			}
		}
	}

	//sort point indices by how many neighbors are within the cutoff
	sort(closeNeighbors.begin(), closeNeighbors.end(), neighSorter);

	clusters.clear();
	vector<int> clusterAssign(N, -1);
	//greedily assign to clusters
	for (unsigned i = 0; i < N; i++)
	{
		unsigned index = closeNeighbors[i].first;
		if (clusterAssign[index] < 0) //unassigned, create new cluster
		{
			unsigned cluster = clusters.size();
			clusters.push_back(vector<unsigned>());
			clusters.back().push_back(index);

			//now consider all close neighbors for inclusion in cluster
			for (unsigned j = 0, m = closeNeighbors[i].second.size(); j < m;
					j++)
			{
				unsigned neigh = closeNeighbors[i].second[j];
				if(clusterAssign[neigh] < 0 && goodForCluster(clusters.back(), neigh, distances, maxdistSq))
				{
					clusters.back().push_back(neigh);
					clusterAssign[neigh] = cluster;
				}
			}

		}
	}
}

//computes a set of solitary grid points that represent the interaction between
//this ligand and the provided receptor in some way
void OBAMolecule::computeInteractionGridPoints(OBAMolecule& receptor,
		MGrid& grid, double interactionDist,
		double maxClusterDist,
		unsigned minClusterPoints, double interactionPointRadius)
{
	grid.clear();
	//first construct a bounding box for the ligand while assembling a
	//vector of atomic coordinates
	BoundingBox ligandBox;
	vector<AtomPoint> points;
	points.reserve(mol.NumAtoms());
	for (OBAtomIterator aitr = mol.BeginAtoms(); aitr != mol.EndAtoms(); ++aitr)
	{
		OBAtom* atom = *aitr;
		points.push_back(AtomPoint(atom->x(), atom->y(), atom->z()));
		ligandBox.update(atom->x(), atom->y(), atom->z());
	}
	ligandBox.extend(interactionDist);

	//then identify all coordinates that are interacting
	double idistSq = interactionDist * interactionDist;
	OBMol& rmol = receptor.getMol();
	for (OBAtomIterator aitr = rmol.BeginAtoms(); aitr != rmol.EndAtoms();
			++aitr)
	{
		OBAtom* a = *aitr;
		if (ligandBox.contains(a->x(), a->y(), a->z()))
		{
			for (unsigned i = 0, n = points.size(); i < n; i++)
			{
				if (points[i].distSq(a->x(), a->y(), a->z()) <= idistSq)
				{
					points[i].interactingCnt++;
				}
			}
		}
	}

	//prune out non-interacting poitns
	vector<AtomPoint> tmp;
	tmp.reserve(points.size());
	for (unsigned i = 0, n = points.size(); i < n; i++)
	{
		if (points[i].interactingCnt > 0)
			tmp.push_back(points[i]);
	}
	points.swap(tmp);
	tmp.clear();

	//cluster these coordinates
	vector<vector<unsigned> > clusters;
	clusterPoints(points, maxClusterDist, clusters);

	//make the cluster centers the interaction grid points
	for (unsigned i = 0, n = clusters.size(); i < n; i++)
	{
		double xtot = 0, ytot = 0, ztot = 0;
		unsigned npts = clusters[i].size();

		if (npts >= minClusterPoints)
		{
			for (unsigned j = 0; j < npts; j++)
			{
				xtot += points[clusters[i][j]].x;
				ytot += points[clusters[i][j]].y;
				ztot += points[clusters[i][j]].z;
			}
			double xave = xtot / (double) npts;
			double yave = ytot / (double) npts;
			double zave = ztot / (double) npts;

			grid.setPoint(xave, yave, zave);
			if(interactionPointRadius > 0)
			{
				grid.markXYZSphere(xave,yave,zave,interactionPointRadius);
			}
			if(Verbose)
			{
				cout << "pseudoatom InteractionPoints, pos=(" << xave << "," << yave << "," << zave << ")\n";
			}
		}
	}
}
