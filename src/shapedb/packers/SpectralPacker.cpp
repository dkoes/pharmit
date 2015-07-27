/*
 * SpectralPacker.cpp
 *
 *  Created on: Oct 24, 2011
 *      Author: dkoes
 */

#include "SpectralPacker.h"
#include "ShapeDistance.h"
#include <boost/shared_ptr.hpp>

#include "CommandLine2/CommandLine.h"

cl::opt<unsigned> KSlice("kslice",cl::Hidden,cl::init(0));

bool SpectralPacker::EigenInd::operator<(const EigenInd& rhs) const
{
	if(eigens == NULL)
		return false;
	if(rhs.eigens == NULL)
		return true;
	const double threshold = 0.00001;
	for (unsigned i = 0, n = eigens->cols(); i < n; i++)
	{
		double lval = (*eigens)(index_of_index, i);
		double rval = (*rhs.eigens)(rhs.index_of_index, i);

		if (lval < rval - threshold)
			return true;
		if (lval > rval + threshold)
			return false;
	}
	return false;
}

//compute distances to next
void SpectralPacker::EigenInd::computeNextDistance(const EigenInd& next)
{
	VectorXd diff = eigens->row(index_of_index)
			- next.eigens->row(next.index_of_index);
	distance_to_next = diff.norm();
}

//create a similarity matrix from graph
//compute log(n) nearest neighbors to determine a good sigma
//or alternatively compute an actual knn graph (todo - a little tricky
//since need to keep the graph connected)
void SpectralPacker::transformDistancesToSimilarity(
		const multi_array<double, 2>& distances, MatrixXd& graph) const
{
	unsigned N = distances.size(); //must be square
	graph = MatrixXd::Zero(N, N);
	unsigned k = max(ceil(log2(N)), 1.0);

	vector<double> kthdist(N, 0);
	double sum = 0;
	double maxdist = 0;
	//compute the k nearest
	vector<double> row(N);
	for (unsigned i = 0; i < N; i++)
	{
		//get distances
		for (unsigned j = 0; j < N; j++)
		{
			if (j == i)
				row[j] = HUGE_VAL;
			else
				row[j] = distances[i][j];
		}
		sort(row.begin(), row.end());
		kthdist[i] = row[k - 1];
		if(row[N-2] > maxdist)
			maxdist = row[N-2];
		sum += row[k - 1];
	}
	double sigma = sum / N / 4; //arbitrary
	double denom = 2 * sigma * sigma;
	//now compute Gaussian similarity
	for (unsigned i = 0; i < N; i++)
	{
		for (unsigned j = 0; j < N; j++)
		{
			if (i == j)
				graph(i, j) = 0;
			else
			{
				//graph(i,j) = maxdist - distances[i][j];
				double d = distances[i][j];
				d *= d; //sqr
				graph(i, j) = exp(-d / denom);
			}
		}
	}
}

//The unnormalized Laplacian is D - W
void SpectralPacker::transformSimilarityToLaplacian(MatrixXd& graph,
		MatrixXd& degrees) const
{
	unsigned N = graph.rows();
	VectorXd ones = VectorXd::Ones(N);
	VectorXd dvec = graph * ones;
	degrees = dvec.asDiagonal();

	graph *= -1;
	graph += degrees;
}

//append a single cluster
void SpectralPacker::createCluster(const DataViewer* dv,
		const vector<EigenInd>& vals, unsigned start, unsigned end,
		vector<Cluster>& clusters) const
{
	if (start >= end)
		return;
	unsigned N = end - start;

	const MappableOctTree *mivs[N];
	const MappableOctTree *msvs[N];

	vector<unsigned> indices(N);

	for (unsigned i = 0; i < N; i++)
	{
		unsigned index = vals[start + i].index;
		indices[i] = index;
		mivs[i] = dv->getMIV(index);
		msvs[i] = dv->getMSV(index);
	}

	clusters.push_back(Cluster(indices, mivs, msvs));
}

//append a single cluster
void SpectralPacker::createCluster(const DataViewer* dv,
		const vector<unsigned>& vals /* index into indices */, const vector<unsigned>& indices,
		vector<Cluster>& clusters) const
{
	unsigned N = vals.size();

	const MappableOctTree *mivs[N];
	const MappableOctTree *msvs[N];

	vector<unsigned> cindices(N);

	for (unsigned i = 0; i < N; i++)
	{
		unsigned index = indices[vals[i]];
		cindices[i] = index;
		mivs[i] = dv->getMIV(index);
		msvs[i] = dv->getMSV(index);
	}

	clusters.push_back(Cluster(cindices, mivs, msvs));
}

//given the ordering in vals between start and end, create packed clusters
//The goal is to have all clusters as dense as possible
void SpectralPacker::createDenseClusters(const DataViewer* dv,
		const vector<EigenInd>& vals, unsigned start, unsigned end,
		SolverPtr solver, vector<Cluster>& clusters) const
{
	if (start >= end)
		return;
	unsigned sz = end - start;

	if (sz % packSize == 0) //evenly divided, go for it
	{
		for (unsigned i = start; i < end; i += packSize)
		{
			createCluster(dv, vals, i, i + packSize, clusters);
		}
	}
	else if (sz < packSize)
	{
		createCluster(dv, vals, start, end, clusters);
	}
	else //peel off the largest chunk we can to get balanced clusters
	{
		unsigned two = (sz % packSize) + packSize;
		unsigned first = two / 2;
		unsigned second = two - first;
		createCluster(dv, vals, start, start + first, clusters);
		createCluster(dv, vals, start + first, start + first + second,
				clusters);
		createDenseClusters(dv, vals, start + first + second, end, solver,
				clusters);
	}
}

//split along the largest difference
void SpectralPacker::createPartitionedClusters(const DataViewer* dv,
		const vector<EigenInd>& vals, unsigned start, unsigned end,
		SolverPtr solver, vector<Cluster>& clusters) const
{
	unsigned sz = end - start;
	if (sz < 4 * packSize)
	{
		createDenseClusters(dv, vals, start, end, solver, clusters);
	}
	else
	{
		double max = 0;
		unsigned pos = start + sz / 2;
		//find largest difference
		for (unsigned i = start + packSize; i < end - packSize; i++)
		{
			if (vals[i].distance_to_next > max)
			{
				max = vals[i].distance_to_next;
				pos = i + 1;
			}
		}

		createPartitionedClusters(dv, vals, start, pos, solver, clusters);
		createPartitionedClusters(dv, vals, pos, end, solver, clusters);
	}
}


static double completeLink(const MatrixXd& D, const vector<unsigned>& a, const vector<unsigned>& b)
{
	//this is the maximum of the minimum distances between cluster members
	double max = 0;
	for (unsigned i = 0, ni = a.size(); i < ni; i++)
	{
		float min = HUGE_VAL;
		for (unsigned j = 0, nj = b.size(); j < nj; j++)
		{
			double dist = D(a[i],b[j]);

			if (dist < min)
				min = dist;
		}
		if (min > max)
			max = min;
	}
	return max;
}

//combine the rows with the closest distances in D, return false if no merging possible
bool SpectralPacker::mergeRowIndices(const MatrixXd& D, vector<vector<unsigned> >& clusters) const
{
	unsigned N = clusters.size();
	if (N == 1)
		return false;

	vector<IntraClusterDist> distances;
	distances.reserve(N * N / 2);

	for (unsigned i = 0; i < N; i++)
	{
		for (unsigned j = 0; j < i; j++)
		{
			double dist = completeLink(D, clusters[i], clusters[j]);
			distances.push_back(IntraClusterDist(i, j, dist));
		}
	}
	sort(distances.begin(), distances.end());

	vector<vector<unsigned> > newclusters;
	newclusters.reserve(N / 2 + 1);
	unsigned merged = 0;
	bool didmerge = false;
	vector<bool> alreadyMerged(N, false);
	for (unsigned d = 0, maxdists = distances.size(); d < maxdists; d++)
	{
		//merge if not already merged
		unsigned i = distances[d].i;
		unsigned j = distances[d].j;
		if (!alreadyMerged[i] && !alreadyMerged[j]
				&& clusters[i].size() + clusters[j].size() <= packSize)
		{
			newclusters.push_back(clusters[i]);
			newclusters.back().insert(newclusters.back().end(), clusters[j].begin(), clusters[j].end());

			alreadyMerged[i] = true;
			alreadyMerged[j] = true;
			merged += 2;
			didmerge = true;
		}

		if (merged >= N - 1)
			break;
	}

	if (merged != N) //find the loner(s), keep it as a singleton (may also be too big)
	{
		for (unsigned i = 0; i < N; i++)
		{
			if (!alreadyMerged[i])
			{
				newclusters.push_back(clusters[i]);
			}
		}
	}

	swap(newclusters, clusters);

	return didmerge;
}

struct UDPair
{
	unsigned val;
	double max;

	UDPair(unsigned r, double m): val(r), max(m) {}
	bool operator<(const UDPair& rhs) const //put largest first
	{
		return max > rhs.max;
	}
};
//greedily construct clusters from the relaxation, relax(i,j) is node i's preference for cluster j
void SpectralPacker::createClustersFromRelaxation(const DataViewer *dv, const MatrixXd& relax,  const vector<unsigned>& indices, vector<Cluster>& clusters) const
{
	vector<bool> packed(relax.rows(), false);
	//sort clusters based on maximum preference
	MatrixXd maxes = relax.colwise().maxCoeff();

	vector<UDPair> Cs;
	for(unsigned i = 0, n = maxes.cols(); i < n; i++)
	{
		Cs.push_back(UDPair(i, maxes(i)));
	}
	sort(Cs.begin(), Cs.end());

	//choose the packSize top indices that prefer this cluster and aren't already taken
	for(unsigned i = 0, n = Cs.size(); i < n; i++)
	{
		vector<UDPair> Is;
		unsigned c = Cs[i].val;
		for(unsigned j = 0, nr = relax.rows(); j < nr; j++)
		{
			if(!packed[j]) {
				Is.push_back(UDPair(j,(double)relax(j, c)));
			}
		}
		sort(Is.begin(), Is.end());
		vector<unsigned> clust;
		for(unsigned j = 0, m = Is.size(); j < packSize && j < m; j++)
		{
			clust.push_back(Is[j].val);
			packed[Is[j].val] = true;
		}
		createCluster(dv, clust, indices, clusters);
	}
}


void SpectralPacker::pack(const DataViewer* dv,
		vector<Cluster>& clusters) const
{
	if (dv->size() == 0)
		return;

	vector<unsigned> indices(dv->size());
	for(unsigned i = 0, n = indices.size(); i < n; i++)
		indices[i] = i;

	//first compute the n^2 distances
	unsigned N = indices.size();
	multi_array<double, 2> distances(extents[N][N]);
	for (unsigned i = 0; i < N; i++)
	{
		distances[i][i] = 0;
		const MappableOctTree *imiv = dv->getMIV(indices[i]);
		const MappableOctTree *imsv = dv->getMSV(indices[i]);

		for (unsigned j = 0; j < i; j++)
		{
			const MappableOctTree *jmiv = dv->getMIV(indices[j]);
			const MappableOctTree *jmsv = dv->getMSV(indices[j]);
			distances[i][j] = distances[j][i] = shapeDistance(imiv, imsv, jmiv,
					jmsv);
		}
	}

	//generate similarity graph
	MatrixXd graph(N, N);
	transformDistancesToSimilarity(distances, graph);

	MatrixXd sim = graph;
	//compute unnormalized graph laplacian
	MatrixXd D(N, N);
	transformSimilarityToLaplacian(graph, D);

	switch (algo)
	{
	case SortDense:
	case SortPartition:
	{
		//find the second eigenvectors (first should be unit)
		SolverPtr solver;
		if (useNormalizedLaplacian)
			solver = SolverPtr(
					new GeneralizedSelfAdjointEigenSolver<MatrixXd>(graph, D));
		else
			solver = SolverPtr(new SelfAdjointEigenSolver<MatrixXd>(graph));

		//associate coordinates of eigenvector with indices and sort
		vector<EigenInd> vals(N);
		for (unsigned i = 0; i < N; i++)
		{
			vals[i] = EigenInd(i, indices[i], solver->eigenvectors());
		}

		//pack as densely as possible
		sort(vals.begin(), vals.end());

		if (algo == SortDense)
		{
			createDenseClusters(dv, vals, 0, vals.size(), solver, clusters);
		}
		else
		{
			for (unsigned i = 0; i < N - 1; i++)
			{
				vals[i].computeNextDistance(vals[i + 1]);
			}
			createPartitionedClusters(dv, vals, 0, vals.size(), solver,
					clusters);
		}
	}
		break;

	case ClusterFullEigen:
	{
		//find the second eigenvectors (first should be unit)
		SolverPtr solver;
		if (useNormalizedLaplacian)
			solver = SolverPtr(
					new GeneralizedSelfAdjointEigenSolver<MatrixXd>(graph, D));
		else
			solver = SolverPtr(new SelfAdjointEigenSolver<MatrixXd>(graph));

		//compute distances between eigen vectors
		MatrixXd D = MatrixXd::Zero(N,N);

		MatrixXd vecs = solver->eigenvectors();
		for(unsigned i = 0; i < N; i++)
		{
			for(unsigned j = 0; j < i; j++)
			{
				MatrixXd d = vecs.row(i) - vecs.row(j);
				if(KSlice > 0)
				{
					d.conservativeResize(1, KSlice);
				}

				double dist = d.norm();
				D(i,j) = dist;
				D(j,i) = dist;
			}
		}

		//cluster rows based on distances
		vector< vector<unsigned> > iclusters(N);
		for(unsigned i = 0; i < N; i++)
		{
			iclusters[i].push_back(i);
		}

		while(mergeRowIndices(D,iclusters))
			;

		for(unsigned i = 0, n = iclusters.size(); i < n; i++)
		{
			createCluster(dv, iclusters[i], indices, clusters);
		}

	}
		break;
	case RelaxationPacking:
		{
			double y = -1/sqrt(N);
			double x = -1.0/(N+sqrt(N));
			MatrixXd V(N,N-1);
			for(unsigned i = 0; i < N-1; i++)
			{
				V(0,i) = y;
			}

			for(unsigned i = 1; i < N; i++)
			{
				for(unsigned j = 0; j < N-1; j++)
				{
					V(i,j) = x;
					if(i == j+1)
						V(i,j)++;
				}
			}

			MatrixXd g2 = V.transpose()*sim*V;
			SelfAdjointEigenSolver<MatrixXd> bar(g2);

			unsigned k = ceil(N/(double)packSize);
			MatrixXd W(k,k-1);
			y = -1/sqrt(k);
			x = -1.0/(k+sqrt(k));
			for(unsigned i = 0; i < k-1; i++)
			{
				W(0,i) = y;
			}

			for(unsigned i = 1; i < k; i++)
			{
				for(unsigned j = 0; j < k-1; j++)
				{
					W(i,j) = x;
					if(i == j+1)
						W(i,j)++;
				}
			}


			MatrixXd Z = bar.eigenvectors().rightCols(k-1).rowwise().reverse();
			Z.conservativeResize(N-1,k-1);

			MatrixXd X = V*Z*W.transpose();
			X = X.array()*sqrt(N/(double)k) + 1.0/(double)k;
			createClustersFromRelaxation(dv, X, indices, clusters);

		}
		break;
	}

}
