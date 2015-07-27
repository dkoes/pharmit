/*
 * SpectralPacker.h
 *
 *  Created on: Oct 24, 2011
 *      Author: dkoes
 *
 *  Packs bottom up using the spectral ordering - the second eigenvector
 *  of a graph Laplacian. Enhances partitioning using additional eigenvectors.
 */

#ifndef SPECTRALPACKER_H_
#define SPECTRALPACKER_H_

#include "Packer.h"
#include "boost/multi_array.hpp"

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

using namespace boost;
using namespace Eigen;


class SpectralPacker: public Packer
{
public:
	enum SpectralAlgEnum {SortDense, SortPartition, ClusterFullEigen,RelaxationPacking};
private:
	//small container for storing and sorting values from eigen vector and index
	struct EigenInd
	{
		const MatrixXd* eigens;
		unsigned index; //value from indices
		unsigned index_of_index; //position within indices, for referencing whole eigenvalue
		double distance_to_next;
		EigenInd(): eigens(NULL), index(0), index_of_index(0), distance_to_next(0) {}
		EigenInd(unsigned i, unsigned i2, const MatrixXd& e): eigens(&e), index(i2), index_of_index(i), distance_to_next(0) {}

		bool operator<(const EigenInd& rhs) const;
		void computeNextDistance(const EigenInd& next);
	};
	typedef shared_ptr< SelfAdjointEigenSolver<MatrixXd> > SolverPtr;

	bool useNormalizedLaplacian;
	SpectralAlgEnum algo;


	void transformDistancesToSimilarity(const multi_array<double, 2>& distances, MatrixXd& graph) const;
	void transformSimilarityToLaplacian(MatrixXd& graph, MatrixXd& degrees) const;

	void createCluster(const DataViewer* dv, const vector<EigenInd>& vals,
			unsigned start, unsigned end, vector<Cluster>& clusters) const;
	void createCluster(const DataViewer* dv,
			const vector<unsigned>& vals /* index into indices */, const vector<unsigned>& indices,
			vector<Cluster>& clusters) const;

	void createClustersFromRelaxation(const DataViewer *dv, const MatrixXd& relax,  const vector<unsigned>& indices, vector<Cluster>& clusters) const;
	void createDenseClusters(const DataViewer* dv, const vector<EigenInd>& vals,
			unsigned start, unsigned end, SolverPtr solver, vector<Cluster>& clusters) const;

	void createPartitionedClusters(const DataViewer* dv, const vector<EigenInd>& vals,
			unsigned start, unsigned end, SolverPtr solver, vector<Cluster>& clusters) const;

	bool mergeRowIndices(const MatrixXd& D, vector<vector<unsigned> >& clusters) const;

public:
	SpectralPacker(unsigned pk, SpectralAlgEnum alg=SortDense, bool norm=true): Packer(pk, NotApplicable), useNormalizedLaplacian(norm), algo(alg) {}
	~SpectralPacker() {}

	void pack(const DataViewer* dv, vector<Cluster>& clusters) const;

};

#endif /* SPECTRALPACKER_H_ */
