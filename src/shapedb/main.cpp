/*
 ShapeDB
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
 * main.cpp
 *
 *  Created on: Sep 12, 2011
 *      Author: dkoes

 The goal of this project is to efficiently store molecular shapes within
 a GSS-tree index for fast (and exact) similarity searching
 http://dx.doi.org/10.1145/304181.304219

 At some point I intend to properly break this out into a library for generic
 shape handling, but at the moment everything assumes we are performing
 molecular shape matching.
 */

#include <iostream>
#include <string>
#include "GSSTypes.h"
#include "boost/array.hpp"

#include "CommandLine2/CommandLine.h"
#include "GSSTreeCreator.h"
#include "GSSTreeSearcher.h"
#include "molecules/Molecule.h"
#include "KSamplePartitioner.h"
#include "packers/Packers.h"
#include <boost/shared_ptr.hpp>
#include "Timer.h"
#include "MiraObject.h"

using namespace std;
using namespace boost;

typedef shared_ptr<Packer> PackerPtr;

enum CommandEnum
{
	Create,
	MiraSearch,
	DCMiraSearch,
	NNSearch,
	DCSearch,
	MolGrid,
	MergeMira,
	BatchSearch,
	BatchMiraSearch,
	BatchDB,
	SearchAllPointCombos,
	CreateTrees,
	CreateFromTrees,
	CreateFromMira
};

cl::opt<CommandEnum> Command(cl::desc("Operation to perform:"), cl::Required,
		cl::values(clEnumVal(Create, "Create a molecule shape index."),
				clEnumVal(CreateTrees,"Create only the shapes with no index."),
				clEnumVal(CreateFromTrees,"Create an index from already created shapes."),
				clEnumVal(CreateFromMira,"Create an index from already created shapes."),
				clEnumVal(MiraSearch, "Nearest neighbor searching of mira objects"),
				clEnumVal(DCMiraSearch, "Distance constraint searching of mira objects"),
				clEnumVal(NNSearch, "Nearest neighbor search"),
				clEnumVal(DCSearch, "Distance constraint search"),
				clEnumVal(MolGrid, "Generate molecule grid and debug output"),
				clEnumVal(MergeMira, "Merge mira files into MIV/MSV"),
				clEnumVal(BatchSearch, "Read in a jobs file for batch processing"),
				clEnumVal(BatchMiraSearch, "Read in a mira jobs file for batch processing"),
				clEnumVal(BatchDB, "Read in a jobs file for batch processing of a list of directories"),
				clEnumVal(SearchAllPointCombos, "Search using all possible subsets of interaction points"),
				clEnumValEnd));

cl::opt<bool> NNSearchAll("nn-search-all",
		cl::desc(
				"Perform exhaustive interaction point search with NNSearch as well as DCSearch"),
		cl::init(false));

enum PackerEnum
{
	FullMerge, Spectral, GreedyMerge, MatchPack
};
cl::opt<PackerEnum> PackerChoice(cl::desc("Packing algorithm:"),
		cl::values(clEnumValN(FullMerge,"full-merge", "Greedy full merge."),
				clEnumValN(GreedyMerge,"greedy-merge", "Greedy iterative merge."),
				clEnumValN(MatchPack,"match-merge", "Optimal matching merging."),
				clEnumValN(Spectral, "spectral", "Spectral packing"),
				clEnumValEnd), cl::init(MatchPack));

cl::opt<unsigned> K("k", cl::desc("k nearest neighbors to find for NNSearch"),
		cl::init(1));
cl::opt<double> NNThreshold("nn-threshold",cl::desc("Similarity cutoff for NNSearch"), cl::init(HUGE_VAL));

cl::opt<unsigned> Knn("knn", cl::desc("K for knn graph creation"), cl::init(8));
cl::opt<unsigned> Sentinals("sentinals",
		cl::desc("Number of sentinals for knn initialization (zero random)"),
		cl::init(32));
cl::opt<SpectralPacker::SpectralAlgEnum> SpectralAlg(
		cl::desc("Spectral packing sub-algorithm:"),
		cl::values(
				clEnumValN(SpectralPacker::SortDense, "sort-dense", "Simple sort followed by dense packing"),
				clEnumValN(SpectralPacker::SortPartition, "sort-partition", "Simple sort followed by largest separator partition packing"),
				clEnumValN(SpectralPacker::ClusterFullEigen, "cluster-eigen", "Cluster eigen values using greedy packer"),
				clEnumValN(SpectralPacker::RelaxationPacking, "relax", "Cluster form relaxation values"),
				clEnumValEnd), cl::init(SpectralPacker::SortDense));

cl::opt<bool> ScanCheck("scancheck",
		cl::desc("Perform a full scan to check results"), cl::Hidden);
cl::opt<bool> ScanOnly("scanonly", cl::desc("Search using only a scan"),
		cl::Hidden);

cl::opt<bool> SingleConformer("single-conformer",
		cl::desc("Output the single best conformer"), cl::init(false));

cl::opt<string> Input("in", cl::desc("Input file"));
cl::opt<string> Output("out", cl::desc("Output file"));
cl::opt<string> Database("db", cl::desc("Database file"));

cl::list<string> Files("files",cl::desc("Files for MiraMerge"));

cl::opt<double> LessDist("less",
		cl::desc(
				"Distance to reduce query mol by for constraint search (default 1A)."),
		cl::init(1.0));
cl::opt<double> MoreDist("more",
		cl::desc(
				"Distance to increase query mol by for constraint search (default 1A)."),
		cl::init(1.0));

cl::opt<double> ProbeRadius("probe-radius",
		cl::desc("Radius of water probe for exmol only"), cl::init(1.4));

cl::opt<double> MaxDimension("max-dim", cl::desc("Maximum dimension."),
		cl::init(64));
cl::opt<double> Resolution("resolution",
		cl::desc("Best resolution for shape database creation."), cl::init(.5));

cl::opt<string> IncludeMol("ligand",
		cl::desc("Molecule to use for minimum included volume"));
cl::opt<string> ExcludeMol("receptor",
		cl::desc("Molecule to use for excluded volume"));

cl::opt<string> MIVShape("mira-miv", cl::desc("Mira shape to use for minimum included volume"));
cl::opt<string> MSVShape("mira-msv", cl::desc("Mira shape to use for maximum surrounding volume"));

cl::opt<bool> useInteractionPoints("use-interaction-points",
		cl::desc(
				"Analyze the ligand-receptor complex and generate a minimum shape centered around different interaction points."),
		cl::init(false));
cl::opt<double> interactionPointRadius("interaction-point-radius",
		cl::desc("Amount to grow interaction points into ligand shape."),
		cl::init(0));
cl::opt<double> interactionDistance("interaction-distance",
		cl::desc(
				"Distance between ligand/receptor atom centers that are considered interacting"),
		cl::init(6));
cl::opt<double> interactionMinCluster("interaction-min-cluster",
		cl::desc("Minimum size of interaction point cluster."), cl::init(3));
cl::opt<double> interactionMaxClusterDist("interaction-max-cluster-dist",
		cl::desc("Maximum span of interaction point cluster."), cl::init(4));

cl::opt<unsigned> KCenters("kcenters",
		cl::desc("number of centers for ksample-split"), cl::init(8));
cl::opt<unsigned> KSampleMult("ksamplex",
		cl::desc("multiplictive factor for ksampling"), cl::init(5));

cl::opt<unsigned> SwitchToPack("switch-to-pack",
		cl::desc("Cutoff to trigger packing in nodes and leaves"),
		cl::init(32768));

cl::opt<unsigned> Pack("pack", cl::desc("Maximum quantities per a node"),
		cl::init(16));

cl::opt<bool> Print("print",cl::desc("Print text summary of output"),cl::init(false));
cl::opt<bool> Verbose("v", cl::desc("Verbose output"));
cl::opt<bool> UseUnnorm("use-unnorm",
		cl::desc("Use unnormalized laplacian in spectral packing"),
		cl::init(false));

cl::opt<Packer::ClusterDistance> ClusterDist(
		cl::desc("Metric for cluster packing distance:"),
		cl::values(
				clEnumValN(Packer::AverageLink, "ave-dist", "Use 'average' metric between MIV/MSV representations of clusters"),
				clEnumValN(Packer::CompleteLink, "complete-dist", "Use complete linkage value between cluster members"),
				clEnumValN(Packer::SingleLink, "single-dist", "Use single linkage value between cluster members"),
				clEnumValN(Packer::TotalLink, "total-dist", "Use total (sum) linkage value between cluster members"),
				clEnumValEnd), cl::init(Packer::AverageLink));

cl::opt<DistanceFunction> ShapeDist(
		cl::desc("Metric for distance between shapes:"),
		cl::values(
				clEnumValN(RelativeVolume,"rel-volume","Relative volume difference"),
				clEnumValN(AbsVolume,"abs-volume", "Absolute volume difference"),
				clEnumValN(Hausdorff,"hausdorff", "Hausdorff distance"),
				clEnumValN(RelativeTriple, "rel-triple", "Triple including selectivity"),
				clEnumValN(AbsoluteTriple, "abs-triple", "Triple including selectivity (absolute)"),
				clEnumValN(IncludeExclude, "include-exclude","For comparing with include/exclude constraints"),
				clEnumValN(RelVolExclude, "relvol-exclude","Volume comparison of ligand, exclusion comparison of receptor"),
				clEnumValEnd), cl::init(RelativeVolume));

cl::opt<unsigned> SuperNodeDepth("superdepth",
		cl::desc("Depth to descend to create aggregrated super root"),
		cl::init(0));

cl::opt<unsigned> TimeTrials("time-trials",
		cl::desc("Number of runs to get average for benchmarking"),
		cl::init(1));
cl::opt<bool> ClearCache("clear-cache",
		cl::desc("Clear file cache between each benchmarking run"),
		cl::init(false));
cl::opt<bool> ClearCacheFirst("clear-cache-first",
		cl::desc("Clear file cache before each benchmarking run"),
		cl::init(false));

cl::opt<string> SproxelColor("sproxel-color",
		cl::desc("Sproxel voxel descriptor"), cl::init("#0000ffff"));

cl::opt<bool> KeepHydrogens("h",
		cl::desc("Retain hydrogens in input molecules"), cl::init(false));

typedef GSSTreeSearcher::ObjectTree ObjectTree;
//create min and max trees from molecular data
//caller takes owner ship of tree memory
static void create_trees(GSSTreeSearcher& gss, const string& includeMol,
		const string& excludeMol, double less, double more,
		GSSTreeSearcher::ObjectTree& smallTree,
		GSSTreeSearcher::ObjectTree& bigTree)
{
	double dimension = gss.getDimension();
	double resolution = gss.getResolution();

	//read query molecule(s)

	//use explicit volumes
	Molecule::iterator imolitr(includeMol, dimension, resolution, KeepHydrogens,
			ProbeRadius);
	Molecule inMol = *imolitr;

	Molecule exMol;

	if(filesystem::path(excludeMol).extension() != ".mira")
	{
		Molecule::iterator exmolitr(excludeMol, dimension, resolution,
				KeepHydrogens, ProbeRadius);
		exMol = *exmolitr;
		//create trees
		//receptor constraint, shrink by more and invert
		bigTree = gss.createTreeFromObject(exMol, more, excludeMol.size() > 0 ); //if exclude empty, make don't invert
	}
	else //is mira object, don't need to compute grid
	{
		MiraObject obj;
		ifstream mirain(excludeMol.c_str());
		obj.read(mirain);
		bigTree = gss.createTreeFromObject(obj, more, true);
	}

	//use the full shape of the ligand
	smallTree = gss.createTreeFromObject(inMol, less);

	if (useInteractionPoints)
	{
		//change small tree to just be interaction points
		MGrid igrid(dimension, resolution);
		inMol.computeInteractionGridPoints(exMol, igrid, interactionDistance,
				interactionMaxClusterDist, interactionMinCluster,
				interactionPointRadius);

		MGrid lgrid;
		smallTree->makeGrid(lgrid, resolution);
		lgrid &= igrid;
		smallTree = shared_ptr<const MappableOctTree>(
				MappableOctTree::createFromGrid(lgrid), free);
	}
}

//do search between include and exclude
static void do_dcsearch(GSSTreeSearcher& gss, ObjectTree smallTree,
		ObjectTree bigTree, const string& output)
{
	ResultMolecules res(SingleConformer);
	bool getres = output.size() > 0 || Print;
	//search
	if (!ScanOnly)
		gss.dc_search(smallTree, bigTree, ObjectTree(), getres, res);

	if (ScanCheck || ScanOnly)
	{
		ResultMolecules res2(SingleConformer);
		gss.dc_scan_search(smallTree, bigTree, ObjectTree(), getres, res2);
		if (res2.size() != res.size())
		{
			cerr << "Scanning found different number: " << res2.size() << "\n";
		}
	}

	res.writeOutput(output, Print);

}

void do_nnsearch(GSSTreeSearcher& gss, const string& input,
		const string& output, unsigned k, double thresh)
{
	//read query molecule(s)
	ResultMolecules res(SingleConformer);
	Molecule::iterator molitr(input, gss.getDimension(), gss.getResolution(),
			KeepHydrogens, ProbeRadius);
	bool getresults = output.size() > 0 || Print;
	for (; molitr; ++molitr)
	{
		const Molecule& mol = *molitr;
		ObjectTree objTree = gss.createTreeFromObject(mol);
		if (k < 1 && !isfinite(thresh)) //scan
		{
			gss.nn_scan(objTree, getresults, res);
		}
		else
		{
			gss.nn_search(objTree, k, thresh, getresults, res);
		}
	}

	res.writeOutput(output, Print);
}

struct QInfo
{
	string str;
	CommandEnum cmd;
	Molecule in;
	Molecule ex;
	unsigned k;
	double less;
	double more;

	QInfo() :
			cmd(Create), k(1), less(0), more(0)
	{
	}

	QInfo(const string& s, const Molecule& i, const Molecule& e, double l,
			double m) :
			str(s), cmd(DCSearch), in(i), ex(e), less(l), more(m)
	{
	}
	QInfo(const string& s, const Molecule& i, unsigned _k) :
			str(s), cmd(NNSearch), in(i), k(_k)
	{
	}
};

//attempt to clear the file system cache,
//this is just for benchmarking purposes and requires the existance of clearfilecache
//which is a setuid (gasp) program I wrote that does precisely that
static void clear_cache()
{
#pragma GCC diagnostic ignored "-Wunused-result"
	std::system("clearfilecache");
}

static void dumpTree(MappableOctTree *tree, ostream& out)
{

	cout << "Size of tree " << tree->bytes() << "\n";
	cout << "Volume of tree " << tree->volume() << "\n";
	cout << "Nodes of tree " << tree->nodes() << "\n";
	vector<unsigned> cnts;
	tree->countLeavesAtDepths(cnts);
	for (unsigned i = 0, n = cnts.size(); i < n; i++)
	{
		cout << i << " : " << cnts[i] << "\n";
	}
	if (out)
	{
		if (filesystem::extension(Output.c_str()) == ".raw")
			tree->dumpRawGrid(out, Resolution);
		if (filesystem::extension(Output.c_str()) == ".map")
			tree->dumpAD4Grid(out, Resolution);
		else if (filesystem::extension(Output.c_str()) == ".mira")
			tree->dumpMiraGrid(out, Resolution);
		else if (filesystem::extension(Output.c_str()) == ".csv")
			tree->dumpSproxelGrid(out, Resolution, SproxelColor);
		else
			tree->dumpGrid(out, Resolution);
	}
}

int main(int argc, char *argv[])
{
	cl::ParseCommandLineOptions(argc, argv);

	bool getres = Output.size() > 0 || Print;

	switch (Command)
	{
	case Create:
		case CreateFromTrees:
		case CreateFromMira:
		{
		//read in all the molecules and calculate the max bounding box
		KSamplePartitioner topdown(KCenters, KSampleMult,
				KSamplePartitioner::AveCenter, SwitchToPack);

		PackerPtr packer;
		switch (PackerChoice)
		{
		case FullMerge:
			packer = PackerPtr(
					new FullMergePacker(Pack, ClusterDist, Knn, Sentinals));
			break;
		case MatchPack:
			packer = PackerPtr(
					new MatcherPacker(Pack, Knn, Sentinals, ClusterDist));
			break;
		case GreedyMerge:
			packer = PackerPtr(
					new GreedyPacker(Pack, ClusterDist, Knn, Sentinals));
			break;
		case Spectral:
			packer = PackerPtr(
					new SpectralPacker(Pack, SpectralAlg, !UseUnnorm));
			break;
		}

		setDistance(ShapeDist);

		GSSLevelCreator leveler(&topdown, packer.get(), SwitchToPack,
				SwitchToPack);

		GSSTreeCreator creator(&leveler, SuperNodeDepth);

		filesystem::path dbpath(Database.c_str());

		if (Command == Create)
		{
			Molecule::iterator molitr(Input, MaxDimension, Resolution,
					KeepHydrogens, ProbeRadius);
			if (!creator.create<Molecule, Molecule::iterator>(dbpath, molitr,
					MaxDimension, Resolution))
			{
				cerr << "Error creating database\n";
				exit(1);
			}
		}
		else if (Command == CreateFromMira)
		{
			MiraObject::iterator miraitr(Input);
			if (!creator.create<MiraObject, MiraObject::iterator>(dbpath,
					miraitr, MaxDimension, Resolution))
			{
				cerr << "Error creating database\n";
				exit(1);
			}
		}
		else if (Command == CreateFromTrees) //trees already created
		{
			filesystem::path treedir(Input.c_str());
			if (!creator.create(dbpath, treedir, MaxDimension, Resolution))
			{
				cerr << "Error creating database\n";
				exit(1);
			}
		}

		if (Verbose)
			creator.printStats(cout);

	}
		break;
	case CreateTrees:
		{
		//read in all the molecules and calculate the max bounding box
		KSamplePartitioner topdown(KCenters, KSampleMult,
				KSamplePartitioner::AveCenter, SwitchToPack);

		PackerPtr packer = PackerPtr(
				new MatcherPacker(Pack, Knn, Sentinals, ClusterDist));

		setDistance(ShapeDist);

		GSSLevelCreator leveler(&topdown, packer.get(), SwitchToPack,
				SwitchToPack);

		GSSTreeCreator creator(&leveler, SuperNodeDepth);

		filesystem::path dbpath(Database.c_str());
		Molecule::iterator molitr(Input, MaxDimension, Resolution,
				KeepHydrogens, ProbeRadius);
		if (!creator.createTreesOnly<Molecule, Molecule::iterator>(dbpath,
				molitr, MaxDimension, Resolution))
		{
			cerr << "Error creating database\n";
			exit(1);
		}

	}
		break;
	case MiraSearch:
		{
		//read in database
		filesystem::path dbfile(Database.c_str());
		GSSTreeSearcher gss(Verbose);
		if (!gss.load(dbfile))
		{
			cerr << "Could not read database " << Database << "\n";
			exit(-1);
		}

		setDistance(ShapeDist);
		ifstream mirain(Input.c_str());
		MiraObject obj;
		obj.read(mirain);

		StringResults res;
		ObjectTree objTree = gss.createTreeFromObject(obj);

		if (K < 1) //scan
		{
			gss.nn_scan(objTree, true, res);
		}
		else
		{
			gss.nn_search(objTree, K, HUGE_VAL, true, res);
		}

		for (unsigned i = 0, n = res.size(); i < n; i++)
		{
			cout << i << "\t" << res.getString(i) << "\t" << res.getScore(i)
					<< "\n";
		}
	}
		break;
	case DCMiraSearch:
		{
		//read in database
		filesystem::path dbfile(Database.c_str());
		GSSTreeSearcher gss(Verbose);

		if (!gss.load(dbfile))
		{
			cerr << "Could not read database " << Database << "\n";
			exit(-1);
		}
		setDistance(ShapeDist);

		double resolution = gss.getResolution();

		setDistance(ShapeDist);

		StringResults res;
		ObjectTree objTree;
		ObjectTree smallTree, bigTree;

		if(Input.size() > 0)
		{
			ifstream mirain(Input.c_str());
			MiraObject obj;
			obj.read(mirain);
			objTree = gss.createTreeFromObject(obj);
		}

		//search with passed shapes, does not invert exclude
		if (MSVShape.size() > 0 && MIVShape.size() > 0)
		{
			MiraObject miv, msv;
			ifstream mivin(MIVShape.c_str());
			ifstream msvin(MSVShape.c_str());
			miv.read(mivin);
			msv.read(msvin);

			smallTree = gss.createTreeFromObject(miv);
			bigTree = gss.createTreeFromObject(msv);
		}
		else if(MSVShape.size() > 0 || MIVShape.size() > 0)
		{
			cerr << "Error, need to specify both MIV and MSV\n";
			exit(-1);
		}
		else //base off of query shape
		{
			MGrid smallgrid, biggrid;
			objTree->makeGrid(smallgrid, resolution);
			biggrid = smallgrid;
			smallgrid.shrink(LessDist);
			biggrid.grow(MoreDist);

			smallTree = ObjectTree(
					MappableOctTree::createFromGrid(smallgrid));
			bigTree = ObjectTree(
					MappableOctTree::createFromGrid(biggrid));
		}

		if (!ScanOnly)
			gss.dc_search(smallTree, bigTree, objTree, Output.size() > 1, res);

		if (ScanCheck || ScanOnly)
		{
			StringResults res2;
			gss.dc_scan_search(smallTree, bigTree, objTree, Output.size() > 1,
					res2);
			if (res2.size() != res.size())
			{
				cerr << "Scanning found different number\n";
			}
		}
		for (unsigned i = 0, n = res.size(); i < n; i++)
		{
			cout << i << "\t" << res.getString(i) << "\t" << res.getScore(i)
					<< "\n";
		}
	}
		break;
	case BatchMiraSearch:
		{
		//read in database
		filesystem::path dbfile(Database.c_str());
		GSSTreeSearcher gss(Verbose);

		if (!gss.load(dbfile))
		{
			cerr << "Could not read database " << Database << "\n";
			exit(-1);
		}
		double resolution = gss.getResolution();

		setDistance(ShapeDist);

		//read in each line of the batch file which should be
		//cmd in_ligand in_receptor(for DC Search)
		ifstream batch(Input.c_str());
		if (!batch)
		{
			cerr << "Could not read batch file " << Input << "\n";
			exit(-1);
		}

		string line;
		while (getline(batch, line))
		{
			stringstream toks(line);
			string cmd, mira;
			toks >> cmd;

			if (!toks)
				break;

			double less = 0, more = 0;
			unsigned k = 1;
			string mivfile, msvfile;

			if (cmd == "DCSearch")
			{
				toks >> mira;
				toks >> less;
				toks >> more;
			}
			else if (cmd == "NNSearch")
			{
				toks >> mira;
				toks >> k;
			}
			else if(cmd == "DCSearch2")
			{
				toks >> mira;
				toks >> mivfile;
				toks >> msvfile;
			}

			vector<double> times;

			if (ClearCacheFirst)
				clear_cache();

			ifstream mirain(mira.c_str());
			MiraObject obj;
			obj.read(mirain);
			ObjectTree objTree = gss.createTreeFromObject(obj);
			ObjectTree smallTree, bigTree;

			if (cmd == "DCSearch")
			{
				MGrid smallgrid, biggrid;
				objTree->makeGrid(smallgrid, resolution);
				biggrid = smallgrid;
				smallgrid.shrink(less);
				biggrid.grow(more);

				smallTree = ObjectTree(
						MappableOctTree::createFromGrid(smallgrid));
				bigTree = ObjectTree(MappableOctTree::createFromGrid(biggrid));
			}
			else if(cmd == "DCSearch2")
			{
				MiraObject miv, msv;
				ifstream msvin(msvfile.c_str());
				ifstream mivin(mivfile.c_str());
				msv.read(msvin);
				miv.read(mivin);

				smallTree = gss.createTreeFromObject(miv);
				bigTree = gss.createTreeFromObject(msv);
			}

			for (unsigned i = 0; i < TimeTrials; i++)
			{
				if (ClearCache)
				{
					clear_cache();
				}
				StringResults res;

				Timer t;

				if (cmd == "DCSearch" || cmd == "DCSearch2")
				{
					gss.dc_search(smallTree, bigTree, objTree, Verbose, res);
				}
				else if (cmd == "NNSearch")
				{
					if (k < 1) //scan
					{
						gss.nn_scan(objTree, Verbose, res);
					}
					else
					{
						gss.nn_search(objTree, k, HUGE_VAL, Verbose, res);
					}
				}
				else
				{
					cerr << "Illegal command " << cmd << " in batch file.\n";
					exit(-1);
				}

				//if verbose, have results
				for (unsigned i = 0, n = res.size(); i < n; i++)
				{
					cout << "RES " << mira << "\t" << res.getString(i) << "\t" << res.getScore(i)
							<< "\n";
				}

				times.push_back(t.elapsed());
			}

			cout << "Batch " << line;
			for (unsigned i = 0, n = times.size(); i < n; i++)
			{
				cout << " " << times[i];
			}
			cout << "\n";
		}
	}
		break;
	case NNSearch:
		{
		//read in database
		filesystem::path dbfile(Database.c_str());
		GSSTreeSearcher gss(Verbose);
		if (!gss.load(dbfile))
		{
			cerr << "Could not read database " << Database << "\n";
			exit(-1);
		}

		setDistance(ShapeDist);

		if (ExcludeMol.size() == 0)
		{
			//actually do a nearest neighbor search intelligently
			//using cutoffs within the tree search - this requires
			//a backed-in distance funciton (ignores shapedist)
			do_nnsearch(gss, IncludeMol, Output, K, NNThreshold);
		}
		else
		{
			//use shapedistance if two shapes are provided
			ResultMolecules res(SingleConformer);
			ObjectTree smallTree, bigTree;
			create_trees(gss, IncludeMol, ExcludeMol, LessDist, MoreDist,
					smallTree, bigTree);
			if (K < 1) //scan
			{
				gss.nn_scan(smallTree, bigTree, getres, res);
			}
			else //actually, also a scan
			{
				gss.nn_search(smallTree, bigTree, K, getres, res);
			}

			res.writeOutput(Output, Print);
		}
	}
		break;
	case DCSearch:
		{
		//read in database
		filesystem::path dbfile(Database.c_str());
		GSSTreeSearcher gss(Verbose);

		if (!gss.load(dbfile))
		{
			cerr << "Could not read database " << Database << "\n";
			exit(-1);
		}
		setDistance(ShapeDist);

		double dimension = gss.getDimension();
		double resolution = gss.getResolution();
		//read query molecule(s)
		if (IncludeMol.size() > 0 || ExcludeMol.size() > 0)
		{
			ObjectTree smallTree, bigTree;
			create_trees(gss, IncludeMol, ExcludeMol, LessDist, MoreDist,
					smallTree, bigTree);
			do_dcsearch(gss, smallTree, bigTree, Output);
		}
		else // range from single molecules
		{
			for (Molecule::iterator inmols(Input, dimension, resolution,
					KeepHydrogens, ProbeRadius); inmols; ++inmols)
			{
				ResultMolecules res(SingleConformer);
				//search
				//create bounding trees

				ObjectTree smallTree = gss.createTreeFromObject(*inmols,
						LessDist);
				ObjectTree bigTree = gss.createTreeFromObject(*inmols,
						-MoreDist);

				if (!ScanOnly)
					gss.dc_search(smallTree, bigTree, ObjectTree(), getres, res);

				if (ScanCheck || ScanOnly)
				{
					ResultMolecules res2(SingleConformer);
					gss.dc_scan_search(smallTree, bigTree, ObjectTree(), getres,
							res2);
					if (res2.size() != res.size())
					{
						cerr << "Scanning found different number\n";
					}
				}

				res.writeOutput(Output, Print);

			}
		}
	}
		break;
	case BatchSearch:
		{
		//read in database
		filesystem::path dbfile(Database.c_str());
		GSSTreeSearcher gss(Verbose);

		if (!gss.load(dbfile))
		{
			cerr << "Could not read database " << Database << "\n";
			exit(-1);
		}
		setDistance(ShapeDist);

		//read in each line of the batch file which should be
		//cmd in_ligand in_receptor(for DC Search)
		ifstream batch(Input.c_str());
		if (!batch)
		{
			cerr << "Could not read batch file " << Input << "\n";
			exit(-1);
		}

		string line;
		while (getline(batch, line))
		{
			stringstream toks(line);
			string cmd, ligand, receptor, output; //output always empty for batch
			toks >> cmd;

			if (!toks)
				break;

			double less = 0, more = 0;
			unsigned k = 1;

			if (cmd == "DCSearch")
			{
				toks >> ligand;
				toks >> receptor;
				toks >> less;
				toks >> more;
			}
			else if (cmd == "NNSearch")
			{
				toks >> ligand;
				toks >> k;
			}

			vector<double> times;
			ObjectTree smallTree, bigTree;
			if (cmd == "DCSearch")
				create_trees(gss, ligand, receptor, less, more, smallTree, bigTree);

			if (ClearCacheFirst)
				clear_cache();

			for (unsigned i = 0; i < TimeTrials; i++)
			{
				if (ClearCache)
				{
					clear_cache();
				}

				Timer t;

				if (cmd == "DCSearch")
				{
					do_dcsearch(gss, smallTree, bigTree, output);
				}
				else if (cmd == "NNSearch")
				{
					do_nnsearch(gss, ligand, output, k, HUGE_VAL);
				}
				else
				{
					cerr << "Illegal command " << cmd << " in batch file.\n";
					exit(-1);
				}

				times.push_back(t.elapsed());
			}

			cout << "Batch " << line;
			for (unsigned i = 0, n = times.size(); i < n; i++)
			{
				cout << " " << times[i];
			}
			cout << "\n";
		}
	}
		break;
	case BatchDB:
		{
		//setup where we run the same set of queries on a list of databases (provided in Database)
		//preproccess the queries to load in receptor/ligand structures

		//read in each line of the batch file which should be
		//cmd in_ligand in_receptor(for DC Search)
		ifstream batch(Input.c_str());
		if (!batch)
		{
			cerr << "Could not read batch file " << Input << "\n";
			exit(-1);
		}

		vector<QInfo> qinfos;
		string line;
		while (getline(batch, line))
		{
			stringstream toks(line);
			string cmd, ligand, receptor, output; //output always empty for batch
			toks >> cmd;

			if (!toks)
				break;
			if (cmd == "DCSearch")
			{
				double less = 0, more = 0;
				toks >> ligand;
				toks >> receptor;
				toks >> less;
				toks >> more;

				Molecule::iterator inmol(ligand, MaxDimension, Resolution,
						KeepHydrogens, ProbeRadius);
				Molecule inMol = *inmol;

				Molecule::iterator exmol(receptor, MaxDimension, Resolution,
						KeepHydrogens, ProbeRadius);
				Molecule exMol = *exmol;

				qinfos.push_back(QInfo(line, inMol, exMol, less, more));
			}
			else if (cmd == "NNSearch")
			{
				unsigned k = 1;
				toks >> ligand;
				toks >> k;

				Molecule::iterator inmol(ligand, MaxDimension, Resolution,
						KeepHydrogens, ProbeRadius);
				Molecule inMol = *inmol;

				qinfos.push_back(QInfo(line, inMol, k));
			}
			else
			{
				cerr << "Illegal command " << cmd << " in batch file.\n";
				exit(-1);
			}
		}

		//read in list of databases
		ifstream dbs(Database.c_str());
		if (!dbs)
		{
			cerr << "Could not read db file list " << Database << "\n";
			exit(-1);
		}

		while (getline(dbs, line))
		{
			filesystem::path dbfile(line.c_str());
			GSSTreeSearcher gss(Verbose);

			if (!gss.load(dbfile))
			{
				cerr << "Could not read database " << line << "\n";
				exit(-1);
			}
			setDistance(ShapeDist);

			if (gss.getResolution() != Resolution
					|| gss.getDimension() != MaxDimension)
			{
				cerr << "Resolution or dimension mismatch\n";
				exit(-1);
			}

			ResultMolecules res(SingleConformer);
			for (unsigned i = 0, n = qinfos.size(); i < n; i++)
			{
				cout << line << " " << qinfos[i].str << "\n";
				if (qinfos[i].cmd == NNSearch)
				{
					gss.nn_search(gss.createTreeFromObject(qinfos[i].in),
							qinfos[i].k, HUGE_VAL, false, res);
				}
				else
				{
					gss.dc_search(
							gss.createTreeFromObject(qinfos[i].in,
									qinfos[i].less),
							gss.createTreeFromObject(qinfos[i].ex,	qinfos[i].more, true),
							ObjectTree(),false, res);
				}
			}
		}

	}
		break;
	case SearchAllPointCombos:
		{
		//this will output a ton of files - it searches all 2^n combinations
		//of n interaction points using both DCSearch and NNSearch scanning

		//read in database
		filesystem::path dbfile(Database.c_str());
		GSSTreeSearcher gss(Verbose);

		if (!gss.load(dbfile))
		{
			cerr << "Could not read database " << Database << "\n";
			exit(-1);
		}
		setDistance(ShapeDist);

		double dimension = gss.getDimension();
		double resolution = gss.getResolution();
		//read query molecule(s) - must have both ligand and receptor
		if (IncludeMol.size() == 0 || ExcludeMol.size() == 0)
		{
			cerr << "SearchAllPointCombos requires both ligand and receptor\n";
			exit(-1);

		}

		if (!useInteractionPoints)
		{
			cerr << "SearchAllPointCombos requires interaction points\n";
			exit(-1);
		}
		double more = MoreDist;
		double less = LessDist;

		//read query molecule(s)
		Molecule::iterator imolitr(IncludeMol, dimension, resolution,
				KeepHydrogens, ProbeRadius);
		Molecule inMol = *imolitr;

		Molecule::iterator exmolitr(ExcludeMol, dimension, resolution,
				KeepHydrogens, ProbeRadius);
		Molecule exMol = *exmolitr;

		//create receptor tree

		//resize receptor constraint
		ObjectTree bigTree = gss.createTreeFromObject(exMol, more, true);

		//get the full shape of the ligand
		MappableOctTree *ligTree = MappableOctTree::create(dimension,
				resolution, inMol);
		MGrid lgrid;
		ligTree->makeGrid(lgrid, resolution);
		free(ligTree);
		ligTree = NULL;

		//shrink if requested
		if (less > 0)
		{
			lgrid.shrink(less);
		}
		//compute interaction points as single zero-radius points
		MGrid igrid(dimension, resolution);
		inMol.computeInteractionGridPoints(exMol, igrid, interactionDistance,
				interactionMaxClusterDist, interactionMinCluster, 0);

		vector<MGrid::Point> ipts;
		igrid.getSetPoints(ipts);
		unsigned max = 1 << ipts.size();

		//for each subset of interaction points (even empty - receptor only)
		for (unsigned i = 0; i < max; i++)
		{
			//compute small tree
			MGrid igrid(dimension, resolution);
			for (unsigned p = 0, np = ipts.size(); p < np; p++)
			{
				if ((1 << p) & i) //use it
				{
					if (interactionPointRadius == 0)
					{
						igrid.setPoint(ipts[p].x, ipts[p].y, ipts[p].z);
					}
					else
					{
						igrid.markXYZSphere(ipts[p].x, ipts[p].y, ipts[p].z,
								interactionPointRadius);
					}
				}
			}
			igrid &= lgrid;
			ObjectTree smallTree(MappableOctTree::createFromGrid(igrid), free);

			//do a scan for nn ranking
			if (NNSearchAll)
			{
				ResultMolecules res(SingleConformer);
				gss.nn_scan(smallTree, bigTree, true, res);
				stringstream outname;
				outname << Output << "_nn_" << i << ".sdf";
				ofstream out(outname.str().c_str());
				for (unsigned r = 0, n = res.size(); r < n; r++)
					res.writeSDF(out, r);

				out.close();
				//gzip output
				stringstream cmd;
				cmd << "gzip " << outname.str();
				::system(cmd.str().c_str());
			}

			//do a filter with DC
			stringstream dcoutname;
			dcoutname << Output << "_dc_" << i << ".sdf";
			do_dcsearch(gss, smallTree, bigTree, dcoutname.str());

		}
	}
		break;
	case MergeMira:
	{
		//read in a bunch of mira files and output the miv/msv
		vector<const MappableOctTree*> trees;
		for(unsigned i = 0, n = Files.size(); i < n; i++)
		{
			string file = Files[i];
			ifstream mirain(file.c_str());
			MiraObject obj;
			obj.read(mirain);
			MappableOctTree *tree = MappableOctTree::create(MaxDimension,
					Resolution, obj);

			trees.push_back(tree);
		}

		MappableOctTree *MIV = MappableOctTree::createFromIntersection(trees.size(), &trees[0]);
		MappableOctTree *MSV = MappableOctTree::createFromUnion(trees.size(), &trees[0]);

		string mivname = Output;
		mivname += "_miv.mira";
		string msvname = Output;
		msvname += "_msv.mira";

		ofstream mivout(mivname.c_str());
		ofstream msvout(msvname.c_str());

		MIV->dumpMiraGrid(mivout, Resolution);
		MSV->dumpMiraGrid(msvout, Resolution);

		//I am knowingly leaking memory here because I am lazy and this should not be in a loop
	}
		break;
	case MolGrid:
		{
		ofstream out(Output.c_str());

		if (useInteractionPoints)
		{
			if (IncludeMol.size() == 0 || ExcludeMol.size() == 0)
			{
				cerr
				<< "Need both ligand and receptor for interaction points.\n";
				exit(-1);
			}
			Molecule::iterator imolitr(IncludeMol, MaxDimension, Resolution,
					KeepHydrogens, ProbeRadius);
			Molecule ligand = *imolitr;

			Molecule::iterator exmolitr(ExcludeMol, MaxDimension, Resolution,
					KeepHydrogens, ProbeRadius);
			Molecule receptor = *exmolitr;

			MGrid igrid(MaxDimension, Resolution);
			ligand.computeInteractionGridPoints(receptor, igrid,
					interactionDistance, interactionMaxClusterDist,
					interactionMinCluster, interactionPointRadius);

			MappableOctTree *tree = MappableOctTree::create(MaxDimension,
					Resolution, ligand);
			MGrid lgrid;
			tree->makeGrid(lgrid, Resolution);
			lgrid.shrink(LessDist);

			lgrid &= igrid;

			free(tree);
			tree = MappableOctTree::createFromGrid(lgrid);

			dumpTree(tree, out);
			free(tree);
		}
		else
		{
			for (Molecule::iterator mitr(Input, MaxDimension, Resolution,
					KeepHydrogens, ProbeRadius); mitr; ++mitr)
			{
				MappableOctTree *tree = MappableOctTree::create(MaxDimension,
						Resolution, *mitr);

				if (LessDist >= 0)
				{
					MGrid grid;
					tree->makeGrid(grid, Resolution);

					grid.shrink(LessDist);
					grid.grow(MoreDist);
					free(tree);
					tree = MappableOctTree::createFromGrid(grid);
				}

				dumpTree(tree, out);
				free(tree);
			}
		}
		break;
	}
	}
	return 0;
}

