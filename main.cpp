/*
 Pharmer: Efficient and Exact 3D Pharmacophore Search
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
 *  Created on: Jul 30, 2010
 *      Author: dkoes
 *
 *  Pharmacophore tools for recognition and search.
 *  Really these would probably more suited as their own executables,
 *  but I'm too lazy to move away from Eclipe's default make process which
 *  requires a single target.
 */

#include "CommandLine2/CommandLine.h"
#include "pharmarec.h"
#include "pharmerdb.h"
#include "PharmerQuery.h"
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include "queryparsers.h"
#include "Timer.h"
#include "MolFilter.h"
#include "PharmerServer.h"
#include <sys/types.h>
#include <sys/wait.h>
#include "ReadMCMol.h"
#include "Excluder.h"

using namespace boost;
using namespace OpenBabel;

cl::opt<bool> SeparateWeight("separate-weight",
		cl::desc("Segregate database based on molecular weight"),
		cl::init(false));
cl::opt<bool> Quiet("q", cl::desc("quiet; suppress informational messages"),
		cl::init(true));
cl::opt<bool> ShowQuery("show-query", cl::desc("print query points"),
		cl::init(false));
cl::opt<bool> Print("print", cl::desc("print results"), cl::init(true));
cl::opt<string> Cmd("cmd",
		cl::desc("command [pharma, dbcreate, dbsearch, server]"),
		cl::Positional);
cl::list<string> Database("dbdir", cl::desc("database directory(s)"));
cl::list<string> inputFiles("in", cl::desc("input file(s)"));
cl::list<string> outputFiles("out", cl::desc("output file(s)"));

cl::opt<string> pharmaSpec("pharmaspec",
		cl::desc("pharmacophore specification"));
cl::opt<unsigned> NThreads("nthreads", cl::desc("utilize n threads; default 1"),
		cl::value_desc("n"), cl::init(1));
cl::opt<double> MaxRMSD("max-rmsd",
		cl::desc("maximum allowed RMSD; default max allowed by query"),
		cl::init(HUGE_VAL));
cl::opt<unsigned> MinWeight("min-weight",
		cl::desc("minimum allowed molecular weight"), cl::init(0));
cl::opt<unsigned> MaxWeight("max-weight",
		cl::desc("maximum allowed molecular weight"), cl::init(UINT_MAX));
cl::opt<unsigned> MinNRot("min-nrot",
		cl::desc("minimum allowed rotatable bonds"), cl::init(0));
cl::opt<unsigned> MaxNRot("max-nrot",
		cl::desc("maximum allowed rotatable bonds"), cl::init(UINT_MAX));

cl::opt<unsigned> ReduceConfs("reduceconfs",
		cl::desc("return at most n conformations for each molecule"),
		cl::value_desc("n"), cl::init(0));
cl::opt<unsigned> MaxOrient("max-orient",
		cl::desc("return at most n orientations of each conformation"),
		cl::value_desc("n"), cl::init(UINT_MAX));
cl::opt<unsigned> MaxHits("max-hits", cl::desc("return at most n results"),
		cl::value_desc("n"), cl::init(UINT_MAX));
cl::opt<unsigned> Port("port", cl::desc("port for server to listen on"));
cl::opt<string> LogDir("logdir", cl::desc("log directory for server"),
		cl::init("."));
cl::opt<bool> ExtraInfo("extra-info",
		cl::desc("Output additional molecular properties.  Slower."),
		cl::init(false));
cl::opt<bool> SortRMSD("sort-rmsd", cl::desc("Sort results by RMSD."),
		cl::init(false));
cl::opt<bool> FilePartition("file-partition",
		cl::desc("Partion database slices based on files"), cl::init(false));

cl::opt<string> MinServer("min-server",cl::desc("minimization server address"));
cl::opt<unsigned> MinPort("min-port",cl::desc("port for minimization server"));

cl::opt<string> Receptor("receptor",
		cl::desc("Receptor file for interaction pharmacophroes"));

typedef void (*pharmaOutputFn)(ostream&, vector<PharmaPoint>&, Excluder& excluder);

static void pharmaNoOutput(ostream&, vector<PharmaPoint>&, Excluder& excluder)
{
}

static void pharmaTxtOutput(ostream& out, vector<PharmaPoint>& points, Excluder& excluder)
{
	for (unsigned i = 0, n = points.size(); i < n; i++)
		out << points[i] << "\n";
}

static void pharmaJSONOutput(ostream& out, vector<PharmaPoint>& points, Excluder& excluder)
{
	Json::Value root;
	convertPharmaJson(root, points);
	excluder.addToJSON(root);
	Json::StyledStreamWriter writer;
	writer.write(out, root);
}

static void pharmaSDFOutput(ostream& out, vector<PharmaPoint>& points, Excluder& excluder)
{
	OBConversion conv;
	conv.SetOutFormat("SDF");
	OBMol mol;

	for (int p = 0, np = points.size(); p < np; p++)
	{
		OBAtom *atom = mol.NewAtom();
		atom->SetAtomicNum(points[p].pharma->atomic_number_label);
		atom->SetVector(points[p].x, points[p].y, points[p].z);
	}

	conv.Write(&mol, &out);
}

//identify all the pharma points within each mol in the input file
static void handle_pharma_cmd(const Pharmas& pharmas)
{

	if (outputFiles.size() > 0 && outputFiles.size() != inputFiles.size())
	{
		cerr << "Number of outputs must equal number of inputs.\n";
		exit(-1);
	}

	for (unsigned i = 0, n = inputFiles.size(); i < n; i++)
	{
		string fname = inputFiles[i];
		ifstream in(fname.c_str());
		if (!in)
		{
			cerr << "Error reading file " << fname << "\n";
			exit(-1);
		}

		ofstream out;
		pharmaOutputFn outfn = pharmaNoOutput;
		//output can be plain text, json, or an sdf file (collection of atoms)
		if (outputFiles.size() > 0)
		{
			out.open(outputFiles[i].c_str());
			if (!out)
			{
				cerr << "Error opening output file " << outputFiles[i] << "\n";
				exit(-1);
			}
			string ext = filesystem::extension(outputFiles[i]);
			if (ext == ".txt" || ext == "")
				outfn = pharmaTxtOutput;
			else if (ext == ".json")
				outfn = pharmaJSONOutput;
			else if (ext == ".sdf")
				outfn = pharmaSDFOutput;
			else
			{
				cerr << "Unsupported output format\n";
				exit(-1);
			}
		}

		//special case - convert  query to points
		if (filesystem::extension(fname) == ".json"
				|| filesystem::extension(fname) == ".ph4"
				|| filesystem::extension(fname) == ".query"
				|| filesystem::extension(fname) == ".txt"
				|| filesystem::extension(fname) == ".pml")
		{
			ifstream in(fname.c_str());
			vector<PharmaPoint> points;
			Excluder excluder;

			if (filesystem::extension(fname) == ".json"
					|| filesystem::extension(fname) == ".query")
			{
				JSonQueryParser parser;
				parser.parse(pharmas, in, points, excluder);
			}
			else if (filesystem::extension(fname) == ".ph4")
			{
				PH4Parser parser;
				parser.parse(pharmas, in, points, excluder);
			}
			else if (filesystem::extension(fname) == ".pml")
			{
				PMLParser parser;
				parser.parse(pharmas, in, points, excluder);
			}
			else
			{
				TextQueryParser parser;
				parser.parse(pharmas, in, points, excluder);
			}
			outfn(out, points, excluder);
		}
		else //pharma recognition
		{
			OBConversion conv;
			OBFormat *format = conv.FormatFromExt(fname.c_str());
			if (format == NULL)
			{
				cerr << "Invalid input file format " << fname << "\n";
				exit(-1);
			}
			conv.SetInFormat(format);
			OBMol mol;
			vector<PharmaPoint> points;

			OBMol receptor;
			if (Receptor.size() > 0)
			{
				OBConversion rconv;
				OBFormat *rformat = rconv.FormatFromExt(Receptor.c_str());
				if (format)
				{
					rconv.SetInFormat(rformat);
					ifstream rin(Receptor.c_str());
					rconv.Read(&receptor, &rin);
				}
			}

        		OBAromaticTyper aromatics;
        		OBAtomTyper atyper;
			while (conv.Read(&mol, &in))
			{
				//perform exactly the same analyses as dbcreate
        			mol.AddHydrogens();

        			mol.FindRingAtomsAndBonds();
				mol.FindChiralCenters();
				mol.PerceiveBondOrders();
				aromatics.AssignAromaticFlags(mol);
				mol.FindSSSR();
				atyper.AssignTypes(mol);
				atyper.AssignHyb(mol);

				if (receptor.NumAtoms() > 0)
				{
					vector<PharmaPoint> screenedout;
					getInteractionPoints(pharmas, receptor, mol, points,
							screenedout);
				}
				else
					getPharmaPoints(pharmas, mol, points);
				Excluder dummy;
				if (!Quiet)
					pharmaTxtOutput(cout, points, dummy);
				outfn(out, points, dummy);
			}
		}
	}
}

struct FilterDBCreate
{
	WeightRangeFilter filter;
	shared_ptr<PharmerDatabaseCreator> db;

	FilterDBCreate()
	{
	}
	FilterDBCreate(double min, double max, PharmerDatabaseCreator* d) :
			filter(min, max), db(d)
	{
	}
};
//create a database
static void handle_dbcreate_cmd(const Pharmas& pharmas)
{
	if (Database.size() == 0)
	{
		cerr
				<< "Need to specify location of database directory to be created.\n";
		exit(-1);
	}

	//create directories
	for (unsigned i = 0, n = Database.size(); i < n; i++)
	{
		if (!filesystem::create_directory(Database[i]))
		{
			cerr << "Unable to create database directory " << Database[i]
					<< "\n";
			exit(-1);
		}
	}

	//check and setup input
	OBConversion conv;
	unsigned long long numBytes = 0;
	for (unsigned i = 0, n = inputFiles.size(); i < n; i++)
	{

		if (!filesystem::exists(inputFiles[i].c_str()))
		{
			cerr << "Invalid input file: " << inputFiles[i] << "\n";
			exit(-1);
		}
		OBFormat *format = conv.FormatFromExt(inputFiles[i].c_str());
		if (format == NULL)
		{
			cerr << "Invalid input format: " << inputFiles[i] << "\n";
			exit(-1);
		}

		numBytes += filesystem::file_size(inputFiles[i]);
	}

	double weightThresholds[] =
	{ 0, 321, 351, 376, 401, 426, 451, 476, 501, HUGE_VAL };
	unsigned nweights = sizeof(weightThresholds) / sizeof(double);

	if (!SeparateWeight)
	{
		nweights = 2;
		weightThresholds[1] = HUGE_VAL;
	}
	//create databases
	//openbabel can't handled multithreaded reading, so we actually have to fork off a process
	//for each database
	for (unsigned d = 0, nd = Database.size(); d < nd; d++)
	{
		if (nd == 1 || fork() == 0)
		{
			//split by weight
			vector<FilterDBCreate> dbs;
			for (unsigned w = 1; w < nweights; w++)
			{
				double min = weightThresholds[w - 1];
				double max = weightThresholds[w];
				WeightRangeFilter filter(min, max);
				string dname = Database[d] + "/w" + lexical_cast<string>(min);
				if (!filesystem::create_directory(dname))
				{
					cerr << "Unable to create database directory " << dname
							<< "\n";
					exit(-1);
				}

				dbs.push_back(
						FilterDBCreate(min, max,
								new PharmerDatabaseCreator(pharmas, dname,
										NThreads)));
			}

			//now read files
			unsigned long readBytes = 0;
			for (unsigned i = 0, n = inputFiles.size(); i < n; i++)
			{
				if (FilePartition)
				{
					if (i % nd != d)
						continue;
				}
				ifstream in(inputFiles[i].c_str());
				OBFormat *format = conv.FormatFromExt(inputFiles[i].c_str());
				if (!Quiet)
					cout << "Adding " << inputFiles[i] << "\n";
				unsigned stride = nd;
				unsigned offset = d;
				if (FilePartition)
				{
					stride = 1;
					offset = 0;
				}
				ReadMCMol reader(in, format, stride, offset, ReduceConfs);
				OBMol mol;

				while (reader.read(mol))
				{
					BOOST_FOREACH(FilterDBCreate& f, dbs)
					{
						if (f.filter.skip(mol))
							continue;
						f.db->addMolToDatabase(mol, mol.GetMolWt());
						break;
					}
				}
				BOOST_FOREACH(FilterDBCreate& f, dbs)
				{
					f.db->writeStats();
				}
				readBytes += filesystem::file_size(inputFiles[i]);
			}

			BOOST_FOREACH(FilterDBCreate& f, dbs)
			{
				f.db->createSpatialIndex();
			}
			exit(0);
		}
	}

	int status;
	while (wait(&status) > 0)
	{
		if (!WIFEXITED(status) && WEXITSTATUS(status) != 0)
			abort();
		continue;
	}

}

//thread class for loading database info
struct LoadDatabase
{
	unsigned totalConf;
	unsigned totalMols;

	LoadDatabase() :
			totalConf(0), totalMols(0)
	{

	}

	void operator()(vector<vector<MolWeightDatabase> >& databases, unsigned i,
			filesystem::path dbpath)
	{
		//look for sub directories, assume weight divided, have to sort first
		vector<unsigned> weights;
		for (filesystem::directory_iterator itr(dbpath), end_itr;
				itr != end_itr; ++itr)
		{
			if (is_directory(itr->status())
					&& filesystem::exists(itr->path() / "info"))
			{
				filesystem::path subdir = itr->path();
				int w;
				sscanf(subdir.filename().c_str(), "w%d", &w);
				weights.push_back(w);
			}
		}

		sort(weights.begin(), weights.end());

		for (unsigned j = 0, nw = weights.size(); j < nw; j++)
		{
			unsigned w = weights[j];
			filesystem::path subdir = dbpath
					/ (string("w" + lexical_cast<string>(w)));
			shared_ptr<PharmerDatabaseSearcher> db(
					new PharmerDatabaseSearcher(subdir));
			if (j > 0)
				databases[i].back().max = w;
			databases[i].push_back(MolWeightDatabase(db, w));

			if (!db->isValid())
			{
				cerr << "Error reading database " << Database[i] << "\n";
				exit(-1);
			}
			totalConf += db->numConformations();
			totalMols += db->numMolecules();
		}
		if(databases[i].size() > 0) //otherwise we were passed invalid dir
			databases[i].back().max = HUGE_VAL;
	}
};
static void loadDatabases(vector<vector<MolWeightDatabase> >& databases,
		unsigned& totalConf, unsigned& totalMols)
{
	totalConf = 0;
	totalMols = 0;
	databases.resize(Database.size());
	vector<LoadDatabase> loaders(Database.size());
	thread_group loading_threads;
	for (unsigned i = 0, n = Database.size(); i < n; i++)
	{
		filesystem::path dbpath(Database[i]);
		if (!filesystem::is_directory(dbpath))
		{
			cerr << "Invalid database directory path: " << Database[i] << "\n";
			exit(-1);
		}

		loading_threads.add_thread(
				new thread(ref(loaders[i]), ref(databases), i, dbpath));
	}
	loading_threads.join_all();

	BOOST_FOREACH(const LoadDatabase& ld, loaders)
	{
		totalConf += ld.totalConf;
		totalMols += ld.totalMols;
	}
}

//search the database
static void handle_dbsearch_cmd()
{
	if (pharmaSpec.size() != 0)
	{
		cerr
				<< "Warning: pharmaspec option not valid for database search; database specification will be used instead\n";
	}

	//databases
	if (Database.size() == 0)
	{
		cerr << "Require database directory path.\n";
		exit(-1);
	}

	vector<vector<MolWeightDatabase> > databases(Database.size());
	unsigned totalC = 0;
	unsigned totalM = 0;
	loadDatabases(databases, totalC, totalM);

	if (!Quiet)
		cout << "Searching " << totalC << " conformations of " << totalM
				<< " compounds.\n";

	//query parameters
	QueryParameters params;
	params.maxRMSD = MaxRMSD;
	params.minWeight = MinWeight;
	params.maxWeight = MaxWeight;
	params.minRot = MinNRot;
	params.maxRot = MaxNRot;
	params.reduceConfs = ReduceConfs;
	params.orientationsPerConf = MaxOrient;
	params.maxHits = MaxHits;
	if (SortRMSD)
		params.sort = SortType::RMSD;

	//data parameters
	DataParameters dparams;
	dparams.extraInfo = ExtraInfo;
	dparams.sort = SortType::RMSD;

	//query file
	if (inputFiles.size() < 1)
	{
		cerr << "Need input pharmacophore query file(s).\n";
		exit(-1);
	}

	Timer timer;

	if (outputFiles.size() > 0 && outputFiles.size() != inputFiles.size())
	{
		cerr << "Number of outputs must equal number of inputs\n";
		exit(-1);
	}

	for (unsigned i = 0, n = inputFiles.size(); i < n; i++)
	{

		if (!PharmerQuery::validFormat(filesystem::extension(inputFiles[i])))
		{
			cerr << "Invalid extension for query file: " << inputFiles[i]
					<< "\n";
			exit(-1);
		}
		cout << "Query " << inputFiles[i] << "\n";
		ifstream qfile(inputFiles[i].c_str());
		if (!qfile)
		{
			cerr << "Could not open query file: " << inputFiles[i] << "\n";
			exit(-1);
		}

		PharmerQuery query(databases, qfile,
				filesystem::extension(inputFiles[i]), params,
				NThreads * databases.size());

		string err;
		if (!query.isValid(err))
		{
			cerr << err << "\n";
			exit(-1);
		}
		if (ShowQuery)
			query.print(cout);

		query.execute(); //blocking

		if (Print) //dump to stdout
		{
			query.outputData(dparams, cout);
		}

		//output file
		if (outputFiles.size() > 0)
		{
			string outname = outputFiles[i];
			string oext = filesystem::extension(outname);
			ofstream out;

			if (oext != ".sdf" && oext != ".txt" && oext != "" && oext != ".gz")
			{
				cerr << "Invalid output format.  Support only .sdf and .txt\n";
				exit(-1);
			}
			out.open(outname.c_str());
			if (!out)
			{
				cerr << "Could not open output file: " << outname << "\n";
				exit(-1);
			}

			if(oext == ".gz") //assumed to be compressed sdf
			{
				boost::iostreams::filtering_ostream gzout;
				gzout.push(boost::iostreams::gzip_compressor());
				gzout.push(out);
				query.outputMols(gzout);
			}
			else if (oext != ".sdf") //text output
			{
				query.outputData(dparams, out);
			}
			else //mol output
			{
				query.outputMols(out);
			}
		}

		cout << "NumResults: " << query.numResults() << "\n";
	}

	cout << "Time: " << timer.elapsed() << "\n";
}

int main(int argc, char *argv[])
{
	cl::ParseCommandLineOptions(argc, argv);
	obErrorLog.StopLogging(); //just slows us down, possibly buggy?
	//if a pharma specification file was given, load that into the global Pharmas
	Pharmas pharmas(defaultPharmaVec);
	if (pharmaSpec.size() > 0)
	{
		ifstream pharmin(pharmaSpec.c_str());
		if (!pharmas.read(pharmin))
		{
			cerr << "Invalid pharmacophore specification file.\n";
			exit(-1);
		}
	}

	if (Cmd == "pharma")
	{
		handle_pharma_cmd(pharmas);
	}
	else if(Cmd == "showpharma")
	{
		pharmas.write(cout);
	}
	else if (Cmd == "dbcreate")
	{
		handle_dbcreate_cmd(pharmas);
	}
	else if (Cmd == "dbsearch")
	{
		handle_dbsearch_cmd();
	}
	else if (Cmd == "server")
	{
		vector<vector<MolWeightDatabase> > databases;
		unsigned totalC = 0;
		unsigned totalM = 0;
                //total hack time - fcgi uses select which can't
                //deal with file descriptors higher than 1024, so let's reserve some
                #define MAXRESERVEDFD (SERVERTHREADS*2) 
                int reservedFD[MAXRESERVEDFD] = {0,};
                for(unsigned i = 0; i < MAXRESERVEDFD; i++)
                {
                        reservedFD[i] = open("/dev/null",O_RDONLY);
                }
                //loadDatabases will open a whole bunch of files
                loadDatabases(databases, totalC, totalM);
                //now free reserved fds
                for(unsigned i = 0; i < MAXRESERVEDFD; i++)
                {
                        close(reservedFD[i]);
                }
		pharmer_server(Port, databases, LogDir, totalC, totalM, MinServer, MinPort);
	}
	else
	{
		cl::PrintHelpMessage();
		if (Cmd.size() == 0)
			cerr << "Command [pharma, dbcreate, dbsearch] required.\n";
		else
			cerr << Cmd << " not a valid command.\n";
		exit(-1);
	}

	return 0;
}

