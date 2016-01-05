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
#include <ShapeConstraints.h>
#include "ReadMCMol.h"
#include "dbloader.h"
#include "GninaConverter.h"

using namespace boost;
using namespace OpenBabel;


cl::opt<bool> Quiet("q", cl::desc("quiet; suppress informational messages"),
		cl::init(true));
cl::opt<bool> ShowQuery("show-query", cl::desc("print query points"),
		cl::init(false));
cl::opt<bool> Print("print", cl::desc("print results"), cl::init(true));
cl::opt<string> Cmd("cmd",
		cl::desc("command [pharma, dbcreate, dbcreateserverdir, dbsearch, server]"),
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
cl::opt<unsigned> Port("port", cl::desc("port for server to listen on"),cl::init(17000));
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

cl::opt<string> Prefixes("prefixes",
		cl::desc("[dbcreateserverdir,server] File of directory prefixes to use for striping."));
cl::opt<string> DBInfo("dbinfo",
		cl::desc("[dbcreateserverdir] JSON file describing database subset"));
cl::opt<string> Ligands("ligs", cl::desc("[dbcreateserverdir] Text file listing locations of molecules"));
cl::opt<bool> NoIndex("noindex",cl::desc("[dbcreateserverdir] Do not create indices"), cl::init(false));

typedef void (*pharmaOutputFn)(ostream&, vector<PharmaPoint>&, ShapeConstraints& excluder);

static void pharmaNoOutput(ostream&, vector<PharmaPoint>&, ShapeConstraints& excluder)
{
}

static void pharmaTxtOutput(ostream& out, vector<PharmaPoint>& points, ShapeConstraints& excluder)
{
	for (unsigned i = 0, n = points.size(); i < n; i++)
		out << points[i] << "\n";
}

static void pharmaJSONOutput(ostream& out, vector<PharmaPoint>& points, ShapeConstraints& excluder)
{
	Json::Value root;
	convertPharmaJson(root, points);
	excluder.addToJSON(root);
	Json::StyledStreamWriter writer;
	writer.write(out, root);
}

static void pharmaSDFOutput(ostream& out, vector<PharmaPoint>& points, ShapeConstraints& excluder)
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
			ShapeConstraints excluder;

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
				ShapeConstraints dummy;
				if (!Quiet)
					pharmaTxtOutput(cout, points, dummy);
				outfn(out, points, dummy);
			}
		}
	}
}


//there was a bug with generating sminaData due to a bug in openbabel's handling
//of multi-conformer OBMols.  This regenerates the sminaData/Index.
//Takes the path to sminaData/molData/sminaIndex
static void handle_fixsmina_cmd()
{
  for (unsigned i = 0, n = inputFiles.size(); i < n; i++)
  {
    filesystem::path p(inputFiles[i]);
    p = p.parent_path();
    filesystem::path sd = p / "sminaData";
    filesystem::path si = p / "sminaIndex";
    filesystem::path md = p / "molData";

    if(!filesystem::exists(sd)) {
      cerr << sd <<" does not exist!\n";
      continue;
    }
    if(!filesystem::exists(si)) {
      cerr << sd <<" does not exist!\n";
      continue;
    }
    if(!filesystem::exists(md)) {
      cerr << sd <<" does not exist!\n";
      continue;
    }


    //read sminaIndex into memory
    vector<  pair<unsigned long, unsigned long> > sminaIndex;
    sminaIndex.resize(filesystem::file_size(si)/sizeof( pair<unsigned long, unsigned long>));
    FILE *sif = fopen(si.c_str(), "r");
    size_t ret = fread(&sminaIndex[0], sizeof(pair<unsigned long, unsigned long>), sminaIndex.size(), sif);
    if(ret == 0) cerr << "Could not read " << si << "\n";
    fclose(sif);

    //memmap molData
    MMappedRegion<unsigned char> molData;
    molData.map(md.string(), true, true);

    //open sminaData for writing (blows away current file)
    FILE *sminaData = fopen(sd.c_str(), "w");

    //setup variables for loop
    PMolReaderSingleAlloc pread;
    OBConversion obconv;
    obconv.SetInAndOutFormats("SDF","SDF");
    OBMol mol;

    //for each mol
    for(unsigned i = 0, n = sminaIndex.size(); i < n; i++)
    {
      //read from molData
      unsigned long molpos = sminaIndex[i].first;
      MolData mdata;
      mdata.read(molData.begin(), molpos, pread);
      stringstream str;
      mdata.mol->writeSDF(str);

      obconv.ReadString(&mol, str.str());
      mol.AddHydrogens();
      mol.SetAutomaticFormalCharge(false);
      DeleteHydrogens(mol); //leaves just polars

      //write out updated sminaData
      stringstream data;
      GninaConverter::MCMolConverter mconv(mol);
      mconv.convertConformer(0, data);
      unsigned long smpos = ftell(sminaData); //location in smina dta
      unsigned sz = data.str().size();

      //output actual data
      fwrite(&sz, sizeof(sz),1, sminaData);
      fwrite(data.str().c_str(), sizeof(char), sz, sminaData);
      //update sminaIndex
      sminaIndex[i].second = smpos; //store location
    }
    //write out updated sminaIndex
    sif = fopen(si.c_str(), "w");
    fwrite(&sminaIndex[0], sizeof(pair<unsigned long, unsigned long>), sminaIndex.size(), sif);
  }
}

//take a file as input and write out the same structurs only with
//embedded pharmacophore data; input and output can be the same files
//temporary result is kept in memory, so file should not be large
static void handle_phogrify_cmd(const Pharmas& pharmas)
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

		stringstream out;
		OBConversion conv;
		OBFormat *format = conv.FormatFromExt(fname.c_str());
		if (format == NULL)
		{
			cerr << "Invalid input file format " << fname << "\n";
			exit(-1);
		}
		conv.SetInFormat(format);
		conv.SetOutFormat(format);

		OBMol mol;
		vector<PharmaPoint> points;


		OBAromaticTyper aromatics;
		OBAtomTyper atyper;
		conv.SetInStream(&in);
		conv.SetOutStream(&out);
		while (conv.Read(&mol))
		{
			OBMol origmol = mol;

			//perform exactly the same analyses as dbcreate
			mol.AddHydrogens();

			mol.FindRingAtomsAndBonds();
			mol.FindChiralCenters();
			mol.PerceiveBondOrders();
			aromatics.AssignAromaticFlags(mol);
			mol.FindSSSR();
			atyper.AssignTypes(mol);
			atyper.AssignHyb(mol);

			getPharmaPoints(pharmas, mol, points);
			stringstream phdata;
			for (unsigned i = 0, n = points.size(); i < n; i++)
				phdata << points[i] << "\n";

			OBPairData* sddata = new OBPairData();
			sddata->SetAttribute("pharmacophore");
			sddata->SetValue(phdata.str());
			origmol.DeleteData("pharmacophore"); //replace
			origmol.SetData(sddata);
			conv.Write(&origmol);
		}

		in.close();
		if(out.str().length() == 0)
		{
			//make sure we don't overwrite original file if we didn't generate output for some reason
			cerr << "Error generating output for: " << fname << "\n";
			continue;
		}

		string outname = outputFiles[i];
		string oext = filesystem::extension(outname);
		ofstream outf;

		outf.open(outname.c_str());
		if (!outf)
		{
			cerr << "Could not open output file: " << outname << "\n";
			exit(-1);
		}


		if(oext == ".gz") //assumed to be compressed sdf
		{
			boost::iostreams::filtering_ostream gzout;
			gzout.push(boost::iostreams::gzip_compressor());
			gzout.push(outf);
			gzout.write(out.str().c_str(), out.str().size());
		}
		else
		{
			outf.write(out.str().c_str(), out.str().size());
		}
	}
}

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

	//create databases
	//openbabel can't handled multithreaded reading, so we actually have to fork off a process
	//for each database
	for (unsigned d = 0, nd = Database.size(); d < nd; d++)
	{
		if (nd == 1 || fork() == 0)
		{
			Json::Value blank;
			PharmerDatabaseCreator db(pharmas, Database[d], blank);
			unsigned long uniqueid = 1;
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
					db.addMolToDatabase(mol, uniqueid*nd+d, mol.GetTitle());
					uniqueid++;
				}

				db.writeStats();
				readBytes += filesystem::file_size(inputFiles[i]);
			}

			db.createSpatialIndex();
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

//read in ligand file names and verify the files exist
struct LigandInfo
{
	filesystem::path file;
	long id;
	string name;

	LigandInfo(): id(0) {}
};

static void signalhandler(int sig)
{
  //ignore
}


//create a database directory within the server framework
//in this framework we provide a file of prefixes where each line is
//a location (on a different hard drive) for creating a strip of the overall database
//for input molecules, we provide a file where each line is the location of the conformers
//if a single molecule in an sdf.gz file, also included on the line are the unique id of the molecule and the
//space delimited possible names for that molecule
//we also specify a database description file which is a json file describing the database
//the json object is indexed by database key; the keys define the subdirectory name to use in prefixes
static void handle_dbcreateserverdir_cmd(const Pharmas& pharmas)
{
	signal(SIGUSR1, signalhandler); //don't let ourselves get interrupted by build signals

	ifstream prefixes(Prefixes.c_str());
	if (Prefixes.size() == 0 || !prefixes)
	{
		cerr << "Problem with prefixes.\n";
		exit(-1);
	}

	ifstream dbinfo(DBInfo.c_str());
	if(DBInfo.size() == 0 || !dbinfo)
	{
		cerr << "Problem with database info.\n";
		exit(-1);
	}

	ifstream ligs(Ligands.c_str());
	if(Ligands.size() == 0 || !filesystem::exists(Ligands))
	{
		cerr << "Need ligand file\n";
		exit(-1);
	}

	//parse dbinfo into json
	Json::Value root; // will contains the root value after parsing.
	Json::Reader reader;
	if(!reader.parse(dbinfo, root)) {
		cerr << "Error reading database info JSON\n";
		exit(-1);
	}


	vector<LigandInfo> liginfos;
	string line;
	while(getline(ligs,line))
	{
		stringstream str(line);
		LigandInfo info;

		str >> info.file;

		if(!filesystem::exists(info.file))
		{
			cerr << "File " << info.file << " does not exist\n";
		}
		str >> info.id;
		if(info.id < 0)
		{
			cerr << "Error in ligand file on line:\n" << line << "\n";
			exit(-1);
		}

		getline(str, info.name); //get rest as name
		liginfos.push_back(info);
	}

	//get key for database, this is the name of the subdir
	if(!root.isMember("subdir"))
	{
		cerr << "Database info needs subdir field.";
	}
	stringstream key;
	key << root["subdir"].asString() << "-" << time(NULL);
	string subset = root["subdir"].asString();

	//check for splits
	vector<string> splits;
	if(root.isMember("splitdirs") && root["splitdirs"].isArray())
	{
		Json::Value sdirs = root["splitdirs"];
		for(unsigned i = 0, n = sdirs.size(); i < n; i++)
		{
			splits.push_back(sdirs[i].asString());
		}
	}

	//create directories
	vector<filesystem::path> directories;
	vector< vector<filesystem::path> > unflatdirectories;
	vector<filesystem::path> topdirectories; //same as directories unless there are splits
	vector<filesystem::path> symlinks;
	if(!prefixes) {
		cerr << "Error with prefixes file\n";
		exit(-1);
	}

	string fname;
	while(getline(prefixes,fname))
	{
		if(fname.length() > 0 && fname[0] == '#') //comment
		{
			continue;
		}
		if(filesystem::exists(fname))
		{
			filesystem::path dbpath(fname);
			dbpath /= key.str();

			topdirectories.push_back(dbpath);
			unflatdirectories.push_back(vector<filesystem::path>()); //needs splits to be dispersed
			if(splits.size() > 0) //create a bunch of subdirs
			{
				for(unsigned i = 0, n = splits.size(); i < n; i++)
				{
					filesystem::path splitdb = dbpath / splits[i];
					filesystem::create_directories(splitdb);
					unflatdirectories.back().push_back(splitdb);
				}
			}
			else
			{
				filesystem::create_directories(dbpath);
				unflatdirectories.back().push_back(dbpath);
			}

			//single link for the whole thing
			filesystem::path link(fname);
			link /= subset;
			symlinks.push_back(link);
		}
		else
		{
			cerr << "Prefix " << fname << " does not exist\n";
		}
	}

	//flatten directories
	while(true)
	{
		unsigned i = 0, n = unflatdirectories.size();
		for(; i < n; i++)
		{
			if(unflatdirectories[i].size() == 0)
				break;
			directories.push_back(unflatdirectories[i].back());
			unflatdirectories[i].pop_back();
		}
		if(i != n)
			break;
	}

	//portion memory between processes
	unsigned long memsz = sysconf (_SC_PHYS_PAGES) * sysconf (_SC_PAGESIZE);
	memsz /= directories.size();
	memsz /= 2; //only take half of available memory

	unsigned maxt = thread::hardware_concurrency();
	unsigned numrunning = 0;
	//multi-thread (fork actually, due to openbabel) across all prefixes
	//create databases
	//openbabel can't handled multithreaded reading, so we actually have to fork off a process
	//for each database
	for (unsigned d = 0, nd = directories.size(); d < nd; d++)
	{
		if (d == (nd-1) || fork() == 0)
		{
			{ //this is an empty scope created to make sure db's destructor is called before exit
				PharmerDatabaseCreator db(pharmas, directories[d], root);
				db.setInMemorySize(memsz);
				OBConversion conv;

				//now read files
				for (unsigned i = 0, n = liginfos.size(); i < n; i++)
				{
					if( (i%nd) == d )
					{ //part of our slice
						const LigandInfo info = liginfos[i];
						string name = info.file.string();
						//openbabel's builtin zlib reader seems to use increasing amounts
						//of memory over time, so use boost's
						ifstream *uncompressed_inmol = new std::ifstream(name.c_str());
						iostreams::filtering_stream<iostreams::input> *inmol = new iostreams::filtering_stream<iostreams::input>();

						std::string::size_type pos = name.rfind(".gz");
						if (pos != std::string::npos)
						{
							inmol->push(iostreams::gzip_decompressor());
						}
						inmol->push(*uncompressed_inmol);


						OBFormat *format = conv.FormatFromExt(info.file.c_str());

						if(format != NULL)
						{
							ReadMCMol reader(*inmol, format, 1, 0, ReduceConfs);
							OBMol mol;

							while (reader.read(mol))
							{
								db.addMolToDatabase(mol, info.id, info.name);
							}
						}

						delete uncompressed_inmol;
						delete inmol;
					}
				}

				liginfos = vector<LigandInfo>(); //free up memory before building indices

				if(!NoIndex)
				{
				  db.createSpatialIndex(); //will write stats
				}
			} //destructors called
			if(d != nd-1)
				exit(0);
		}
		else //parent
		{
			numrunning++;
			if(numrunning >= maxt) //using all threads
			{
				int status;
				wait(&status);
				numrunning--;
			}
		}
	}

	int status;
	while (wait(&status) > 0)
	{
		if (!WIFEXITED(status) && WEXITSTATUS(status) != 0)
			abort();
		continue;
	}

	//all done, create symlinks to non-timestamped directories
	assert(symlinks.size() == topdirectories.size());
	for (unsigned d = 0, nd = topdirectories.size(); d < nd; d++)
	{
		if(filesystem::exists(symlinks[d]))
		{
			//remove preexisting symlink
			if(filesystem::is_symlink(symlinks[d]))
			{
				filesystem::remove(symlinks[d]);
			}
			else
			{
				cerr << "Trying to replace a non-symlink: " << symlinks[d] << "\n";
				exit(-1);
			}
		}
		filesystem::create_directory_symlink(topdirectories[d], symlinks[d]);
		filesystem::path dbipath = topdirectories[d] / "dbinfo.json";
		ofstream dbfile(dbipath.c_str());
		Json::StyledStreamWriter jwrite;
		jwrite.write(dbfile, root);
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

	StripedSearchers databases;
	unsigned totalC = 0;
	unsigned totalM = 0;

	vector<filesystem::path> dbpaths;
	for(unsigned i = 0, n = Database.size(); i < n; i++)
	{
		dbpaths.push_back(filesystem::path(Database[i]));
	}
	loadDatabases(dbpaths, databases);

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

		PharmerQuery query(databases.stripes, qfile,
				filesystem::extension(inputFiles[i]), params,
				NThreads * databases.stripes.size());

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
	isotab.Init(); //avoid race conditions! stupid openbabel
	etab.Init();
	ttab.Init();
	resdat.Init();

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
	else if (Cmd == "phogrify")
	{
		handle_phogrify_cmd(pharmas);
	}
	else if(Cmd == "showpharma")
	{
		pharmas.write(cout);
	}
	else if (Cmd == "dbcreate")
	{
		handle_dbcreate_cmd(pharmas);
	}
	else if(Cmd == "dbcreateserverdir")
	{
		handle_dbcreateserverdir_cmd(pharmas);
	}
	else if (Cmd == "dbsearch")
	{
		handle_dbsearch_cmd();
	}
	else if(Cmd == "fixsmina")
	{
	  handle_fixsmina_cmd();
	}
	else if (Cmd == "server")
	{
		vector<filesystem::path> prefixpaths;
		unordered_map<string, StripedSearchers > databases;

                //total hack time - fcgi uses select which can't
		//deal with file descriptors higher than 1024, so let's reserve some
		#define MAXRESERVEDFD (SERVERTHREADS*2)
		int reservedFD[MAXRESERVEDFD] = {0,};
		for(unsigned i = 0; i < MAXRESERVEDFD; i++)
		{
				reservedFD[i] = open("/dev/null",O_RDONLY);
		}
		//loadDatabases will open a whole bunch of files
		if(Prefixes.length() > 0 && Database.size() > 0)
		{
			cerr << "Cannot specify both dbdir and prefixes\n";
			exit(-1);
		}
		else if(Database.size() > 0)
		{
			//only one subset
			vector<filesystem::path> dbpaths;
			for(unsigned i = 0, n = Database.size(); i < n; i++)
			{
				dbpaths.push_back(filesystem::path(Database[i]));
			}
			loadDatabases(dbpaths, databases[""]);
		}
		else
		{
			//use prefixes
			ifstream prefixes(Prefixes.c_str());
			string line;
			while(getline(prefixes, line))
			{
				if(filesystem::exists(line))
				{
					prefixpaths.push_back(filesystem::path(line));
				}
				else
					cerr << line << " does not exist\n";
			}
			if(prefixpaths.size() == 0)
			{
				cerr << "No valid prefixes\n";
				exit(-1);
			}
			loadFromPrefixes(prefixpaths, databases);
		}
		//now free reserved fds
		for(unsigned i = 0; i < MAXRESERVEDFD; i++)
		{
				close(reservedFD[i]);
		}
		pharmer_server(Port, prefixpaths, databases, LogDir, MinServer, MinPort);
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

