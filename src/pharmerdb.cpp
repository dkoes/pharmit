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
 * pharmerdb.cpp
 *
 *  Created on: Aug 2, 2010
 *      Author: dkoes
 */

#include "tripletmatching.h"
#include <boost/lexical_cast.hpp>
#include <fcntl.h>
#include <openbabel/mol.h>
#include <openbabel/descriptor.h>
#include "CommandLine2/CommandLine.h"
#include "Triplet.h"
#include "BitSetTree.h"
#include "ReadMCMol.h"
#include "Timer.h"
#include "pharmerdb.h"
#include "TripleIndexer.h"
#include "PMol.h"
#include "MolProperties.h"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "SminaConverter.h"
#include "ShapeObj.h"
#include "ShapeResults.h"
#include "ShapeConstraints.h"
#include "pharminfo.h"

using namespace boost;
using namespace std;
using namespace OpenBabel;

extern cl::opt<bool> Quiet;
extern cl::opt<unsigned> ReduceConfs;
extern cl::opt<bool> ComputeThresholds;

//location comparison functions for pointdata
bool comparePointDataX(const ThreePointData& lhs, const ThreePointData& rhs)
{
	return lhs.l1 < rhs.l1;
}

//location comparison functions for pointdata
bool comparePointDataY(const ThreePointData& lhs, const ThreePointData& rhs)
{
	return lhs.l2 < rhs.l2;
}

//location comparison functions for pointdata
bool comparePointDataZ(const ThreePointData& lhs, const ThreePointData& rhs)
{
	return lhs.l3 < rhs.l3;
}

//write out to buf, return number of bytes written
unsigned MolDataCreator::ConfCreator::write(unsigned char *buf)
{
	unsigned pos = 0;
	memcpy(buf, &header, sizeof(header));
	pos += sizeof(header);

	memcpy(buf + pos, molData.data(), molData.size());
	pos += molData.size();

	return pos;
}

//initialize data for a conformation
MolDataCreator::ConfCreator::ConfCreator(unsigned mid, float mw, unsigned cid,
		OBMol& m) :
		header(mid, mw, cid)
{
	PMolCreator mol(m, true);
	stringstream mdata;
	if (mol.writeBinary(mdata))
		molData = mdata.str();
	else
		molData = "";
}

unsigned MolDataCreator::maxIndex = 0;

//generate pharma points and database infor for mol
void MolDataCreator::processMol(OBMol& mol, MolProperties& props, unsigned mid)
{
	//identify pharma
	mcpoints.size();
	getPharmaPointsMC(pharmas, mol, mcpoints);


	if (mcpoints.size() > 0)
		props.setHB(mcpoints[0]); //hb feature cnts are not conformer specific

	double weight = mol.GetMolWt();
	//count number of rotatable bounds
	unsigned numrot = countRotatableBonds(mol);

	//store each conformation seperately, this results in about 50% more space used,
	//but as much as 50% performance improvement when retrieving mols

	//generate information for each conformation
	for (unsigned c = 0, nc = mol.NumConformers(); c < nc; c++)
	{
		unsigned confidx = c;
		mol.SetConformer(confidx);
		confs.push_back(ConfCreator(mid, weight, confidx, mol));
		ConfCreator& conf = confs.back();

		if (!conf.isValid())
		{
			cerr << mol.GetTitle() << " is too large.  Omitting.\n";
			confs.pop_back();
			continue;
		}
		//setup pharma points
		const vector<PharmaPoint>& points = mcpoints[confidx];

		if (maxIndex < points.size())
			maxIndex = points.size();

		//compute and canonicalize all triple points
		unsigned n = points.size();
		if (n > (1 << TPD_INDEX_BITS))
		{
			cerr << "Too many pharma points (" << n << ") in "
					<< mol.GetTitle() << "\n";
			n = (1 << TPD_INDEX_BITS);
		}

		for (unsigned i = 0; i < n; i++)
		{
			for (unsigned j = i + 1; j < n; j++)
			{
				for (unsigned k = j + 1; k < n; k++)
				{
					unsigned pos = tindex(points[i].pharma->index,
							points[j].pharma->index, points[k].pharma->index);
					assert(pos < pdatas.size());
					ThreePointData trip(bufferSize, weight, numrot, points, i,
							j, k);
					pdatas[pos].push_back(trip);
				}
			}
		}

		confOffsets.push_back(bufferSize);
		bufferSize += conf.byteSize();
		numConfs++;
	}
}

//flatten all the conf info into a buffer
void MolDataCreator::createBuffer()
{
	buffer = new unsigned char[bufferSize];
	unsigned pos = 0;
	for (unsigned i = 0, n = confs.size(); i < n; i++)
	{
		pos += confs[i].write(&buffer[pos]);
	}
	assert(pos == bufferSize);

	//clear mem
	confs.clear();
}

//write out the data with proper offsets to files (assumed positioned properly and single threaded)
unsigned MolDataCreator::write(FILE *molData,
		vector<PointDataFile>& pointDataFiles)
{
	unsigned long location = ftell(molData);

	location = ftell(molData);

	unsigned mid = location >> (TPD_MOLDATA_BITS - TPD_MOLID_BITS);

	//update the point datas

	for (unsigned p = 0, np = pdatas.size(); p < np; p++)
	{
		for (unsigned i = 0, n = pdatas[p].size(); i < n; i++)
		{
			pdatas[p][i].molPos += location;
		}
	}

	//write out
	fwrite(buffer, bufferSize, 1, molData);
	for (unsigned p = 0, np = pdatas.size(); p < np; p++)
	{
		pointDataFiles[p].write(&pdatas[p][0], sizeof(pdatas[p][0]),
				pdatas[p].size());
	}

	//align the start of every mol so the upper MOLID_BITS correspond to an mol id
	//(assuming the conformers fit within MOLDATA-MOLID bits)
	//we assume (hope) the filesystem has support for holes
	location = ftell(molData);
	location >>= (TPD_MOLDATA_BITS - TPD_MOLID_BITS);
	location++;
	location <<= (TPD_MOLDATA_BITS - TPD_MOLID_BITS);
	fseek(molData, location, SEEK_SET);

	if (ftell(molData) >= (1L << TPD_MOLDATA_BITS))
	{
		cerr
				<< "Input database is too large (moldata more than 1TB). Split and retry.\n";
		abort();
	}
	return mid;
}

void MolData::clear()
{
	mol = NULL;
	confidx = 0;
	mid = 0;
	molWeight = 0;
}

bool MolData::read(FILE *molData, PMolReader& areader)
{
	clear();
	MolDataHeader header;
	//get fixed sized header
	int ret = fread(&header, sizeof(header), 1, molData);
	assert(ret == 1);

	//copy to self
	confidx = header.confidx;
	molWeight = header.molWeight;
	mid = header.molID;

	//now get mol itself
	mol = areader.readPMol(molData);

	return mol != NULL;
}

//read in mol data at position location from file molData
bool MolData::read(FILE *molData, unsigned long location, PMolReader& areader)
{
	fseek(molData, location, SEEK_SET);
	return read(molData, areader);
}

//read in mol data at position location in array molData
bool MolData::read(unsigned char *molData, unsigned long location,
		PMolReader& areader)
{
	clear();
	MolDataHeader header;
	//get fixed sized header
	unsigned long loc = location;
	memcpy(&header, &molData[loc], sizeof(header));
	loc += sizeof(header);

	//copy to self
	confidx = header.confidx;
	molWeight = header.molWeight;
	mid = header.molID;

	//now get mol itself
	mol = areader.readPMol(&molData[loc]);

	return mol != NULL;
}

//read just metadata, don't parse full mol
void MolData::readDataOnly(unsigned char *molData, unsigned long location)
{
	clear();
	MolDataHeader header;
	//get fixed sized header
	unsigned long loc = location;
	memcpy(&header, &molData[loc], sizeof(header));
	loc += sizeof(header);

	//copy to self
	confidx = header.confidx;
	molWeight = header.molWeight;
	mid = header.molID;
}

void PharmerDatabaseCreator::initializeDatabases()
{
	//initialize self-managed flat file databases
	//open for writing
	//info
	filesystem::path ipath = dbpath;
	ipath /= "info";
	assert(!filesystem::exists(ipath));
	info = fopen(ipath.string().c_str(), "w+"); //create
	assert(info);

	//human readable pharmacophore classes
	filesystem::path lpath = dbpath;
	lpath /= "lookup";
	lookup = fopen(lpath.string().c_str(), "w"); //create
	assert(lookup);

	//mid index (to better support lots of conformers)
	filesystem::path mpath = dbpath;
	mpath /= "mids";
	midList = fopen(mpath.string().c_str(), "w"); //create
	assert(midList);

	//moldata
	filesystem::path mdpath = dbpath;
	mdpath /= "molData";
	molData = fopen(mdpath.string().c_str(), "w+");
	assert(molData);

	//smina data
	filesystem::path smipath = dbpath / "sminaIndex";
	sminaIndex = fopen(smipath.string().c_str(), "w+");
	assert(sminaIndex);
	filesystem::path smdpath = dbpath / "sminaData";
	sminaData = fopen(smdpath.string().c_str(), "w+");
	assert(sminaData);


	filesystem::path pipath = dbpath / "pharmInfo";
	pharmInfoData = fopen(pipath.string().c_str(), "w+");
	assert(pharmInfoData);


	//bincnts
	filesystem::path binpath = dbpath;
	binpath /= "binCnts";
	binData = fopen(binpath.string().c_str(), "w+");
	assert(binData);

	//output pharmas
	filesystem::path phpath = dbpath / "pharmas";
	ofstream pharmout(phpath.string().c_str());
	pharmas.write(pharmout);

	//pointData
	pointDataFiles.resize(tindex.size());
	for (int i = 0, n = tindex.size(); i < n; i++)
	{
		string pname = string("pointData_") + lexical_cast<string>(i);
		filesystem::path pdpath = dbpath / pname;
		pointDataFiles[i] = PointDataFile(pdpath.string(), i);
	}

	//geoData
	geoDataFiles.resize(tindex.size(), NULL);

	//property files
	MolProperties::createFiles(dbpath, propFiles);

	//initialize tempory files (used for out-of-memory partitioning)
	for (unsigned i = 0; i < NUMTMPFILES; i++)
	{
		filesystem::path tmppath = dbpath
				/ (lexical_cast<string>(i) + string(".tmp"));
		tmpFiles[i] = fopen(tmppath.c_str(), "w+");
		unlink(tmppath.c_str()); //basically anonymous disk-backed storage
	}
}

//write out index into "correct" mid
void PharmerDatabaseCreator::writeMIDs()
{
	if (midList)
	{
		fseek(midList, 0, SEEK_SET);
		unsigned cnt = 0;
		for (unsigned i = 0, n = mids.size(); i < n; i++)
		{
			fwrite(&mids[i], sizeof(unsigned), 1, midList);
			cnt++;
			while (i + 1 < n && cnt < mids[i + 1])
			{
				fwrite(&mids[i], sizeof(unsigned), 1, midList);
				cnt++;
			}
		}
	}
}

void PharmerDatabaseCreator::writeStats()
{
	if (info)
	{
		fseek(info, 0, SEEK_SET);
		fwrite(&stats, sizeof(stats), 1, info);
		fflush(info);
	}
	if (lookup)
	{
		tindex.dump(lookup);
	}
	writeMIDs();

	//annotate db json and output
	dbinfo["numMols"] = numMolecules();
	dbinfo["numConfs"] = numConformations();

	posix_time::ptime time(posix_time::second_clock::local_time());
	dbinfo["updated"] = boost::posix_time::to_simple_string(time);

	filesystem::path dbipath = dbpath / "dbinfo.json";
	ofstream dbfile(dbipath.c_str());
	Json::StyledStreamWriter jwrite;
	jwrite.write(dbfile, dbinfo);
}

//add a multiconformer mol to the database
void PharmerDatabaseCreator::addMolToDatabase(OBMol& mol, long uniqueid,
		const string& name)
{
	static OBAromaticTyper aromatics;
	static OBAtomTyper atyper;

	if (mol.NumAtoms() == 0) //skipped
		return;

	if (ReduceConfs > 0)
	{
		while (mol.NumConformers() > (int) ReduceConfs)
		{
			mol.DeleteConformer(mol.NumConformers() - 1);
		}
	}
	unsigned nc = mol.NumConformers();

	mol.AddHydrogens();

	aromatics.AssignAromaticFlags(mol);
	mol.FindSSSR();
	atyper.AssignTypes(mol);
	atyper.AssignHyb(mol);

	//calculate properties
	MolProperties props;
	props.calculate(mol, uniqueid);

	//store conformers aligned to inertial moments
	for(unsigned i = 0; i < nc; i++)
	{
		mol.SetConformer(i);
		ShapeObj::normalizeMol(mol);
	}

	if (mol.NumConformers() != (int) nc)
	{
		//there was a bug in openbabel that blew away conformers when adding hydrogens,
		//make sure we're using afixed version
		abort();
	}

	mol.SetTitle(name.c_str());
	//generate moldata
	MolDataCreator mdc(pharmas, tindex, mol, props, stats[NumMols]);

	unsigned mid = mdc.write(molData, pointDataFiles);
	mids.push_back(mid);
	props.write(mid, propFiles);

	//output smina and shape data here
	SminaConverter::MCMolConverter mcsmina(mol);
	const vector<unsigned>& confOffsets = mdc.ConfOffsets();
	unsigned long mloc = mid;
	mloc <<= (TPD_MOLDATA_BITS - TPD_MOLID_BITS);
	ShapeObj::MolInfo minfo(mol, mloc);

	for (unsigned i = 0, n = confOffsets.size(); i < n; i++)
	{
		stringstream data;
		mcsmina.convertConformer(i, data);
		unsigned sz = data.str().size();
		unsigned long pos = mloc + confOffsets[i]; //location in molData
		unsigned long smpos = ftell(sminaData); //location in smina dta

		//output index mapping
		fwrite(&pos, sizeof(pos), 1, sminaIndex);
		fwrite(&smpos, sizeof(smpos), 1, sminaIndex);

		//output actual data
		fwrite(&sz, sizeof(sz), 1, sminaData);
		fwrite(data.str().c_str(), sizeof(char), sz, sminaData);

		//shape data
		minfo.molPos = pos;
		minfo.pharmPos = writePharmacophoreInfo(pharmInfoData, mdc.getConfFeatures(i),pharmas);

		mol.SetConformer(i);
		ShapeObj shobj(mol, minfo, shapedb.getDimension(), shapedb.getResolution());
		shapedb.addObject(shobj);
	}

	stats[NumMols]++;
	stats[NumConfs] += mdc.NumConfs();
	stats[NumDbPoints] += mdc.NumPoints();
}

//coalesce very small collections and
//mmap pointData files, closing file pointer
void PharmerDatabaseCreator::initPointDataArrays()
{
	assert(pointDataArrays == NULL);
	unsigned n = pointDataFiles.size();
	pointDataArrays = new MMappedRegion<ThreePointData> [n];
	for (unsigned i = 0; i < n; i++)
	{
		if (pointDataFiles[i].isValid())
		{
			pointDataFiles[i].close();
			pointDataArrays[i].map(pointDataFiles[i].name, false, true);
		}
	}
}

//split along the the longest access
//fill in split info
static void chooseSplit(const BoundingBox& box, SplitInfo& info,
		unsigned long long numPoints, unsigned long long numUnique,
		unsigned pos)
{

	if (numPoints == 0)
	{
		info.type = NoSplit;
	}
	else if ((numUnique <= MAX_UNIQUE_POINTS_PER_LEAF || numPoints
			< MAX_POINTDATAS_PER_LEAF))
	{
		//not worth further dividing
		info.type = NoSplit;
	}
	else if (2 * pos >= SPLITS_PER_GEOPAGE && numPoints
			< MIN_POINTDATAS_PER_PAGE)
	{
		info.type = NoSplit;
		//don't create a whole new geo page
	}
	else //do a split
	{
		int xdiff = box.maxx - box.minx;
		int ydiff = box.maxy - box.miny;
		int zdiff = box.maxz - box.minz;

		if (xdiff >= ydiff && xdiff >= zdiff)
		{
			info.type = SplitXFixed;
			info.min = box.minx;
			info.max = box.maxx;
			info.func = comparePointDataX;
		}
		else if (ydiff >= xdiff && ydiff >= zdiff)
		{
			info.type = SplitYFixed;
			info.min = box.miny;
			info.max = box.maxy;
			info.func = comparePointDataY;
		}
		else //z
		{
			info.type = SplitZFixed;
			info.min = box.minz;
			info.max = box.maxz;
			info.func = comparePointDataZ;
		}
	}
}

//find the median value in one sequential scan using counts
static unsigned short findMedianValue(ThreePointData *start,
		ThreePointData *end,
		const SplitInfo& info, unsigned short min, unsigned short max)
{
	unsigned long total = 0;
	int left = min;
	int right = left + (max - min);
	unsigned long counts[right - left + 1];
	memset(counts, 0, sizeof(counts));

	for (ThreePointData *itr = start; itr != end; itr++)
	{
		counts[info.splitVal(*itr) - left]++;
		total++;
	}
	unsigned long mid = total / 2;
	unsigned long cnt = 0;
	for (int i = left; i <= right; i++)
	{
		cnt += counts[i - left];
		if (cnt > mid)
			return i;
	}
	abort();
}

class PivotCompareRnd
{
	short pivot;
	const SplitInfo& info;
	public:
	PivotCompareRnd(short p, const SplitInfo& i) :
			pivot(p), info(i)
	{
	}

	//flip a coint as to which way equal vals go
	bool operator()(const ThreePointData& x)
	{
		int diff = info.splitVal(x) - pivot;
		if (diff < 0)
			return true;
		if (diff > 0)
			return false;
		return ::random() % 2;
	}
};

// in or out of memory partitioning
// divides up input (start to end) into three groups - less than, equal to, and greater than median
//then outputs them in order and returns a middle position
ThreePointData* PharmerDatabaseCreator::partitionData(ThreePointData *start,
		ThreePointData *end, SplitInfo& info, unsigned short median)
{
	unsigned long num = end - start;

	if (num < pdatasFitInMemory)
	{
		PivotCompareRnd rcmp(median, info);
		ThreePointData *ret = partition(start, end, rcmp);
		if (ret == end || ret == start) //we got unlucky and did not partition well
		{
			sort(start, end, info.func);
			ret = start + (end - start) / 2;
		}
		return ret;
	}
	else
	{
		assert(tmpFiles[0]); //less
		assert(tmpFiles[1]); //equal
		assert(tmpFiles[2]); //more

		FILE *less = tmpFiles[0];
		FILE *equal = tmpFiles[1];
		FILE *more = tmpFiles[2];

		rewind(less);
		rewind(equal);
		rewind(more);

		unsigned long nLess = 0, nEqual = 0, nMore = 0;

		for (ThreePointData *itr = start; itr != end; itr++)
		{
			unsigned short val = info.splitVal(*itr);
			if (val < median)
			{
				nLess++;
				fwrite(itr, sizeof(ThreePointData), 1, less);
			}
			else if (val == median)
			{
				nEqual++;
				fwrite(itr, sizeof(ThreePointData), 1, equal);
			}
			else
			{
				nMore++;
				fwrite(itr, sizeof(ThreePointData), 1, more);
			}
		}

		//copy back in order
		rewind(less);
		size_t num = fread(start, sizeof(ThreePointData), nLess, less);
		assert(num == nLess);

		rewind(equal);
		num = fread(start + nLess, sizeof(ThreePointData), nEqual, equal);
		assert(num == nEqual);

		rewind(more);
		num = fread(start + nLess + nEqual, sizeof(ThreePointData), nMore,
				more);
		assert(num == nMore);

		return start + (end - start) / 2;
	}
}

//create an internal kd tree
//each node splits points data in two along the axis with the largest spread
//every node has an explicit bounding box (so in that sense this is more like an R-tree)
void PharmerDatabaseCreator::doSplitInPage(unsigned pharma, FILE *geoFile,
		GeoKDPage& page,
		unsigned pos, ThreePointData *start, ThreePointData *end,
		ThreePointData *begin,
		unsigned depth)
{
	//first compute the bounding box and count points
	unsigned long numPoints = end - start;
	unsigned long numUnique = 0;
	TripleInt last(0, 0, 0);
	BoundingBox box;
	for (ThreePointData *itr = start; itr != end; itr++)
	{
		TripleInt coords(itr->l1, itr->l2, itr->l3);
		box.update(coords);
		if (coords != last)
		{
			last = coords;
			numUnique++;
		}
	}

	//choose splitting dimension
	SplitInfo info;
	chooseSplit(box, info, numPoints, numUnique, pos);

	page.nodes[pos].box = box;
	page.nodes[pos].splitType = info.type;
	if (info.type == NoSplit)
	{
		//this is the end
	}
	else //need to split
	{
		//choose a median value, points with this value can be on either side
		unsigned short medianVal = findMedianValue(start, end, info,
				info.getMin(box),
				info.getMax(box));
		ThreePointData *median = partitionData(start, end, info, medianVal);

		page.nodes[pos].splitVal = medianVal;
		page.nodes[pos].splitData = median - begin;

		if (2 * pos + 1 < SPLITS_PER_GEOPAGE) //stay on internal page
		{
			doSplitInPage(pharma, geoFile, page, 2 * pos, start, median, begin,
					depth);
			doSplitInPage(pharma, geoFile, page, 2 * pos + 1, median, end,
					begin,
					depth);
		}
		else //children pushed to new page
		{
			assert(2*pos >= SPLITS_PER_GEOPAGE);
			unsigned lpos = 2 * pos - SPLITS_PER_GEOPAGE;
			unsigned rpos = lpos + 1;
			assert(rpos < SPLITS_PER_GEOPAGE);
			page.nextPages[lpos] = doSplitNewPage(pharma, geoFile, start,
					median,
					begin, depth);
			page.nextPages[rpos] = doSplitNewPage(pharma, geoFile, median, end,
					begin,
					depth);
		}
	}
}

//create a geokdpage that indexes the data between start and end in geoFile
//writes out page to geoFile and returns its location as an index into an array
//of geokdpages
unsigned long PharmerDatabaseCreator::doSplitNewPage(unsigned pharma,
		FILE *geoFile,
		ThreePointData *start, ThreePointData *end, ThreePointData *begin,
		unsigned depth)
{
	GeoKDPage page;

	//reserve space
	unique_lock<shared_mutex> lock(fileAccessLock);
	fseek(geoFile, 0, SEEK_END);
	unsigned long location = ftell(geoFile);
	assert(location % sizeof(GeoKDPage) == 0);
	int ret = fwrite(&page, sizeof(page), 1, geoFile);
	assert(ret == 1);
	lock.unlock();

	depth++;
	//create the page
	doSplitInPage(pharma, geoFile, page, 1, start, end, begin, depth);

	stats[NumInternalPages]++;

	//write it out
	lock.lock();
	fseek(geoFile, location, SEEK_SET);
	location /= sizeof(GeoKDPage); //indexify
	ret = fwrite(&page, sizeof(page), 1, geoFile);
	assert(ret == 1);

	return location;
}

//saturated binning
static unsigned lengthbin(unsigned len)
{
	unsigned ret = len / LENGTHDIV;
	if (ret >= LENGTH_BINS)
		ret = LENGTH_BINS - 1;

	return ret;
}

void PharmerDatabaseCreator::incrementBinCnt(const ThreePointData& t,
		unsigned pclass)
{
	binnedCnts[pclass][lengthbin(t.l1)][lengthbin(t.l2)][lengthbin(t.l3)]++;
}

//create spatial index for all ijk triples
void PharmerDatabaseCreator::createIJKSpatialIndex(int p)
{
	//create spatial
	ThreePointData *start = pointDataArrays[p].begin();
	ThreePointData *begin = start;
	ThreePointData *end = pointDataArrays[p].end();

	if (start == end)
		return;

	//count binned frequencies
	for (ThreePointData *itr = start; itr != end; itr++)
	{
		incrementBinCnt(*itr, p);
	}

	//open
	string gname = string("geoData_") + lexical_cast<string>(p);
	filesystem::path gpath = dbpath / gname;
	FILE * geoFile = fopen(gpath.string().c_str(), "w");
	geoDataFiles[p] = geoFile;
	assert(geoFile);

	//write out blank page, there's nothing at index 0
	fseek(geoFile, 0, SEEK_SET);
	GeoKDPage blank;
	fwrite(&blank, sizeof(blank), 1, geoFile);

	//top-down recursively create kd/r tree spatial data structure
	doSplitNewPage(p, geoFile, start, end, begin, 0);
}

/* Create spatial index. */
void PharmerDatabaseCreator::createSpatialIndex()
{
	Timer t;
	cout << "Creating spatial index..." << endl;
	//close and then mmap pointinfo file
	initPointDataArrays();

	for (unsigned p = 0, n = tindex.size(); p < n; p++)
		createIJKSpatialIndex(p);

	for (unsigned i = 0, n = binnedCnts.size(); i < n; i++)
		fwrite(binnedCnts[i].c_array(), LENGTH_BINS * LENGTH_BINS * LENGTH_BINS,
				sizeof(unsigned), binData);

	shapedb.createIndex();

	cout << stats[NumConfs] << "\tconformations\n";
	cout << stats[NumMols] << "\tmolecules\n";
	cout << stats[NumInternalPages] << "\tinternal pages\n";
	cout << "Index creation time: " << t.elapsed() << "s ("
			<< t.elapsedProcess() << "s)" << "\n";

	writeStats();

}

//Searcher

void PharmerDatabaseSearcher::initializeDatabases()
{
	//initialize self-managed flat file databases
	//open for reading
	//info - and read it in
	filesystem::path ipath = dbpath;
	ipath /= "info";
	info = fopen(ipath.string().c_str(), "r");
	assert(info);
	int read = fread(&stats, sizeof(stats), 1, info);
	if (read != 1)
	{
		cerr << ipath.string() << "\n";
		abort();
	}

	//json info
	filesystem::path dbipath = dbpath / "dbinfo.json";
	ifstream json(dbipath.c_str());
	Json::Reader reader;
	if (!reader.parse(json, dbinfo))
	{
		cerr << "Error reading database info JSON\n";
		exit(-1);
	}

	//pharmas
	filesystem::path pharmpath = dbpath / "pharmas";
	ifstream pharmain(pharmpath.string().c_str());
	if (!pharmas.read(pharmain))
	{
		cerr << "Invalid pharma specification in database.\n";
		return;
	}

	tindex.set(pharmas.size());

	//moldata
	filesystem::path mdpath = dbpath;
	mdpath /= "molData";
	molData.map(mdpath.string(), true, true);

	//smina index and data - do not need to exist
	filesystem::path smIndex = dbpath / "sminaIndex";
	if (filesystem::exists(smIndex))
	{
		sminaIndex.map(smIndex.string(), true, true);
		filesystem::path smData = dbpath / "sminaData";
		sminaData.map(smData.string(), true, true);
	}

	//pharmacophore info
	filesystem::path phInfo = dbpath / "pharmInfo";
	if (filesystem::exists(phInfo))
	{
		pharmInfoData.map(phInfo.string(), true, true);
	}


	//mids
	filesystem::path mpath = dbpath;
	mpath /= "mids";
	if (filesystem::exists(mpath)) //back-wards compat
		midList.map(mpath.string(), true, true, true);

	//length histogram
	filesystem::path binpath = dbpath / "binCnts";
	binnedCnts.map(binpath.string(), true, false);
	if (binnedCnts.length() == 0)
	{
		cerr << "Missing binnedCnts " << binpath << "\n";
		return;
	}
	//pointData
	unsigned n = tindex.size();
	tripletDataArrays = new MMappedRegion<ThreePointData> [n];
	for (unsigned i = 0; i < n; i++)
	{
		string pname = string("pointData_") + lexical_cast<string>(i);
		filesystem::path pdpath = dbpath / pname;
		if (filesystem::exists(pdpath))
			tripletDataArrays[i].map(pdpath.string(), true, true);
	}

	//geoData
	geoDataArrays = new MMappedRegion<GeoKDPage> [n];
	for (unsigned i = 0; i < n; i++)
	{
		string gname = string("geoData_") + lexical_cast<string>(i);
		filesystem::path gpath = dbpath / gname;
		if (filesystem::exists(gpath))
			geoDataArrays[i].map(gpath.string(), true, true);
	}

	//property data
	MolProperties::initializeReader(dbpath, props);
	valid = true;

	filesystem::path shape = dbpath / "shape";
	shapesearch.load(shape);
}

unsigned PharmerDatabaseSearcher::getBinCnt(unsigned pclass, unsigned i,
		unsigned j, unsigned k)
{
	unsigned index = pclass;
	index *= LENGTH_BINS;
	index += i;
	index *= LENGTH_BINS;
	index += j;
	index *= LENGTH_BINS;
	index += k;

	return binnedCnts[index];
}

//generate ranking and return index of best triplet
//use binnedCnts histogram
unsigned PharmerDatabaseSearcher::rankTriplets(
		const vector<QueryTriplet>& triplets, vector<double>& ranking)
{
	double minval = HUGE_VAL;
	unsigned besttrip = 0;
	ranking.resize(triplets.size(), 0);

	//add up values from histogram
	for (unsigned i = 0, n = triplets.size(); i < n; i++)
	{
		double val = 0;
		const QueryTriplet& trip = triplets[i];
		unsigned pclass = tindex(trip.getPharma(0), trip.getPharma(1),
				trip.getPharma(2));

		for (unsigned x = lengthbin(trip.getRanges()[0].min), mx =
				lengthbin(trip.getRanges()[0].max); x <= mx; x++)
			for (unsigned y = lengthbin(trip.getRanges()[1].min), my =
					lengthbin(trip.getRanges()[1].max); y <= my; y++)
				for (unsigned z = lengthbin(trip.getRanges()[2].min),
						mz = lengthbin(trip.getRanges()[2].max); z
						<= mz; z++)
				{
					val += getBinCnt(pclass, x, y, z);
				}

		ranking[i] = val * triplets[i].numOrderings(); //symmetries will increases matches

		if (val < minval)
		{
			minval = val;
			besttrip = i;
		}
	}
	return besttrip;
}

//get values along appropriate axis for split
static void getMinMax(const Triplet& t, SplitType split, int& min, int& max)
{
	min = max = 0;
	switch (split)
	{
		case SplitXFixed:
			min = t.getRanges()[0].min;
			max = t.getRanges()[0].max;
			break;
		case SplitYFixed:
			min = t.getRanges()[1].min;
			max = t.getRanges()[1].max;
			break;
		case SplitZFixed:
			min = t.getRanges()[2].min;
			max = t.getRanges()[2].max;
			break;
		default:
			;
	}
}

inline void dump(ostream& out, const ThreePointData& tpd)
{
	out << tpd.molPos << ": " << tpd.l1 << " " << tpd.l2 << " " << tpd.l3
			<< "\n";
}

void PharmerDatabaseSearcher::queryProcessPoints(QueryInfo& t,
		unsigned long startLoc, unsigned long endLoc)
{
	const ThreePointData *start = t.data + startLoc;
	const ThreePointData *end = t.data + endLoc;

	//put the range of point datas into the priorty queue,
	//let the PQ filter out bad values
	unsigned cnt = 0;
	for (const ThreePointData *itr = start; itr != end; itr++)
	{
		unsigned mid = getBaseMID(itr->molID());
		if (t.M.add(mid, *itr, t.triplet, t.which))
		{
			cnt++;
		}
	}
	if (cnt == 0)
		emptyCnt++;
	matchedCnt += cnt;
	processCnt += (end - start);
}

cl::opt<unsigned> searchCutoff("search-cutoff", cl::desc(
		"Number of triplets to stop kd-tree splitting at"), cl::Hidden,
		cl::init(2048));

static unsigned pageCnt = 0;
void PharmerDatabaseSearcher::queryIndex(QueryInfo& t, const GeoKDPage *page,
		unsigned pos, unsigned long startLoc, unsigned long endLoc)
{
	if (startLoc == endLoc)
		return;
	if (pos >= SPLITS_PER_GEOPAGE)
	{
		if (t.stopEarly)
			return;
		//need another page
		unsigned long index = page->nextPages[pos - SPLITS_PER_GEOPAGE];
		pageCnt++;
		queryIndex(t, &t.pages[index], 1, startLoc, endLoc);
	}
	else
	{
		//internal to this page
		const GeoKDPageNode& node = page->nodes[pos];
		//are there potentially points within this region?
		if (!t.triplet.inRange(node.box))
		{
			return;
		}
		if (node.splitType == NoSplit)
		{
			queryProcessPoints(t, startLoc, endLoc);
		}
		else if (endLoc - startLoc < searchCutoff)
		{
			//don't spend a lot of time at the edges, just read and process
			queryProcessPoints(t, startLoc, endLoc);
		}
		else if (t.triplet.allInRange(node.box))
		{
			//this worth doing both when accessing disk and, surprisingly,
			//when accessing memory
			//all points beneath are valid (maybe)
			queryProcessPoints(t, startLoc, endLoc);
		}
		else //consider splitter
		{
			int min, max;
			getMinMax(t.triplet, node.splitType, min, max);

			//some sanity checks
			assert(startLoc <= node.splitData);
			assert(endLoc >= node.splitData);

			//always do left and then right so we access data sequentially
			//if min or max are exactly equal, then we must consider both sides of the split
			//in case the split is in the middle of a range of equal values
			if ((max > node.splitVal && min < node.splitVal)
					|| max == node.splitVal
					|| min == node.splitVal)
			{
				//have to do both sides
				//first go left
				queryIndex(t, page, 2 * pos, startLoc, node.splitData);
				queryIndex(t, page, 2 * pos + 1, node.splitData, endLoc);
			}
			else if (max <= node.splitVal) //just the left
			{
				queryIndex(t, page, 2 * pos, startLoc, node.splitData);
			}
			else if (min >= node.splitVal) //just right
			{
				queryIndex(t, page, 2 * pos + 1, node.splitData, endLoc);
			}
		}
	}
}


void PharmerDatabaseSearcher::generateShapeMatches(const ShapeConstraints& constraints,
		ShapeResults& results)
{
	//create tree version of constraints
	GSSTreeSearcher::ObjectTree small = shared_ptr<const MappableOctTree>(
					MappableOctTree::createFromGrid(constraints.getInclusiveGrid()), free);
	shared_ptr<MappableOctTree> big = shared_ptr<MappableOctTree>(
					MappableOctTree::createFromGrid(constraints.getExclusiveGrid()), free);
	GSSTreeSearcher::ObjectTree lig = shared_ptr<const MappableOctTree>(
					MappableOctTree::createFromGrid(constraints.getLigandGrid()), free);

	big->invert();

	shapesearch.dc_search(small, big, lig, true, results);
}

//put all matching triplets into Q
//each triplet is expanded to get all overlapping trips
//this is intentionally not multi-threaded - performance-wise it's better to
//split the database up and parralel match
void PharmerDatabaseSearcher::generateTripletMatches(const vector<vector<
		QueryTriplet> >& triplets, TripletMatches& M, bool& stopEarly)
{
	pageCnt = 0;
	emptyCnt = 0;
	processCnt = 0;
	matchedCnt = 0;
	for (unsigned i = 0, n = triplets.size(); i < n; i++)
	{
		for (unsigned t = 0, nt = triplets[i].size(); t < nt; t++)
		{
			const QueryTriplet& trip = triplets[i][t];
			unsigned pclass = tindex(trip.getPharma(0), trip.getPharma(1),
					trip.getPharma(2));
			const GeoKDPage *pages = geoDataArrays[pclass].begin();
			const ThreePointData *data = tripletDataArrays[pclass].begin();
			QueryInfo qinfo(trip, t, pages, data, i, M, stopEarly);
			queryIndex(qinfo, &pages[1], 1, 0,
					tripletDataArrays[pclass].length());
		}

		if (!Quiet)
		{
			cout << i << " MatchedCnt " << matchedCnt << " ProcessedCnt "
					<< processCnt << "  " << 100.0 * double(matchedCnt)
							/ processCnt << "% ";
			cout << "Pages " << pageCnt << " ";
			cout << "Empty " << emptyCnt << " ";
			cout << "Chunks: " << M.numChunks() << "\n";
		}
		if (!M.nextIndex())
			break;

		if (stopEarly)
			break;
	}
	if (!Quiet)
	{
		cout << "PageCnt :" << pageCnt << "\n";
		cout << "EmptyCnt : " << emptyCnt << "\n";
		cout << "MatchedCnt : " << matchedCnt << "\n";
		cout << "ProcessedCnt : " << processCnt << "\n";
	}
}

typedef pair<unsigned long, unsigned long> ulong_pair;
static bool first_pair_cmp(const ulong_pair& lhs, const ulong_pair& rhs)
{
	return lhs.first < rhs.first;
}
//find the smina data location given the moldata location and copy the smina
//data to out
void PharmerDatabaseSearcher::getSminaData(unsigned long molloc, ostream& out)
{
	ulong_pair findit(molloc, 0);

	//molloc is the conformer location in molData, do binary search
	ulong_pair *pos = lower_bound(sminaIndex.begin(), sminaIndex.end(), findit,
			first_pair_cmp);

	if (pos == sminaIndex.end())
	{
		pos = sminaIndex.end() - 1;
	}
	else if (pos->first > molloc)
	{
		pos--;
	}
	assert(pos >= sminaIndex.begin());

	//pos is now the correct spot
	unsigned long sminaloc = pos->second;
	unsigned sz = 0;
	memcpy(&sz, sminaData.begin() + sminaloc, sizeof(unsigned)); //size of smina data
	out.write(sminaData.begin() + sminaloc + sizeof(unsigned), sz);
}

//given a location in pharmInfoData and a query, check to see if the pharmacophore at phlocation
//matches the query (which is assumed to be aligned to already match the pharmacophore)
bool PharmerDatabaseSearcher::alignedPharmasMatch(unsigned long phlocation, const vector<PharmaPoint>& query)
{
	if(pharmInfoData.length() == 0) return true; //gracefully fail
	if(phlocation >= pharmInfoData.length()) abort(); //not graceful

	return pharmacophoreMatchesQuery(pharmInfoData.begin()+phlocation, query, pharmas);

}


