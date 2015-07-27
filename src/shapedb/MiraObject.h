/*
 * MiraObject.h
 *
 * A MiraObject is created from
 * Given a DIRECTORY containing *.mira files will iterate through all the files
 * and view them as objects suitable for GSSTree creation.
 *
 *  Created on: May 31, 2013
 *      Author: dkoes
 */

#ifndef MIRAOBJECT_H_
#define MIRAOBJECT_H_

#include "Cube.h"
#include "MGrid.h"
#include <boost/filesystem.hpp>
#include <boost/multi_array.hpp>
#include <cstring>
#include <arpa/inet.h>

class MiraIterator;

struct MiraHeader
{
	char fileid[5]; // must be VOXEL
	unsigned char control_z; // must be 26 (0x1a)
	unsigned short version; // must be 1 (?)
	unsigned short xres;
	unsigned short yres;
	unsigned short zres;
	unsigned short flag; // 0x4 if RGB data, 0x8 if RGBA data, 0 = byte data
	unsigned int map_offset; // must be 256
	unsigned int voxel_offset; // equal to 256 + (xres + yres + zres) * sizeof(double)
	char unused[104]; //we put resolution = num here
	char text[128]; // information text
};

class MiraObject {
	MiraHeader header;
	boost::multi_array<bool, 3> data;
	unsigned maxdim; //require cube
	double resolution; //mira objects are resolution-less, so may to set manually
						//our mira objects have the resolution in the unused field
public:
	typedef MiraIterator iterator;

	MiraObject(): maxdim(0), resolution(1)
	{
	}

	~MiraObject()
	{
	}

	//get/set the resolution (scaling factor)
	void setResolution(double r) { resolution = r; }
	double getResolution() const { return resolution; }

	unsigned getDimension() const { return maxdim; }

	//return the number of "on" voxels
	unsigned numSetBits() const {
		unsigned cnt = 0;
		for(unsigned i = 0; i < maxdim; i++)
			for(unsigned j = 0; j < maxdim; j++)
				for(unsigned k = 0; k < maxdim; k++)
					cnt += data[i][j][k];
		return cnt;
	}

	//NOTE: since I'm too lazy to implement a hierarchy of grids for fast
	//intersection testing, I'm just going to return true since the tree
	//building code will always eventually check with containsPoint
	//Should MiraObjects actually matter for something other than getting
	//this paper published in a CS journal, this should be improved
	bool intersects(const Cube& cube) const
	{
		return true;
	}

	//return true if point is within object
	//grids are unit resolution anchored at the origin, but gss trees
	//are centered at origin so translate
	bool containsPoint(float x, float y, float z) const
	{
		x /= resolution;
		y /= resolution;
		z /= resolution;
		unsigned X = round(-0.5+x+maxdim/2.0);
		unsigned Y = round(-0.5+y+maxdim/2.0);
		unsigned Z = round(-0.5+z+maxdim/2.0);

		if(X >= maxdim || Y >= maxdim || Z >= maxdim)
			return false;

		return data[X][Y][Z];
	}

	//just writes out whatever is in the text field
	void write(ostream& out) const
	{
		out.write(header.text, strlen(header.text)+1);
	}

	//read contents of mira stream into object, replacing existing object
	//sets text field to filename if specified
	bool read(istream& in, std::string fname = "")
	{
		in.read((char*)&header,sizeof(header));

		if(header.fileid[0] != 'V') return false; //simplest check
		if(header.flag != 0) {
			std::cerr << "Error: only support byte data MIRA format.\n";
			return false;
		}
		unsigned maxx = ntohs(header.xres);
		unsigned maxy = ntohs(header.yres);
		unsigned maxz = ntohs(header.zres);

		if(maxx != maxy || maxy != maxz)
		{
			std::cerr << "Error: only support cubic volumes.\n";
			return false;
		}
		maxdim = maxx;
		header.text[127] = 0; //ensure null termination
		if(fname.size() > 0)
			strncpy(header.text, fname.c_str(), sizeof(header.text)-1);

		const char* resprefix = "resolution = ";
		unsigned prefixn = strlen(resprefix);
		if(strncmp(header.unused,resprefix,prefixn) == 0)
		{
			resolution = atof(&header.unused[prefixn]);
		}

		//absorb map
		unsigned voxeloff = header.voxel_offset;
		unsigned mapoff = header.map_offset;
		while(voxeloff > mapoff)
		{
			in.get();
			voxeloff--;
		}
		//only resize if necessary
		if(data.strides()[0] != maxdim || data.strides()[1] != maxdim || data.strides()[2] != maxdim)
		{
			data.resize(boost::extents[maxdim][maxdim][maxdim]);
		}
		//zero out
		std::fill(data.origin(), data.origin()+data.size(), 0);

		for(unsigned i = 0, n = maxdim*maxdim*maxdim; i < n; i++)
		{
			unsigned c = in.get();
			data.origin()[i] = c;
		}
		return true;
	}
};

class MiraIterator
{
	vector<string> files;
	unsigned file_pos;
	MiraObject current;
public:
	virtual ~MiraIterator() {}


	MiraIterator(const string& dname): file_pos(0)
	{
		//get all the mira files form the directory
		using namespace boost::filesystem;
		for(directory_iterator ditr = directory_iterator(dname), end; ditr != end; ditr++)
		{
			path name = *ditr;
			if(name.extension() == ".mira")
			{
				files.push_back(name.string());
			}
		}

		if(files.size() > 0) {
			ifstream in(files[0].c_str());
			current.read(in, files[0]);
		}
	}

	//validity check
	operator bool() const
	{
		return file_pos < files.size();
	}

	//current mol
	const MiraObject& operator*() const
	{
		return current;
	}

	void operator++()
	{
		file_pos++;
		if(file_pos < files.size())
		{
			ifstream in(files[file_pos].c_str());
			current.read(in, files[file_pos]);
		}
	}

};



#endif /* MIRAOBJECT_H_ */
