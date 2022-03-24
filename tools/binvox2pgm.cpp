/* This program converts a 3D discrete object from the binvox format
 * to PGM (3D) format.
 *
 * This code is adapted from
 *   https://www.patrickmin.com/binvox/read_binvox.cc
 *
 * Compile with:
 *   g++ binvox2pgm.cpp -O3 -o binvox2pgm
 */

#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>

typedef unsigned char byte;

static int version;
static int depth, height, width;
static int size;
static byte *voxels = 0;
static float tx, ty, tz;
static float scale;


int read_binvox(std::string filespec)
{
	std::ifstream *input = new std::ifstream(filespec.c_str(), std::ios::in | std::ios::binary);
	//
	// read header
	//
	std::string line;
	*input >> line; // #binvox
	if (line.compare("#binvox") != 0) {
		std::cerr << "Error: first line reads [" << line << "] instead of [#binvox]" << std::endl;
		delete input;
		return 0;
	}
	*input >> version;
	std::cout << "reading binvox version " << version << std::endl;

	depth = -1;
	int done = 0;
	while (input->good() && !done) {
		*input >> line;
		if (line.compare("data") == 0) done = 1;
		else if (line.compare("dim") == 0) {
			*input >> depth >> height >> width;
		}
		else if (line.compare("translate") == 0) {
			*input >> tx >> ty >> tz;
		}
		else if (line.compare("scale") == 0) {
			*input >> scale;
		}
		else {
			std::cerr << " unrecognized keyword [" << line << "], skipping" << std::endl;
			char c;
			do { // skip until end of line
				c = input->get();
			} while(input->good() && (c != '\n'));
		}
	}
	if (!done) {
		std::cerr << " error reading header" << std::endl;
		return 0;
	}
	if (depth == -1) {
		std::cout << " missing dimensions in header" << std::endl;
		return 0;
	}
	size = width * height * depth;
	voxels = new byte[size];
	if (!voxels) {
		std::cerr << " error allocating memory" << std::endl;
		return 0;
	}

	//
	// read voxel data
	//
	byte value;
	byte count;
	int index = 0;
	int end_index = 0;
	int nr_voxels = 0;

	input->unsetf(std::ios::skipws); // need to read every byte now (!)
	*input >> value; // read the linefeed char
	while ((end_index < size) && input->good()) {
		*input >> value >> count;
		if (input->good()) {
			end_index = index + count;
			if (end_index > size)
				return 0;
			for (int i=index; i < end_index; i++)
				voxels[i] = value;
			if (value)
				nr_voxels += count;
			index = end_index;
		} // if file still ok
	} // while
	input->close();
	std::cout << " " << nr_voxels << " voxels read" << std::endl;
	return 1;
}


int main(int argc, char **argv)
{
	if (argc != 3) {
		std::cerr << "Usage: binvox2pgm <binvox filename> <pgm filename> " << std::endl << std::endl;
		exit(1);
	}
	if (!read_binvox(argv[1])) {
		std::cerr << "Error reading [" << argv[1] << "]" << std::endl << std::endl;
		exit(1);
	}
	//
	// now write the data to as ASCII
	//
	std::ofstream *out = new std::ofstream(argv[2]);
	if(!out->good()) {
		std::cerr << "Error opening [" << argv[2] << "]" << std::endl << std::endl;
		exit(1);
	}
	std::cout << "Writing voxel data to ASCII file... ";

	*out << "P2" << std::endl;
	*out << width+2 << " " << height+2 << " " << depth+2 << std::endl;
	*out << "255" << std::endl;
	int index = 0;
	for (int k = 0; k < depth+2; k++)
	{
		for (int j = 0; j < height+2; j++)
		{
			for (int i = 0; i < width+2; i++)
			{
				if (i==0 || i==width+1 || j==0 || j==height+1 || k==0 || k==depth+1)
					*out << "0 ";
				else
				{
					*out << (char) (voxels[index] + '0') << " ";
					index++;
				}
			}
			*out << std::endl;
		}
	}
	out->close();
	std::cout << "Done" << std::endl;
	return 0;
}
