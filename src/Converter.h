#ifndef CONVERTER_H_
#define CONVERTER_H_

#include "Voxel.h"
#include <R.h>
#include <Rinternals.h>

class Converter
{
	private:
		Voxel* voxels;
		
	public:
		Converter(double*, double*, int*, int, int, int);
		Converter(double*, double*, int*, int*,int, int, int, int);
		Voxel& getVoxels();
};

#endif /*CONVERTER_H_*/
