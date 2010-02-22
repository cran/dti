#ifndef FIBERTRACKING_H_
#define FIBERTRACKING_H_

#include "Fiber.h"
#include "VectorList.h"

using namespace std;

class Fibertracking
{
	private:
		
		// index of current voxel
		int cur_voxel_index;
		
		// index of the voxel where the last fiber started
		int last_start_voxel;
		
		// index of the last plane direction
		int last_plane_dir;
		
		// number of fibers
		int num_fibers;
		
		Fiber currentFiber;
		VectorList curVectorList;
		VectorList allVectors;
	
		// dimension of the voxel cube
		int dim_x, dim_y, dim_z;
		
		// array of the voxels
		Voxel *voxels;
		// voxelextentions
		double voxelext_x, voxelext_y, voxelext_z;
	
		// intersection angle of two vectors
		double intersec_angle;
		
		// minimum anisotropy
		double min_anisotropy;
		
		// maximum intersection angle
		double max_intersec_angle;
	
		// position vector of the start points or cutting points
		Vector start_o;
	
		// normal vectors from the planes E1 to E6 (inward-looking)
		Vector n_e1;
		Vector n_e2;
		Vector n_e3;
		Vector n_e4;
		Vector n_e5;
		Vector n_e6;
		
		// true := changing direction of next voxel 
		bool change_dir;
		
		// tracks the next voxel in forward direction
		void nextVoxel_forward();
		
		// tracks the next voxel in backward direction
		void nextVoxel_backward();
		
		// tracks the hole fiber in forward direction
		void trackFiber_forward();
		
		// tracks the hole fiber in backward direction
		void trackFiber_backward();
		
	public:
		/**  constructors & destructor  **/
		Fibertracking();
		Fibertracking(Voxel&, int , int, int, double, double, double, double, double);
		
//		~Fibertracking();
		
		void findAllFibers();
		void findMarkedFibers(int*);
		int getLength();
		double* convertToDouble();
};

#endif /*FIBERTRACKING_H_*/
