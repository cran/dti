#ifndef FIBER_H_
#define FIBER_H_

#include "Voxel.h"

/*
 * The class 'Fiber' is a LinkedList of Voxel objects.
 */

class Fiber
{
	private:
		Voxel* start;
		Voxel* end;
	
		int length;
		
	public:
		/**  constructors & destructors  **/
		Fiber();
		Fiber(Voxel*);
		
//		~Fiber();
		
		/** getters  **/
		int getLength();
		
		/**  LinkedList-methods  **/
		void add_at_end(Voxel&);
		void add_at_start(Voxel&);
		void del_at_start();
		
		// sets all voxels within the fiber to visited = false
		void unvisit();
		
		//  print-method
		void print();
};

#endif /*FIBER_H_*/
