#ifndef VOXEL_H_
#define VOXEL_H_

#include "Vector.h"

using namespace std;

/*
 * The class 'Voxel' contains the direction vectors and the anisotropy
 * from the given R variables.
 */

class Voxel
{
	private:
	
		/* x: 	x-location
 		 * y: 	y-location
 		 * z: 	z-location
 		 * direction:	direction-vector
		 * anisotropy:	anisotropy 		*/
	
		int x, y, z, order, dir_index;
		
		Vector position;
		Vector *directions;
		double anisotropy;
	
		bool startable;
		
		// LinkedList-Variables
		Voxel* next;
		Voxel* prev;
		
		bool visited;
		
//		Voxel *subvoxels;
	
	public:
		/**  constructors & destructor  **/
		Voxel();
		Voxel(int, int, int, int, Vector&, double);
		
//		~Voxel();
		
		/**  getters  **/
		bool isVisited();
		bool isStartable();
		double getAnisotropy();
		Vector* getDirections();
		Vector& getPosition();
		int getX();
		int getY();
		int getZ();
		int getOrder();
		int getDir_Index();
		Voxel* getNext();
		Voxel* getPrev();
//		Voxel* getSubvoxels();
		
		/**  setters  **/
		void setDir_Index(int);
		void setVisited(bool);
		void setStartable(bool);
		void setNext(Voxel*);
		void setPrev(Voxel*);
		
		/**  print-methode  **/
		void print();
};

#endif /*VOXEL_H_*/
