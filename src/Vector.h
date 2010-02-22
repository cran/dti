#ifndef VECTOR_H_
#define VECTOR_H_
#include <iostream>
#include <string>
#include <sstream>
#include <math.h>

using namespace std;

/*
 * The class 'Vector' offers typical vector operations.
 */

class Vector
{
	private:
		// length of Vector 
		int n;
		
		// components of Vector
		double* components;
		
		// LinkedList-Variables
		Vector* next;
		Vector* prev;

	public:
		/**  constructors & destructor  **/
		Vector() {};
		Vector(int);
		Vector(double, double, double);
		Vector(double*, int);
		
//		~Vector();
		
		// print-method
		void print();
		
		/**  getters  **/
		double* getComponents();
		int 	getN();
		Vector* getNext();
		Vector* getPrev();
		
		/**  setters  **/
		void setNext(Vector*);
		void setPrev(Vector*);
		
		/**  mathematical functions  **/
		
		// norm of a vector
		static double norm(Vector&);
		
		//  divide vector with scalar
		Vector& operator/(double);
		
		// multiply vector with scalar
		Vector& operator*(double);
		
		// multiply vector with vector
		double operator*(Vector&);
		
		// add vector to vector
		Vector& operator+(Vector&);

		// subtract vector from vector
		Vector& operator-(Vector&);

		// cross product of two vectors
		Vector& cross(Vector&);
};

#endif /*VECTOR_H_*/
