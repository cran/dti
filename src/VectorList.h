#ifndef VECTORLIST_H_
#define VECTORLIST_H_

#include "Vector.h"

class VectorList{
	
	private:
			
		Vector* start;
		Vector* end;
		
		int length;
		int minLength;
		int num_nan;
		
	public:
		/**  constructors & destructor  **/
		VectorList();
		VectorList(Vector&);
		
//		~VectorList();
		
		/** getters  **/
		int getLength();
		int getMinLength();
		int getNum_Nan();
		Vector& getStart();
		Vector& getEnd();
		
		/**  LinkedList-methods  **/
		void add_list(VectorList&);
		
		void add_at_end(Vector&);
		void add_at_start(Vector&);
		void del_at_start();
		
		//  print-method
		void print();
		void print(int until);
};

#endif /*VECTORLIST_H_*/
