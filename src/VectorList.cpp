#include "VectorList.h"
#include <R.h>

using namespace std;

VectorList::VectorList()
{
	this->start = NULL;
	this->end   = NULL;
	length = 0;
	this->minLength = 6;
	num_nan = 0;
}

VectorList::VectorList(Vector &set_start)
{
	Vector *temp = new Vector(set_start);
	
	this->start = temp;
	this->end   = temp;
	this->start->setPrev(NULL);
	this->  end->setNext(NULL);
	length = 1;
	this->minLength = 6;
	
	if (isnan(set_start.getComponents()[1]))
	{
		num_nan = 1;
	}
}

//VectorList::~VectorList()
//{
//	while (this.length > 0)
//	{
//		this.del_at_start();
//	}
//}

int VectorList::getLength()			{return length;}
int VectorList::getMinLength()		{return minLength;}
int VectorList::getNum_Nan()		{return num_nan;}
Vector& VectorList::getStart()		{return *start;}
Vector& VectorList::getEnd()		{return *end;}

void VectorList::add_list(VectorList &addL)
{
	if (addL.getLength() >= 2*minLength)
	{
		Vector *nan = new Vector(0., R_NaN, 0.);
		addL.add_at_start(*nan);
		
		int i;
		for (i = 0; i < addL.getLength();)
		{
			this->add_at_end(addL.getStart());
			addL.del_at_start();
		}
	}
	else
	{
//		Rprintf("Fiber to short and will be ignored.\n");
	}
}

void VectorList::add_at_end(Vector &add)
{
	Vector *temp = new Vector(add.getComponents()[0], add.getComponents()[1], add.getComponents()[2]);
	
//	Rprintf("\n################\n#  add_at_end  #\n################\n");
//	
//	Rprintf("prev<-temp:\n");
//	temp->getPrev()->print();
//	Rprintf("temp:\n");
//	temp->print();
//	Rprintf("temp->next:\n");
//	temp->getNext()->print();
//	
//	Rprintf("=====================\n");
//	
//	Rprintf("prev<-end:\n");
//	end->getPrev()->print();
//	Rprintf("end:\n");
//	end->print();
	
	if (end == NULL)
	{
		start = temp;
		end   = temp;
		temp->setPrev(NULL);
		temp->setNext(NULL);
	}
	else
	{
		temp->setPrev(end);
		end->setNext(temp);
		end = temp;
		end->setNext(NULL);
	}
	
	if (isnan(add.getComponents()[1]))
	{	
		num_nan++;
	}
	
	length++;
}

void VectorList::add_at_start(Vector &add)
{
	Vector *temp = new Vector(add.getComponents()[0], add.getComponents()[1], add.getComponents()[2]);
	
	if (start == NULL)
	{
		start = temp;
		end   = temp;
		temp->setPrev(NULL);
		temp->setNext(NULL);
	}
	else
	{
		start->setPrev(temp);
		temp->setNext(start);
		start = temp;
		start->setPrev(NULL);
	}
	
	if (isnan(add.getComponents()[1]))
	{
		num_nan++;
	}
	
	length++;
}

void VectorList::del_at_start()
{
	if (isnan(start->getComponents()[1]))
	{
		num_nan--;
	}
	
	if (length > 1)
	{
		start = start->getNext();
		start->getPrev()->setNext(NULL);
		start->setPrev(NULL);
	}
	else
	{
// 		Rprintf("%f, %f, %f\n", start->getComponents()[0],
// 		start->getComponents()[1], start->getComponents()[2]);
		delete start;
	}
	
	length--;
}

void VectorList::print()
{
	if (length == 0)
	{
		return;
	}

//	Rprintf("========== Fiber start ==========\n");
	
	Vector* aktuell = start;
	
	while (aktuell != NULL)
	{
		aktuell->print();
		aktuell = aktuell->getNext();
	}

//	Rprintf("=========== Fiber  end ===========\n");	
}

void VectorList::print(int until)
{
	int counter = 0;
	
	if (length == 0)
	{
		return;
	}
	
	Vector* aktuell = start;
	
	while (aktuell != NULL && counter < until)
	{
		aktuell->print();
		aktuell = aktuell->getNext();
		counter++;
	}
}
