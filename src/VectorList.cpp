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

VectorList::VectorList(Vector &start)
{
	Vector *temp = new Vector(start);
	
	this->start = temp;
	this->end   = temp;
	this->start->setPrev(NULL);
	this->  end->setNext(NULL);
	length = 1;
	this->minLength = 6;
	
	if (isnan(start.getComponents()[1]))
	{
		num_nan = 1;
	}
}

//VectorList::~VectorList()
//{
//	delete start;
//	delete end;
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
//		printf("Fiber to short and will be ignored.\n");
	}
}

void VectorList::add_at_end(Vector &add)
{
	Vector *temp = new Vector(add.getComponents()[0], add.getComponents()[1], add.getComponents()[2]);
	
//	printf("\n################\n#  add_at_end  #\n################\n");
//	
//	printf("prev<-temp:\n");
//	temp->getPrev()->print();
//	printf("temp:\n");
//	temp->print();
//	printf("temp->next:\n");
//	temp->getNext()->print();
//	
//	printf("=====================\n");
//	
//	printf("prev<-end:\n");
//	end->getPrev()->print();
//	printf("end:\n");
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

//	printf("========== Fiber start ==========\n");
	
	Vector* aktuell = start;
	
	while (aktuell != NULL)
	{
		aktuell->print();
		aktuell = aktuell->getNext();
	}

//	printf("=========== Fiber  end ===========\n");	
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
