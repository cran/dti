#include "Fiber.h"

using namespace std;

Fiber::Fiber()
{
	start = NULL;
	end   = NULL;
	length = 0;
}

Fiber::Fiber(Voxel* start)
{
	this->start = start;
	this->end   = start;
	length = 1;
}

//Fiber::~Fiber()
//{
//	delete start;
//	delete end;
//	
//	fclose(voxels);
//	
//	delete voxels;
//}

int Fiber::getLength()		{return length;}

void Fiber::add_at_end(Voxel &add)
{
	if (end == NULL)
	{
		add.setPrev(NULL);
		add.setNext(NULL);
		start = &add;
		end   = &add;
	}
	else
	{
		add.setPrev(end);
		end->setNext(&add);
		end = &add;
		end->setNext(NULL);
	}
	
	length++;
}

void Fiber::add_at_start(Voxel &add)
{
	if (start == NULL)
	{
		add.setPrev(NULL);
		add.setNext(NULL);
		start = &add;
		end   = &add;
	}
	else
	{
		start->setPrev(&add);
		add.setNext(start);
		start = &add;
		start->setPrev(NULL);
	}
	
	length++;
}

void Fiber::del_at_start()
{
	start->getNext()->setPrev(NULL);
	start->setNext(NULL);
	
	delete start;
	
	length--;
}

void Fiber::unvisit()
{
	if (length == 0)
	{
		return;
	}
	
	Voxel* aktuell = start;
	
	int i;
	
	for (i = 0; i < length; i++)
	{
		aktuell->setVisited(false);
		aktuell = aktuell->getNext();
		
//		if (aktuell == NULL)
//		{
//			break;
//		}	
	}
}

void Fiber::print()
{
	if (length == 0)
	{
		return;
	}
	
	Voxel* aktuell = start;
	
	int i;
	
	for (i = 0; i < length; i++)
	{
		aktuell->print();
		aktuell = aktuell->getNext();
		
//		if (aktuell == NULL)
//		{
//			break;
//		}	
	}
}
