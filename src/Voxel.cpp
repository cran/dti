#include "Voxel.h"

using namespace std;

Voxel::Voxel()
{
	this->order = 0;
	this->dir_index = 0;
	this->position = Vector(-1, -1, -1);
	this->directions = new Vector[order];
	this->x = -1;
	this->y = -1;
	this->z = -1;
	
	this->anisotropy = 0.;
	
	this->visited = false;
	this->startable = false;
	
	this->next = NULL;
	this->prev = NULL;
}

Voxel::Voxel(int x_in, int y_in, int z_in, int order_in, Vector& directions_in, double anisotropy_in)
{
	this->order = order_in;
	this->dir_index = 0;
	this->directions = &directions_in;
	this->anisotropy = anisotropy_in;
	this->x = x_in;
	this->y = y_in;
	this->z = z_in;

	this->position = Vector((double)this->x,(double)this->y,(double)this->z);
	
	visited = false;
	startable = false;
	
	next = NULL;
	prev = NULL;
}

//Voxel::Voxel(int x, int y, int z, int order, Vector* directions, double anisotropy)
//{
//	this->order = order;
//	this->directions = directions;
//	this->anisotropy = anisotropy;
//	this->x = x;
//	this->y = y;
//	this->z = z;
//
//	this->position = Vector((double)this->x,(double)this->y,(double)this->z);
//	
//	visited = false;
//	startable = true;
//	
//	next = NULL;
//	prev = NULL;
//}

//Voxel::~Voxel()
//{
//	delete next;
//	delete prev;
//	
//	 position.~Vector();
//	direction.~Vector();
//}

/**  getters  **/

bool Voxel::isVisited()			{return visited;}
bool Voxel::isStartable()		{return startable;}
double Voxel::getAnisotropy()	{return anisotropy;}
Vector* Voxel::getDirections()	{return directions;}
Vector& Voxel::getPosition()	{return position;}
int Voxel::getX()				{return x;}
int Voxel::getY()				{return y;}
int Voxel::getZ()				{return z;}
int Voxel::getOrder()			{return order;}
int Voxel::getDir_Index()		{return dir_index;}
Voxel* Voxel::getNext()			{return next;}
Voxel* Voxel::getPrev()			{return prev;}

/**  setters  **/
void Voxel::setDir_Index(int set_dir_index)	{this->dir_index = set_dir_index;}
void Voxel::setVisited(bool set_visited)	{this->visited = set_visited;}
void Voxel::setStartable(bool set_startable){this->startable = set_startable;}
void Voxel::setNext(Voxel *set_next)		{this->next = set_next;}
void Voxel::setPrev(Voxel *set_prev)		{this->prev = set_prev;}

/**  print-methode  **/
void Voxel::print()
{
  //    return;
  //  Rprintf("%d %d %d, visited = %d, aniso = %f\n", x, y, z, this->visited, this->anisotropy);
}
