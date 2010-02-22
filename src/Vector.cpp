
#include "Vector.h"

using namespace std;

Vector::Vector(int numbArgs)
{
	this->n = numbArgs;
	
	double *input = new double[numbArgs];
		
	for (int i = 0; i < this->n; i++)
	{
		input[i] = 0.;
	}
	
	this->components = input;
}

Vector::Vector(double x, double y, double z)
{
	double *input = new double[3];
	input[0] = x;
	input[1] = y;
	input[2] = z;
	
	this->components = input;
	this->n = 3;
}

Vector::Vector(double* d_input, int length)
{
	this->n = length;
	this->components = d_input;
}

//Vector::~Vector()
//{
//	delete components;
//	delete prev;
//	delete next;
//}

void Vector::print()
{
	if (this == NULL)
	{
//		printf("Vector is NULL\n");
		return;
	}
	
//	printf(" _   _ \n");	
		
	for (int i = 0; i < this->n; i++)
	{
		if (i == n-1)
		{
//			printf("|_%f_|\n", this->components[i]);
		}
		else
		{
//			printf("| %f |\n", this->components[i]);
		}
	}
}

double* Vector::getComponents()	{return components;}
int Vector::getN()				{return n;}

Vector* Vector::getNext()		{return next;}
Vector* Vector::getPrev()		{return prev;}

void Vector::setNext(Vector *next)		{this->next = next;}
void Vector::setPrev(Vector *prev)		{this->prev = prev;}

double Vector::norm(Vector& v)
{
	double norm = 0.;
		
	for (int i = 0; i < v.n; i++)
	{
		norm += (v.components[i]*v.components[i]);
	}
	
	return sqrt(norm);
}

Vector& Vector::operator /(double d)
{
	static Vector out(this->n);
	
	for (int i = 0; i < out.n; i++)
	{
		out.components[i] = (this->components[i] / d);
	}
	
	return out;
}

Vector& Vector::operator *(double d)
{
	static Vector out(this->n);
	
	for (int i = 0; i < out.n; i++)
	{
		out.components[i] = this->components[i] *d;
	}

	return out;
}

double Vector::operator *(Vector& v2)
{
	if (this->n != v2.n)
	{
		cout << "error in multiplication: vectors must have the same length!" << endl;
		return 0.;
	}
	else
	{
		double out = 0.;
		
		for (int i = 0; i < this->n; i++)
		{
			out += this->components[i] * v2.components[i];
		}
		
		return out;
	}
}

Vector& Vector::operator +(Vector& v2)
{
	if (this->n != v2.n)
	{
		cout << "error in addition: vectors must have the same length!" << endl;
		
		Vector * out = NULL;
		return *out;
	}
	else
	{
		static Vector out(this->n);
		
		for (int i = 0; i < out.n; i++)
		{
			out.components[i] = this->components[i] + v2.components[i];
		}
		
		return out;
	}
}


Vector& Vector::operator -(Vector& v2)
{
	if (this->n != v2.n)
	{
		cout << "error in addition: vectors must have the same length!" << endl;
		
		static Vector out;
		return out;
	}
	else
	{
		static Vector out(this->n);
		
		for (int i = 0; i < out.n; i++)
		{
			out.components[i] = this->components[i] - v2.components[i];
		}
		
		return out;
	}
}

Vector& Vector::cross(Vector& v2)
{
	if (this->n != v2.n)
	{
		cout << "error im cross product: vectors must have the same length!" << endl;
		
		static Vector out;
		return out;
	}
	else
	{
		static Vector out(this->n);
		
		for (int i = 0; i < out.n; i++)
		{
			int j = i+1;
			int k = i+2;
			
			if (k >= out.n)
			{
				k -= out.n;
			}
			
			if (j >= out.n)
			{
				j -= out.n;
			}
			
			out.components[i] = this->components[j]*v2.components[k] - this->components[k]*v2.components[j];
		}
		
		return out;
	}
}

