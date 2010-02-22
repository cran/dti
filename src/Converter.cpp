
#include "Converter.h"

Voxel& Converter::getVoxels()
{
	return *voxels;
}

Converter::Converter(double* data_dir_coords, double* data_FA_values, int dim_x, int dim_y, int dim_z)
{
//	printf("%d, %d, %d\n", dim_x, dim_y, dim_z);
	
	voxels = new Voxel[dim_x*dim_y*dim_z];
	
//	printf("converting doubles to Voxels...\n");

	Vector *dir;
    
	double dir_x, dir_y, dir_z, FA;
	int i = 0, j = 0, k = 0, num_vector = 0, counter = 0;
	
    for (k = 0; k < dim_z; k++)
    {
	    for (j = 0; j < dim_y; j++)
	    {
	    	for (i = 0; i < dim_x; i++)
		   	{
		   		dir_x = data_dir_coords[num_vector];num_vector++;
		   		dir_y = data_dir_coords[num_vector];num_vector++;
		   		dir_z = data_dir_coords[num_vector];				num_vector++;
		   		FA 	  =  data_FA_values[i + j * dim_x + k * dim_x * dim_y];
		   		
				dir = new Vector(dir_x, dir_y, dir_z);
				
				Voxel v(i, j, k, 1, *dir, FA);
				
				voxels[i + j * dim_x + k * dim_x * dim_y] = v;
				
				counter++;
		   	}
	    }
    }

//    printf("Converted %d Voxels.\n", counter);
}


Converter::Converter(double* data_dir_coords, double* data_FA_values, int* data_order, int maxorder, int dim_x, int dim_y, int dim_z)
{
//	printf("%d, %d, %d\n", dim_x, dim_y, dim_z);
	
	voxels = new Voxel[dim_x*dim_y*dim_z];
	
//	printf("converting doubles to Voxels...\n");

	double dir_x, dir_y, dir_z, FA;
	int i = 0, j = 0, k = 0, l = 0, num_vector = 0, counter = 0, order = 0;
	
    for (k = 0; k < dim_z; k++)
    {
	    for (j = 0; j < dim_y; j++)
	    {
	    	for (i = 0; i < dim_x; i++)
		   	{
		   		FA 	  = data_FA_values[i + j * dim_x + k * dim_x * dim_y];
		   		order = data_order[i + j * dim_x + k * dim_x * dim_y];
		   		
		   		Vector *dir = new Vector[order];
		   		
		   		for (l = 0; l < order; l++)
		   		{
			   		dir_x = data_dir_coords[num_vector];num_vector++;
			   		dir_y = data_dir_coords[num_vector];num_vector++;
			   		dir_z = data_dir_coords[num_vector];num_vector++;
			   		
			   		dir[l] = Vector(dir_x, dir_y, dir_z);
		   		}
				
		   		num_vector += (maxorder - order)*3;
		   				   		
				Voxel v(i, j, k, order, *dir, FA);
				
				voxels[i + j * dim_x + k * dim_x * dim_y] = v;
				
				counter++;
		   	}
	    }
    }

//    printf("Converted %d Voxels.\n", counter);
}
