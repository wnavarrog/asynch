#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "structs.h"

int main(int argc,char* argv[])
{
	unsigned int i,j,size;
	unsigned int N = pow(4,atoi(argv[1])-1);
	unsigned int done = 1;
	short unsigned int split_size = 3;
	double tol = .001;

	//Initialize sys
	Link** sys = malloc(N*sizeof(Link*));
	for(i=0;i<N;i++)
	{
		sys[i] = malloc(sizeof(Link));
		sys[i]->location = i;
		sys[i]->numparents = 0;
		sys[i]->parents = malloc(split_size*sizeof(Link*));
		sys[i]->c = NULL;
	}

	//Split each link into three links repeatedly
	while(done < N)
	{
		size = done;
		for(i=0;i<size;i++)
		{
			if(sys[i]->numparents == 0)
				sys[i]->numparents = split_size;
			else
			{
				sys[done]->numparents = split_size;
				for(j=0;j<split_size;j++)
					sys[done]->parents[j] = sys[i]->parents[j];
			}

			for(j=0;j<split_size;j++)
			{
				sys[i]->parents[j] = sys[done];
				sys[done]->c = sys[i];
				done++;
			}
		}
	}

	//Create .rvr file
	FILE* file = fopen("peano.rvr","w");
	fprintf(file,"%u\n\n",N);

	for(i=0;i<N;i++)
	{
		fprintf(file,"%u\n%u ",sys[i]->location,sys[i]->numparents);
		for(j=0;j<sys[i]->numparents;j++)
			fprintf(file,"%u ",sys[i]->parents[j]->location);
		fprintf(file,"\n\n");
	}

	//Create .ini file
	file = fopen("peano.ini","w");
	fprintf(file,"0\n%u\n0.0\n\n",N);
	
	for(i=0;i<N;i++)
		fprintf(file,"%u\n%f\n\n",i,1.0);

	//Create .prm file
	file = fopen("peano.prm","w");
	fprintf(file,"%u\n\n",N);

	for(i=0;i<N;i++)
		fprintf(file,"%u\n%f %f %f %f %f %f %f %f %f %f %f %f\n\n",i,.500,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0);
/*
	//Create .rkd file
	file = fopen("peano.rkd","w");
	fprintf(file,"0\n%u\n\n",N);

	for(i=0;i<N;i++)
	{
		fprintf(file,"%u\n%f\n0.0\n%f\n0.0\n",i,tol,tol);
		if(sys[i]->numparents == 0)	fprintf(file,"2\n.1\n\n");
		else				fprintf(file,"1\n.1\n\n");
	}
*/
	//Cleanup
	fclose(file);
	for(i=0;i<N;i++)
	{
		free(sys[i]->parents);
		free(sys[i]);
	}
	free(sys);

	return 0;
}

