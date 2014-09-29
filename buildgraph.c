#include <stdio.h>
#include <stdlib.h>

int main(int argc,char* argv[])
{
	if(argc < 2)
	{
		printf("Command line arguments required: An input .rvr file, and an output .gra file.\n");
		return 1;
	}

	unsigned int i,j,id,numparents,N;
	int offset;
	FILE* input = fopen(argv[1],"r");
	FILE* output = fopen(argv[2],"w");

	//Check if files opened
	if(input == NULL)
	{
		printf("Error opening file %s.\n",argv[1]);
		return 1;
	}
	if(output == NULL)
	{
		printf("Error opening file %s.\n",argv[2]);
		return 1;
	}

	fscanf(input,"%u",&N);

	for(i=0;i<N;i++)
	{
		fscanf(input,"%u\n%u",&id,&numparents);
		if(i == 0)	offset = -id + 1;

		for(j=0;j<numparents;j++)
		{
			fscanf(input,"%u",&id);
			fprintf(output,"%u ",id + offset);
		}
		fprintf(output,"\n");
	}

	return 0;
}

