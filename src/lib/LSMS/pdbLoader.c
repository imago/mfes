#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "molTypes.h"

Molecule loadCATrace(char *fileName)
{
	FILE *protFile;
	int i = 0;
	char tmp[9];
	char add;

	Molecule p = (Molecule)malloc(sizeof(struct prot));
	char *line = (char *)malloc(sizeof(char)*85);

	protFile = fopen(fileName,"r");
	if (!protFile)
	{
		printf("Cannot open input molecule!\n");
		exit(1);
	}


	p->npoints = 0;
	while(fgets(line,90,protFile)!=NULL)
	{
		if (strncmp(line,"ATOM",4)==0)
		{
			switch (line[13])
			{
				case	'C' : if (line[14]=='A') p->npoints++; break;
				case	'c' : if (line[14]=='a') p->npoints++; break;
			}
		}
	}
	fclose(protFile);

	protFile = fopen(fileName,"r");

	p->xpoints = (float *)malloc(p->npoints*sizeof(float));
	p->ypoints = (float *)malloc(p->npoints*sizeof(float));
	p->zpoints = (float *)malloc(p->npoints*sizeof(float));
	p->rpoints = (float *)malloc(p->npoints*sizeof(float));


	while(fgets(line,90,protFile)!=NULL)
	{
		if (strncmp(line,"ATOM",4)==0)
		{
			add = 0;
			switch (line[13])
			{
				case	'C' : if (line[14]=='A') { add=1; p->rpoints[i] = 1.872f; } break;
				case	'c' : if (line[14]=='a') { add=1; p->rpoints[i] = 1.872f; } break;
			}
			if (add)
			{
				strncpy(tmp,line+30,8);
				tmp[8]='\0';
				p->xpoints[i] = (float)atof(tmp);
				strncpy(tmp,line+38,8);
				tmp[8]='\0';
				p->ypoints[i] = (float)atof(tmp);
				strncpy(tmp,line+46,8);
				tmp[8]='\0';
				p->zpoints[i] = (float)atof(tmp);
				i++;
			}
		}
	}

	fclose(protFile);
	free(line);
	return p;
}

Molecule loadMolecule(char *fileName)
{
	FILE *protFile;
	int i = 0;
	char tmp[9];
	char add;


	Molecule p = (Molecule)malloc(sizeof(struct prot));
	char *line = (char *)malloc(sizeof(char)*85);

	protFile = fopen(fileName,"r");
	if (!protFile)
	{
		printf("Cannot open input molecule!\n");
		exit(1);
	}


	p->npoints = 0;
	while(fgets(line,90,protFile)!=NULL)
	{
		if (strncmp(line,"ATOM",4)==0)
		{
			switch (line[13])
			{
				case	'C' : p->npoints++; break;
				case	'c' : p->npoints++; break;
				case	'N' : p->npoints++; break;
				case	'n' : p->npoints++; break;
				case	'O' : p->npoints++; break;
				case	'o' : p->npoints++; break;
				case	'S' : p->npoints++; break;
				case	's' : p->npoints++; break;
				case	'P' : p->npoints++; break;
				case	'p' : p->npoints++; break;
			}
		}
	}
	fclose(protFile);

	protFile = fopen(fileName,"r");

	p->xpoints = (float *)malloc(p->npoints*sizeof(float));
	p->ypoints = (float *)malloc(p->npoints*sizeof(float));
	p->zpoints = (float *)malloc(p->npoints*sizeof(float));
	p->rpoints = (float *)malloc(p->npoints*sizeof(float));

	while(fgets(line,90,protFile)!=NULL)
	{
		if (strncmp(line,"ATOM",4)==0)
		{
			add = 0;
			switch (line[13])
			{
				case	'C' : add=1; p->rpoints[i] = 1.872f; break;
				case	'c' : add=1; p->rpoints[i] = 1.872f; break;
				case	'N' : add=1; p->rpoints[i] = 1.507f; break;
				case	'n' : add=1; p->rpoints[i] = 1.507f; break;
				case	'O' : add=1; p->rpoints[i] = 1.400f; break;
				case	'o' : add=1; p->rpoints[i] = 1.400f; break;
				case	'S' : add=1; p->rpoints[i] = 1.848f; break;
				case	's' : add=1; p->rpoints[i] = 1.848f; break;
				case	'P' : add=1; p->rpoints[i] = 1.880f; break;
				case	'p' : add=1; p->rpoints[i] = 1.880f; break;
			}
			/*
			switch (line[13])
			{
				case	'C' : add=1; p->rpoints[i] = 1.71; break;
				case	'c' : add=1; p->rpoints[i] = 1.71; break;
				case	'N' : add=1; p->rpoints[i] = 1.71; break;
				case	'n' : add=1; p->rpoints[i] = 1.71; break;
				case	'O' : add=1; p->rpoints[i] = 1.71; break;
				case	'o' : add=1; p->rpoints[i] = 1.71; break;
				case	'S' : add=1; p->rpoints[i] = 1.71; break;
				case	's' : add=1; p->rpoints[i] = 1.71; break;
				case	'P' : add=1; p->rpoints[i] = 1.71; break;
				case	'p' : add=1; p->rpoints[i] = 1.71; break;
			}*/
			if (add)
			{
				strncpy(tmp,line+30,8);
				tmp[8]='\0';
				p->xpoints[i] = (float)atof(tmp);
				strncpy(tmp,line+38,8);
				tmp[8]='\0';
				p->ypoints[i] = (float)atof(tmp);
				strncpy(tmp,line+46,8);
				tmp[8]='\0';
				p->zpoints[i] = (float)atof(tmp);
				i++;
			}
		}
	}

	fclose(protFile);

	return p;
}

