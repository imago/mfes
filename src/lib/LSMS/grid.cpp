#include <stdlib.h>
#include <stdio.h>
#include "gridTypes.h"
#include "molTypes.h"
#include "grid.h"

typedef struct nb
{
	short x;
	short y;
	short z;
	struct nb* next;
} nbNode;

lGrid createGrid(short i)
{
	short s;
	short m,n,k;
	lGrid g = (lGrid)malloc(sizeof(struct g));
	g->N = i;
	s = 512/i;
	g->stepSize = s;

	g->matrix = (GridPoint***)malloc(i*sizeof(GridPoint**));
	for (m=0;m<i;m++)
		g->matrix[m] = (GridPoint**)malloc(i*sizeof(GridPoint*));
	for (m=0;m<i;m++)
		for (n=0;n<i;n++)
			g->matrix[m][n] = (GridPoint*)malloc(i*sizeof(GridPoint));

	for(m=0;m<i;m++)
		for (n=0;n<i;n++)
			for (k=0;k<i;k++)
	{
		g->matrix[m][n][k].point.x = -256+m*s;
		g->matrix[m][n][k].point.y = -256+n*s;
		g->matrix[m][n][k].point.z = -256+k*s;
		g->matrix[m][n][k].phi = 1; // everything is outside surface
		g->matrix[m][n][k].from=-2;
		g->matrix[m][n][k].dist= 0;
	}
	return g;
}

void signDistanceGridMol(lGrid g, Molecule mol, double PR)
{
	int i;
	int a,b,c;
	int len = mol->npoints;
	int icx,icy,icz;
	int r;
	float s = 512/g->N;
	float dx,dy,dz,fr;
	
	for (i=0;i<len;i++){
		r = (int)((mol->rpoints[i]+PR)/s)+1;
		//printf("r = %d\n",r);
		icx = (mol->xpoints[i]+256)/s;
		icy = (mol->ypoints[i]+256)/s;
		icz = (mol->zpoints[i]+256)/s;
		//printf("s gx gy gz mx my mz %f %d %d %d %f %f %f\n",s,g->matrix[icx][icy][icz].point.x,g->matrix[icx][icy][icz].point.y,g->matrix[icx][icy][icz].point.z,mol->xpoints[i],mol->ypoints[i],mol->zpoints[i]);
		for (a=icx-r;a<=icx+r;a++)
		for (b=icy-r;b<=icy+r;b++)
		for (c=icz-r;c<=icz+r;c++){
			dx = g->matrix[a][b][c].point.x - mol->xpoints[i];
			dy = g->matrix[a][b][c].point.y - mol->ypoints[i];
			dz = g->matrix[a][b][c].point.z - mol->zpoints[i];
			fr = mol->rpoints[i]+PR;
			if (dx*dx + dy*dy + dz*dz <= fr*fr) 
				g->matrix[a][b][c].phi = 0;
		}
	}
}

int findProbesMol(lGrid g, float PR)
{
	int i,j,k;
	int a,b,c;
	int co = 0;
	float fx,fy,fz,fr;
	float dx,dy,dz;
	int icx,icy,icz,r;
	int l = g->N;
	float s = 512/g->N;
	
	co = 0;

	for (i=1;i<l-1;i++)
	for (j=1;j<l-1;j++)
	for (k=1;k<l-1;k++)
		if (g->matrix[i][j][k].phi==1 &&
		   (g->matrix[i][j][k-1].phi==0 ||
		    g->matrix[i][j][k+1].phi==0 ||
		    g->matrix[i][j-1][k].phi==0 ||
		    g->matrix[i][j+1][k].phi==0 ||
		    g->matrix[i+1][j][k].phi==0 ||
		    g->matrix[i-1][j][k].phi==0))
			co++;

	Molecule m = (Molecule)malloc(co*sizeof(struct prot));	
	m->npoints = co;
	m->xpoints = (float *)malloc(co*sizeof(float));
	m->ypoints = (float *)malloc(co*sizeof(float));
	m->zpoints = (float *)malloc(co*sizeof(float));
	m->rpoints = (float *)malloc(co*sizeof(float));
	co = 0;
	for (i=1;i<l-1;i++)
	for (j=1;j<l-1;j++)
	for (k=1;k<l-1;k++)
		if (g->matrix[i][j][k].phi==1 &&
		   (g->matrix[i][j][k-1].phi==0 ||
		    g->matrix[i][j][k+1].phi==0 ||
		    g->matrix[i][j-1][k].phi==0 ||
		    g->matrix[i][j+1][k].phi==0 ||
		    g->matrix[i+1][j][k].phi==0 ||
		    g->matrix[i-1][j][k].phi==0)){
			m->xpoints[co] = g->matrix[i][j][k].point.x;
			m->ypoints[co] = g->matrix[i][j][k].point.y;
			m->zpoints[co] = g->matrix[i][j][k].point.z;
			m->rpoints[co] = PR;
			co++;
		    }

	for (i=0;i<m->npoints;i++){
		fx = m->xpoints[i];
		fy = m->ypoints[i];
		fz = m->zpoints[i];
		fr = m->rpoints[i];
		r = (int)(fr/s)+1;
		icx = (fx+256)/s;
		icy = (fy+256)/s;
		icz = (fz+256)/s;
		for (a=icx-r;a<=icx+r;a++)
		for (b=icy-r;b<=icy+r;b++)
		for (c=icz-r;c<=icz+r;c++){	
			if (a>=g->N || b>=g->N || c>=g->N || a<0 || b<0 || c<0)
				continue;
			dx = g->matrix[a][b][c].point.x - fx;
			dy = g->matrix[a][b][c].point.y - fy;
			dz = g->matrix[a][b][c].point.z - fz;
			if (dx*dx + dy*dy + dz*dz <= fr*fr) 
				g->matrix[a][b][c].phi = 1;
		}
	}
	free(m);
	return co;
}

void shrink(lGrid g, double PR)
{
	int xp,xn,yp,yn,zp,zn;
	int i,j,k;
	int co = 0,p;
	float fr;
	float dx,dy,dz;
	int l = g->N;
	float s = 512/g->N;
	float temp_dist;
	PPoint* nb_head, *nb_orig, *nb_next_head,*nb_next_orig;
	
	PR=PR*PR+0.5*s;	//be careful
	co = 0;

	for (i=1;i<l-1;i++)
	for (j=1;j<l-1;j++)
	for (k=1;k<l-1;k++)
		if (g->matrix[i][j][k].phi==1 &&
		   (g->matrix[i][j][k-1].phi==0 ||
		    g->matrix[i][j][k+1].phi==0 ||
		    g->matrix[i][j-1][k].phi==0 ||
		    g->matrix[i][j+1][k].phi==0 ||
		    g->matrix[i+1][j][k].phi==0 ||
		    g->matrix[i-1][j][k].phi==0))
			co++;
	
	nb_head = (PPoint*)malloc(co*sizeof(PPoint));	//head line
	nb_orig = (PPoint*)malloc(co*sizeof(PPoint));	//where does this point come from
	
	co=0;
	for (i=1;i<l-1;i++)
	for (j=1;j<l-1;j++)
	for (k=1;k<l-1;k++){
		//should we start from first "1" or from last "0"
		if (g->matrix[i][j][k].phi==1){
			if(g->matrix[i][j][k-1].phi==0 ||
		    		g->matrix[i][j][k+1].phi==0 ||
		    		g->matrix[i][j-1][k].phi==0 ||
		    		g->matrix[i][j+1][k].phi==0 ||
		    		g->matrix[i+1][j][k].phi==0 ||
		   		g->matrix[i-1][j][k].phi==0){
					nb_head[co].x = i;
					nb_head[co].y = j;
					nb_head[co].z = k;
					nb_orig[co].x = i;
					nb_orig[co].y = j;
					nb_orig[co].z = k;
					co++;
			}
		}
	}
	
	float abs = 0.35*s;

	while (co!=0){
//printf("%d\n",co);
		nb_next_head = (PPoint*)malloc(co*6*sizeof(PPoint));	//next level (at most 6 times larger)
		nb_next_orig = (PPoint*)malloc(co*6*sizeof(PPoint));	//it is possible to be larger than nb_head

		p = 0;
		for (i=0;i<co;i++){
			xp=nb_head[i].x+1;
			xn=nb_head[i].x-1;
			yp=nb_head[i].y+1;
			yn=nb_head[i].y-1;
			zp=nb_head[i].z+1;
			zn=nb_head[i].z-1;
			if (g->matrix[xp][nb_head[i].y][nb_head[i].z].phi==0 || 
			   (g->matrix[xp][nb_head[i].y][nb_head[i].z].phi==1&&g->matrix[xp][nb_head[i].y][nb_head[i].z].from!=-1)){
				if(g->matrix[xp][nb_head[i].y][nb_head[i].z].phi==0)	fr=PR;
				else	fr=g->matrix[xp][nb_head[i].y][nb_head[i].z].dist;
				
				dx = g->matrix[xp][nb_head[i].y][nb_head[i].z].point.x - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.x;
				dy = g->matrix[xp][nb_head[i].y][nb_head[i].z].point.y - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.y;
				dz = g->matrix[xp][nb_head[i].y][nb_head[i].z].point.z - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.z;
				temp_dist=dx*dx+dy*dy+dz*dz;
			
				if(temp_dist<=fr && (fr-temp_dist >  abs)){
					nb_next_head[p].x=xp;
					nb_next_head[p].y=nb_head[i].y;
					nb_next_head[p].z=nb_head[i].z;
					nb_next_orig[p].x=nb_orig[i].x;
					nb_next_orig[p].y=nb_orig[i].y;
					nb_next_orig[p].z=nb_orig[i].z;
					g->matrix[xp][nb_head[i].y][nb_head[i].z].from=i;//index of nb_head
					g->matrix[xp][nb_head[i].y][nb_head[i].z].dist=temp_dist;
					p++;
					
					if(g->matrix[xp][nb_head[i].y][nb_head[i].z].phi==1 && temp_dist==fr)
						p--;	//doesn't count
					g->matrix[xp][nb_head[i].y][nb_head[i].z].phi=1;
				}
			}

			if (g->matrix[xn][nb_head[i].y][nb_head[i].z].phi==0 || 
			   (g->matrix[xn][nb_head[i].y][nb_head[i].z].phi==1&&g->matrix[xn][nb_head[i].y][nb_head[i].z].from!=-1)){
				if(g->matrix[xn][nb_head[i].y][nb_head[i].z].phi==0)	fr=PR;
				else	fr=g->matrix[xn][nb_head[i].y][nb_head[i].z].dist;
				
				dx = g->matrix[xn][nb_head[i].y][nb_head[i].z].point.x - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.x;
				dy = g->matrix[xn][nb_head[i].y][nb_head[i].z].point.y - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.y;
				dz = g->matrix[xn][nb_head[i].y][nb_head[i].z].point.z - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.z;
				temp_dist=dx*dx+dy*dy+dz*dz;
				
				if(temp_dist<=fr && (fr-temp_dist >  abs)){
					nb_next_head[p].x=xn;
					nb_next_head[p].y=nb_head[i].y;
					nb_next_head[p].z=nb_head[i].z;
					nb_next_orig[p].x=nb_orig[i].x;
					nb_next_orig[p].y=nb_orig[i].y;
					nb_next_orig[p].z=nb_orig[i].z;
					g->matrix[xn][nb_head[i].y][nb_head[i].z].from=i;//index of nb_head
					g->matrix[xn][nb_head[i].y][nb_head[i].z].dist=temp_dist;
					p++;
					
					if(g->matrix[xn][nb_head[i].y][nb_head[i].z].phi==1 && temp_dist==fr)
						p--;
					g->matrix[xn][nb_head[i].y][nb_head[i].z].phi=1;
				}
			}

			if (g->matrix[nb_head[i].x][yp][nb_head[i].z].phi==0 || 
			   (g->matrix[nb_head[i].x][yp][nb_head[i].z].phi==1&&g->matrix[nb_head[i].x][yp][nb_head[i].z].from!=-1)){
				if(g->matrix[nb_head[i].x][yp][nb_head[i].z].phi==0)	fr=PR;
				else	fr=g->matrix[nb_head[i].x][yp][nb_head[i].z].dist;

				dx = g->matrix[nb_head[i].x][yp][nb_head[i].z].point.x - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.x;
				dy = g->matrix[nb_head[i].x][yp][nb_head[i].z].point.y - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.y;
				dz = g->matrix[nb_head[i].x][yp][nb_head[i].z].point.z - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.z;
				temp_dist=dx*dx+dy*dy+dz*dz;
				
				if(temp_dist<=fr && (fr-temp_dist >  abs)){
					nb_next_head[p].x=nb_head[i].x;
					nb_next_head[p].y=yp;
					nb_next_head[p].z=nb_head[i].z;
					nb_next_orig[p].x=nb_orig[i].x;
					nb_next_orig[p].y=nb_orig[i].y;
					nb_next_orig[p].z=nb_orig[i].z;
					g->matrix[nb_head[i].x][yp][nb_head[i].z].from=i;//index of nb_head
					g->matrix[nb_head[i].x][yp][nb_head[i].z].dist=temp_dist;
					p++;
					if(g->matrix[nb_head[i].x][yp][nb_head[i].z].phi==1 && temp_dist==fr)
						p--;
					g->matrix[nb_head[i].x][yp][nb_head[i].z].phi=1;
				}
			}

			if (g->matrix[nb_head[i].x][yn][nb_head[i].z].phi==0 || 
			   (g->matrix[nb_head[i].x][yn][nb_head[i].z].phi==1&&g->matrix[nb_head[i].x][yn][nb_head[i].z].from!=-1)){
				if(g->matrix[nb_head[i].x][yn][nb_head[i].z].phi==0)	fr=PR;
				else	fr=g->matrix[nb_head[i].x][yn][nb_head[i].z].dist;
				
				dx = g->matrix[nb_head[i].x][yn][nb_head[i].z].point.x - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.x;
				dy = g->matrix[nb_head[i].x][yn][nb_head[i].z].point.y - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.y;
				dz = g->matrix[nb_head[i].x][yn][nb_head[i].z].point.z - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.z;
				temp_dist=dx*dx+dy*dy+dz*dz;
				
				if(temp_dist<=fr && (fr-temp_dist >  abs)){
					nb_next_head[p].x=nb_head[i].x;
					nb_next_head[p].y=yn;
					nb_next_head[p].z=nb_head[i].z;
					nb_next_orig[p].x=nb_orig[i].x;
					nb_next_orig[p].y=nb_orig[i].y;
					nb_next_orig[p].z=nb_orig[i].z;
					g->matrix[nb_head[i].x][yn][nb_head[i].z].from=i;//index of nb_head
					g->matrix[nb_head[i].x][yn][nb_head[i].z].dist=temp_dist;
					p++;
					if(g->matrix[nb_head[i].x][yn][nb_head[i].z].phi==1 && temp_dist==fr)
						p--;
					g->matrix[nb_head[i].x][yn][nb_head[i].z].phi=1;
				}
			}
			
			if (g->matrix[nb_head[i].x][nb_head[i].y][zp].phi==0 || 
			   (g->matrix[nb_head[i].x][nb_head[i].y][zp].phi==1&&g->matrix[nb_head[i].x][nb_head[i].y][zp].from!=-1)){
				if(g->matrix[nb_head[i].x][nb_head[i].y][zp].phi==0)	fr=PR;
				else	fr=g->matrix[nb_head[i].x][nb_head[i].y][zp].dist;
				
				dx = g->matrix[nb_head[i].x][nb_head[i].y][zp].point.x - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.x;
				dy = g->matrix[nb_head[i].x][nb_head[i].y][zp].point.y - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.y;
				dz = g->matrix[nb_head[i].x][nb_head[i].y][zp].point.z - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.z;
				temp_dist=dx*dx+dy*dy+dz*dz;
				
				if(temp_dist<=fr && (fr-temp_dist >  abs)){
					nb_next_head[p].x=nb_head[i].x;
					nb_next_head[p].y=nb_head[i].y;
					nb_next_head[p].z=zp;
					nb_next_orig[p].x=nb_orig[i].x;
					nb_next_orig[p].y=nb_orig[i].y;
					nb_next_orig[p].z=nb_orig[i].z;
					g->matrix[nb_head[i].x][nb_head[i].y][zp].from=i;//index of nb_head
					g->matrix[nb_head[i].x][nb_head[i].y][zp].dist=temp_dist;
					p++;
					if(g->matrix[nb_head[i].x][nb_head[i].y][zp].phi==1 && temp_dist==fr)
						p--;
					g->matrix[nb_head[i].x][nb_head[i].y][zp].phi=1;
				}
			}
			
			if (g->matrix[nb_head[i].x][nb_head[i].y][zn].phi==0 || 
			   (g->matrix[nb_head[i].x][nb_head[i].y][zn].phi==1&&g->matrix[nb_head[i].x][nb_head[i].y][zn].from!=-1)){
				if(g->matrix[nb_head[i].x][nb_head[i].y][zn].phi==0)	fr=PR;
				else	fr=g->matrix[nb_head[i].x][nb_head[i].y][zn].dist;
				
				dx = g->matrix[nb_head[i].x][nb_head[i].y][zn].point.x - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.x;
				dy = g->matrix[nb_head[i].x][nb_head[i].y][zn].point.y - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.y;
				dz = g->matrix[nb_head[i].x][nb_head[i].y][zn].point.z - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.z;
				temp_dist=dx*dx+dy*dy+dz*dz;
				
				if(temp_dist<=fr && (fr-temp_dist >  abs)){
					nb_next_head[p].x=nb_head[i].x;
					nb_next_head[p].y=nb_head[i].y;
					nb_next_head[p].z=zn;
					nb_next_orig[p].x=nb_orig[i].x;
					nb_next_orig[p].y=nb_orig[i].y;
					nb_next_orig[p].z=nb_orig[i].z;
					g->matrix[nb_head[i].x][nb_head[i].y][zn].from=i;//index of nb_head
					g->matrix[nb_head[i].x][nb_head[i].y][zn].dist=temp_dist;
					p++;
					if(g->matrix[nb_head[i].x][nb_head[i].y][zn].phi==1 && temp_dist==fr)
						p--;
					g->matrix[nb_head[i].x][nb_head[i].y][zn].phi=1;
				}
			}
		}
		free(nb_head);
		free(nb_orig);
		nb_head = nb_next_head;
		nb_orig = nb_next_orig;
		co=p;
	}
	free(nb_head);
	free(nb_orig);
}
	
int fastMarching(lGrid g, char inner)
{
	int len = g->N;
	int i;
	int p;
	int x,y,z;
	PPoint* nb_head;
	PPoint* tmp;
	int co;
	int vol = 0;
	
	i=0;
	p=len-1;
	//initialize outer surface and narrow band
	for (x=i;x<=p;x++)		
		for (y=i;y<=p;y++){
			g->matrix[x][y][i].phi = 3;
			g->matrix[x][y][p].phi = 3;
			g->matrix[x][i][y].phi = 3;
			g->matrix[x][p][y].phi = 3;
			g->matrix[i][x][y].phi = 3;
			g->matrix[p][x][y].phi = 3;
		}
	i++; p--;
//because we use + 1, we directly set a bend of boundary
	co = 0;
	for (x=i;x<=p;x++)		
		for (y=i;y<=p;y++){
			g->matrix[x][y][i].phi = 3;
			g->matrix[x][y][p].phi = 3;
			g->matrix[x][i][y].phi = 3;
			g->matrix[x][p][y].phi = 3;
			g->matrix[i][x][y].phi = 3;
			g->matrix[p][x][y].phi = 3;
			co+=6;
		}

	nb_head = (PPoint*)malloc(co*sizeof(PPoint));

	co = 0;
	for (x=i;x<=p;x++)		
		for (y=i;y<=p;y++){
			nb_head[co].x = x;
			nb_head[co].y = y;
			nb_head[co++].z = i;

			nb_head[co].x = x;
			nb_head[co].y = y;
			nb_head[co++].z = p;

			nb_head[co].x = x;
			nb_head[co].y = i;
			nb_head[co++].z = y;

			nb_head[co].x = x;
			nb_head[co].y = p;
			nb_head[co++].z = y;

			nb_head[co].x = i;
			nb_head[co].y = x;
			nb_head[co++].z = y;

			nb_head[co].x = p;
			nb_head[co].y = x;
			nb_head[co++].z = y;
		}

	while (co!=0){
		p = 0;
		tmp = (PPoint*)malloc(co*6*sizeof(PPoint));	//at most 6 times larger than nb_head
		for (i=0;i<co;i++){
			if (g->matrix[nb_head[i].x+1][nb_head[i].y][nb_head[i].z].phi==1){
				g->matrix[nb_head[i].x+1][nb_head[i].y][nb_head[i].z].phi=3;
				tmp[p].x=nb_head[i].x+1;
				tmp[p].y=nb_head[i].y;
				tmp[p++].z=nb_head[i].z;
			}
			if (g->matrix[nb_head[i].x-1][nb_head[i].y][nb_head[i].z].phi==1){
				g->matrix[nb_head[i].x-1][nb_head[i].y][nb_head[i].z].phi=3;
				tmp[p].x=nb_head[i].x-1;
				tmp[p].y=nb_head[i].y;
				tmp[p++].z=nb_head[i].z;
			}
			if (g->matrix[nb_head[i].x][nb_head[i].y+1][nb_head[i].z].phi==1){
				g->matrix[nb_head[i].x][nb_head[i].y+1][nb_head[i].z].phi=3;
				tmp[p].x=nb_head[i].x;
				tmp[p].y=nb_head[i].y+1;
				tmp[p++].z=nb_head[i].z;
			}
			if (g->matrix[nb_head[i].x][nb_head[i].y-1][nb_head[i].z].phi==1){
				g->matrix[nb_head[i].x][nb_head[i].y-1][nb_head[i].z].phi=3;
				tmp[p].x=nb_head[i].x;
				tmp[p].y=nb_head[i].y-1;
				tmp[p++].z=nb_head[i].z;
			}
			if (g->matrix[nb_head[i].x][nb_head[i].y][nb_head[i].z+1].phi==1){
				g->matrix[nb_head[i].x][nb_head[i].y][nb_head[i].z+1].phi=3;
				tmp[p].x=nb_head[i].x;
				tmp[p].y=nb_head[i].y;
				tmp[p++].z=nb_head[i].z+1;
			}
			if (g->matrix[nb_head[i].x][nb_head[i].y][nb_head[i].z-1].phi==1){
				g->matrix[nb_head[i].x][nb_head[i].y][nb_head[i].z-1].phi=3;
				tmp[p].x=nb_head[i].x;
				tmp[p].y=nb_head[i].y;
				tmp[p++].z=nb_head[i].z-1;
			}
		}
		free(nb_head);
		nb_head = tmp;
		co = p;
	}
	free(nb_head);

	// outer surface
	if (!inner){
		for (x=0;x<len;x++)
		  for (y=0;y<len;y++)
		    for (z=0;z<len;z++)
			if (g->matrix[x][y][z].phi==3)
				g->matrix[x][y][z].phi=1;
			else{
				g->matrix[x][y][z].phi=0;
				vol++;
			}
	}else{
	// inner cave
		for (x=0;x<len;x++)
		  for (y=0;y<len;y++)
		    for (z=0;z<len;z++)
			if (g->matrix[x][y][z].phi!=1)
				g->matrix[x][y][z].phi=1;
			else{
				g->matrix[x][y][z].phi=0;
				vol++;
			}
	}
	return vol;
}
