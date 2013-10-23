typedef struct p
{
	short x;
	short y;
	short z;
} PPoint;

typedef struct gp
{
	PPoint point;
	char phi;
	int from;
	float dist;
} GridPoint;

typedef struct g
{
	short N;
	GridPoint*** matrix;
	short stepSize;
} *lGrid;
