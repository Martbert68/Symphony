#define X_SIZE 1920 
#define STRIDE 5760
#define Y_SIZE  1080 
#define Y_SIZ2  540 
#define MEM 100000
#define AP 35000 

void paint(unsigned char *,int);
void ripple(unsigned char *,int, int);
void apaint(unsigned char *,int);
void gaus(unsigned char *,int);
void blurt(unsigned char *,int);
void eartht(unsigned char *, int);
void plott(unsigned char *,int,int,int,int,int);
void linet( unsigned char *, int *, int , int , int , float , int , float ,float, float, float);
void lineu( unsigned char *, float *, float *,int r, int g, int b, float thick);
void merge(unsigned char *, unsigned char *, int);
void pline(unsigned char *, int );
int gline(int *, int *, int , int , int );
struct floc{float *x;float *y;float *z;float *v; float *mv;float *dx; float *dy; float *dz; float *dist; int *r; int *g; int *b; float *n; 
	float *tx; float *ty; float *tz; float *tdx; float *tdy; float *tdz; int *nn;};
struct ap{float *x;float *y;float *z;};
struct try{float *pi;float *pj;float *pk; float dist;};
struct triangle{float *pi;float *pj;float *pk;};
struct plate{float *pi;float *pj;float *pk; float *pl; float l; int r; int g; int b;};
struct real_plate{float *ri;float *rj;float *rk; float *rl; float *n; int r; int g; int b;};
struct view{float x; float y; float z; float pan; float tilt; float cp; float sp; float ct; float st; float l;};
int read_JPEG_file (char *, unsigned char *, int *);
int jayit(unsigned char *,int, int, char *);
void load_image( unsigned char *, char *, int );
void edge_detect (unsigned char *, struct ap *);
void draw_triangle (unsigned char *,struct triangle *,int,int,int);
void draw_plate (unsigned char *,struct plate *);
void measure_plate(struct real_plate *, struct plate *,struct view *,struct view *);
void shadow_plate(struct real_plate *, struct real_plate *,struct view *);
void draw_biod (unsigned char *,float *,float,int,int,int,float);
void plane_intersect(float *,float *,float *,float,float,float);
void three_two(float *,float *,float ,float);
int all_one(float *,float *, struct view *);
int t_plate(struct real_plate *, struct plate *, struct view *);
int make_box (float , float , float , float , float , float , float , float , float , struct real_plate *, int,int,int);
void twist (float *, float *,float);

//int rend_tri(struct try *,struct trf *, float *,float ,float , float, float,float );
void matrix(struct view *);


void measure_flock(struct floc *, float *, int, float,float,float,float);
int load_stl_flock(char *, struct floc *, int);
int obj_to_flock(struct floc *, struct floc *,int,int);
void sort_flock(struct floc *, int *, int);
int move_flock(unsigned char *,struct floc *, struct ap *, int *);
void calc_point(unsigned char *,struct floc *, struct ap *, int *);
void load_point(struct ap *,int);
void show_flock(unsigned char *,struct floc *, int);
void clear(unsigned char *);
void dist_tri(struct try *, float *, int);
void sort_tri(struct try *, int *, int);
void split (char *, char *,int , int , int );
long load_wav(short *samp,char *name);
