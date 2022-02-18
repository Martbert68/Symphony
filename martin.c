#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "martin.h"


long load_wav(short *samp,char *name)
{
        FILE *loadfile;


        int *fhead,chans,sample_rate,bits_pers,byte_rate,size;
        fhead=(int *)malloc(sizeof(int)*11);
        long ba;
	ba=0;

        chans=2;
        sample_rate=44100;
        bits_pers=16;
        byte_rate=(sample_rate*chans*bits_pers)/8;

        loadfile=fopen(name,"rb");
        if (loadfile !=NULL){
        	ba=fread(fhead,sizeof(int),11,loadfile);
        	size=2*fhead[1];
		int i;
        	printf ("loading %s header %ld %d\n",name,ba,size);
		for (i=0;i<11;i++){ printf ("I %d Val %d \n",i,fhead[i]);}
        	ba=fread(samp,sizeof(short),size,loadfile);
        	printf ("loading %s data  %ld\n",name,ba);
        	fclose (loadfile); } else { printf ("Can't load %s\n",name);}
        free(fhead);
        return ba ;
}



void split (char *image2, char *image4,int ang, int pha, int depth)
{
        float xrmove,xgmove,xbmove;
        float yrmove,ygmove,ybmove;
        int x,y;

        xrmove=(float)depth*(sin((2*M_PI*(float)ang)/360));
        yrmove=(float)depth*(cos((2*M_PI*(float)ang)/360));

        xgmove=(float)depth*(sin((2*M_PI*(float)(ang+pha))/360));
        ygmove=(float)depth*(cos((2*M_PI*(float)(ang+pha))/360));

        xbmove=(float)depth*(sin((2*M_PI*(float)(ang+(2*pha)))/360));
        ybmove=(float)depth*(cos((2*M_PI*(float)(ang+(2*pha)))/360));

        for (y=0;y<Y_SIZE;y++)
        {
                int p,r,ym;
                p=y*STRIDE;
                ym=y+yrmove;
                if (ym<0){ym=Y_SIZE+ym;} if (ym>=Y_SIZE){ym=Y_SIZE-ym;}
                r=ym*STRIDE;
                for (x=0;x<STRIDE;x+=3)
                {
                        int xm;
                        xm=x+(3*(int)xrmove);
                        if (xm<0){xm=STRIDE+xm;} if (xm>=STRIDE){xm=xm-STRIDE;}
                        image4[p+x]=image2[r+xm];
                }
        }
        for (y=0;y<Y_SIZE;y++)
        {
                int p,r,ym;
                p=y*STRIDE;
                ym=y+ygmove;
                if (ym<0){ym=0;} if (ym>=Y_SIZE){ym=Y_SIZE-1;}
                r=ym*STRIDE;
                for (x=1;x<STRIDE;x+=3)
                {
                        int xm;
                        xm=x+(3*(int)xgmove);
                        if (xm<0){xm=0;} if (xm>=STRIDE){xm=STRIDE-2;}
                        image4[p+x]=image2[r+xm];
                }
        }
        for (y=0;y<Y_SIZE;y++)
        {
                int p,r,ym;
                p=y*STRIDE;
                ym=y+ybmove;
                if (ym<0){ym=0;} if (ym>=Y_SIZE){ym=Y_SIZE-1;}
                r=ym*STRIDE;
                for (x=2;x<STRIDE;x+=3)
                {
                        int xm;
                        xm=x+(3*(int)xbmove);
                        if (xm<0){xm=0;} if (xm>=STRIDE){xm=STRIDE-1;}
                        image4[p+x]=image2[r+xm];
                }
        }
        memcpy(image2,image4,STRIDE*Y_SIZE);
}






int obj_to_flock(struct floc *obj, struct floc *f, int tri_count, int all)
{
	int along;
	int *done;
	done=(int *)malloc(sizeof(int)*all);

        // Lets preallocate
	//
	
        for (along=0;along<all;along++){
                int sel,un,aa;
                un=1;
		if (tri_count<all){ un=0; sel=rand()%tri_count;}
                while (un)
                {
                        un=0;
                        sel=rand()%tri_count;
                        for (aa=0;aa<along;aa++)
                        {
                                if (sel==done[aa]){ un=1; printf("eh\n");}
                        }
                }
                done[along]=sel;

                f->tx[along]=obj->x[sel];
                f->ty[along]=obj->y[sel];
                f->tz[along]=obj->z[sel];
                f->tdx[along]=obj->dx[sel];
                f->tdy[along]=obj->dy[sel];
                f->tdz[along]=obj->dz[sel];
        }
	free (done);
}


int load_stl_flock(char *name, struct floc *obj, int shrink)  
{

	FILE *stl;

       stl=fopen(name,"r");
       if (stl==NULL){exit(1);}

       short attrib;
       float normal[3];
       float xyz[3];
       int ll,loop,along;
       loop=0;
       while ( fread(normal,sizeof(float),3,stl))
       {
                ll=loop/shrink;
		float nn;
		nn=sqrt((normal[0]*normal[0])+(normal[1]*normal[1])+(normal[2]*normal[2]));
                obj->dx[ll]=normal[0]/nn;
                obj->dy[ll]=normal[1]/nn;
                obj->dz[ll]=normal[2]/nn;
                fread(xyz,sizeof(float),3,stl); // 3 tri cords
		//printf ("%f %f %f \n",xyz[0],xyz[1],xyz[2]);
		obj->x[ll]=xyz[0];
		obj->y[ll]=xyz[2];
		obj->z[ll]=xyz[1];
                fread(xyz,sizeof(float),3,stl);
                fread(xyz,sizeof(float),3,stl);
                fread(&attrib,sizeof(short),1,stl);
                loop++;
        }

        printf ("I reckon there are %d trianges and I read %d\n",loop,ll);
        fclose(stl);

        // We have loaded the file!

        float xmax,ymax,zmax; float xmin,ymin,zmin;
        xmax=-1000000;ymax=-1000000;zmax=-1000000; xmin=1000000;ymin=1000000;zmin=1000000;
        for (along=0;along<ll;along++)
        {
                if (obj->x[along]>xmax){xmax=obj->x[along];} if (obj->x[along]<xmin){xmin=obj->x[along];} 
                if (obj->y[along]>ymax){ymax=obj->y[along];} if (obj->y[along]<ymin){ymin=obj->y[along];} 
                if (obj->z[along]>zmax){zmax=obj->z[along];} if (obj->z[along]<zmin){zmin=obj->z[along];} 
	}

        //printf ("xmax %f xmin %f ymax %f ymin %f zmax %f zmin %f \n",xmax,xmin,ymax,ymin,zmax,zmin);

	// centre everything round the origin +/-500 in the longest direction.
	float xspan,yspan,zspan,span,scale;
	span=1000;
	xspan=xmax-xmin; yspan=ymax-ymin; zspan=zmax-zmin;
	if (xspan>yspan && xspan>zspan){  scale=span/xspan;}
	if (yspan>xspan && yspan>zspan){  scale=span/yspan;}
	if (zspan>yspan && zspan>xspan){  scale=span/zspan;}



        // as a guess lets scale xmax to scale and centre .
	// take away X min and you are from the origin tak eaway half the span you are centred aroud the origin
        for (along=0;along<ll;along++)
        {
                obj->x[along]=scale*(obj->x[along]-xmin-(xspan/2));
                obj->y[along]=scale*(obj->y[along]-ymin-(yspan/2));
                obj->z[along]=scale*(obj->z[along]-zmin-(zspan/2));
        }

	return ll;
}



void sort_flock(struct floc *in, int *order, int count)
{
        // first thing lets bubble sort the z values.
        int sorting;
        sorting=1;
        while (sorting)
        {
                sorting=0;
		int m;
                for (m=0;m<count-1;m++)
                {
                        if (in->dist[order[m]] < in->dist[order[m+1]])
                        {
                        	int temp;
                                sorting=1;
                                temp=order[m];
                                order[m]=order[m+1];
                                order[m+1]=temp;
                        }
                }
        }
}	


void measure_flock(struct floc *in,float *view, int count, float cp, float sp, float ct, float st)
{
	int along;

	for (along=0;along<count;along++)
	{
		float x,y,z,dx,dy,dz,l,nx,ny,nz,dot;
		x=in->x[along];
		y=in->y[along];
		z=in->z[along];
		dx=x-view[0];
		dy=y-view[1];
		dz=z-view[2];
		in->dist[along]=sqrt((dx*dx)+(dy*dy)+(dz*dz));
		//l=sqrt(((x*x))+((y*y))+(z*z));
		//in->dx[along]=(x/l);in->dy[along]=(y/l);in->dz[along]=(z/l);	
		nx=ct*cp;
		ny=ct*sp;
		nz=st;
		dot=((nx*(in->dx[along]))+(ny*(in->dy[along]))+(nz*in->dz[along]));
		if (dot>0){ in->n[along]=dot;}else{in->n[along]=-dot;}
	}
}

void sort_tri(struct try *in,int *order,int count)
{
	// bubble sort the order.
	// order[0] furtherst away
	int complete,n;
	complete=0;
	n=0;
	while ( !complete)
	{
		int along;
		complete=1;
		for (along=0;along<count-1;along++)
		{
			int pp,pq,pop;	
			float dp,dq;
			pp=order[along];
			pq=order[along+1];
			dp=(in+pp)->dist;
			dq=(in+pq)->dist;
			//printf ("%f %f \n",dp,dq);
			if (dp<dq){ order[along]=pq;order[along+1]=pp;complete=0;}
		}
		printf ("eeek! %d\n",n++);
	}
}

void dist_tri(struct try *in,float *view,int count)
{
	int along;
	for (along=0;along<count;along++)
	{
		float x,y,z;
		x=(in+along)->pi[0]-view[0];
		y=(in+along)->pj[0]-view[1];
		z=(in+along)->pk[0]-view[2];
		(in+along)->dist=sqrt((x*x)+(y*y)+(z*z));
	}
}

/*
int rend_tri(struct try *inn, struct trf *out,float *view, float cp,float sp, float ct, float st, float l)
{
	float in[3];
	float xy[2];
	in[0]=inn->pi[0];in[1]=inn->pi[1];in[2]=inn->pi[2];
        if (! all_one(in,xy,view,cp,sp,ct,st,l)){ printf("fail \n");return 0;}
        out->pi[0]=xy[0]; out->pi[1]=xy[1];

	in[0]=inn->pj[0];in[1]=inn->pj[1];in[2]=inn->pj[2];
        if (! all_one(in,xy,view,cp,sp,ct,st,l)){ printf("fail \n");return 0;}
        out->pj[0]=xy[0]; out->pj[1]=xy[1];

	in[0]=inn->pk[0];in[1]=inn->pk[1];in[2]=inn->pk[2];
        if (! all_one(in,xy,view,cp,sp,ct,st,l)){ printf("fail \n");return 0;}
        out->pk[0]=xy[0]; out->pk[1]=xy[1];

	return 1;
}*/	

void shadow_plate( struct real_plate *rp, struct real_plate *sp,struct view *light)
{
	sp->ri[0]=0; sp->ri[1]=0; sp->ri[2]=0;
	sp->rj[0]=200; sp->rj[1]=0; sp->rj[2]=0;
	sp->rk[0]=200; sp->rk[1]=200; sp->rk[2]=0;
	sp->rl[0]=0; sp->rl[1]=200; sp->rl[2]=0;

	sp->r=0;
	sp->g=0;
	sp->b=0;
}


void measure_plate( struct real_plate *rp, struct plate *p,struct view *v, struct view *i)
{
	float mx,my,mz,d,dx,dy,dz,dt,nx,ny,nz,ix,iy,iz,rx,ry,rz,in,bright;
	dx=rp->ri[0]-v->x; dy=rp->ri[1]-v->y; dz=rp->ri[2]-v->z;
	mx=rp->ri[0]; my=rp->ri[1]; mz=rp->ri[2];
	d=(dx*dx)+(dy*dy)+(dz*dz);

	dx=rp->rj[0]-v->x; dy=rp->rj[1]-v->y; dz=rp->rj[2]-v->z;
	mx+=rp->rj[0]; my+=rp->rj[1]; mz+=rp->rj[2];
	dt=(dx*dx)+(dy*dy)+(dz*dz);
       	if (dt>d){ d=dt;}

	dx=rp->rk[0]-v->x; dy=rp->rk[1]-v->y; dz=rp->rk[2]-v->z;
	mx+=rp->rk[0]; my+=rp->rk[1]; mz+=rp->rk[2];
	dt=(dx*dx)+(dy*dy)+(dz*dz);
       	if (dt>d){ d=dt;}

	dx=rp->rl[0]-v->x; dy=rp->rl[1]-v->y; dz=rp->rl[2]-v->z;
	mx+=rp->rl[0]; my+=rp->rl[1]; mz+=rp->rl[2];
	dt=(dx*dx)+(dy*dy)+(dz*dz);
       	if (dt>d){ d=dt;}
	p->l=sqrt(d);

	mx/=4;my/=4;mz/=4;
	// normal is the plate normal
	nx=rp->n[0];
	ny=rp->n[1];
	nz=rp->n[2];

	/* PARALLEL INCIDENT
	ix=cos(i->tilt)*cos(i->pan);
	iy=cos(i->tilt)*sin(i->pan);
	iz=sin(i->tilt); */

	/* point incident */
	/* mid plate */

	dx=mx-i->x;
	dy=my-i->y;
	dz=mz-i->z;
	dt=sqrt((dx*dx)+(dy*dy)+(dz*dz));

	ix=-dx/dt;
	iy=-dy/dt;
	iz=-dz/dt;
	//printf ("ix %f iy %f iz %f\n",ix,iy,iz);

	// ref = incident - 2(i.n)n
	in=2*((ix*nx)+(iy*ny)+(iz*nz));
	rx=ix-(in*nx);
	ry=iy-(in*ny);
	rz=iz-(in*nz);

	// normal is now normal to view 
	nx=cos(v->tilt)*cos(v->pan);
	ny=cos(v->tilt)*sin(v->pan);
	nz=sin(v->tilt); 


	bright=(rx*nx)+ (ry*ny)+ (rz*nz);

	//printf ("nx %f ny %f nz %f ix %f iy %f iz %f rx %f ry %f rz %f bright %f \n",nx,ny,nz,ix,iy,iz,rx,ry,rz,bright);

	bright=-bright;
	if (bright<0){bright=0;}
	bright=0.2+(4*bright/5);
	if (bright>1){ printf ("Bright %f \n",bright);bright=1;}
	p->r=(float)rp->r*(bright);
	p->g=(float)rp->g*bright;
	p->b=(float)rp->b*bright;
}

void twist ( float *x, float *y, float ang)
{
	float temp;
	temp=x[0];
	x[0]=(x[0]*cos(ang))-(y[0]*sin(ang));
	y[0]=(y[0]*cos(ang))+(temp*sin(ang));
}



int make_box (float x, float y, float z, float theta, float delta, float phi, float h, float w, float d, struct real_plate *rp, int r, int g, int b)
{
	float px,py,pz,l;


	px=w;
	py=0;
	pz=0;

	twist(&py,&pz,-theta);
	twist(&pz,&px,-delta);
	twist(&py,&px,-phi);

	float *xp,*yp,*zp;
	xp=(float *)malloc(sizeof(float)*4);
	yp=(float *)malloc(sizeof(float)*4);
	zp=(float *)malloc(sizeof(float)*4);

	xp[0]=0;xp[1]=0;xp[2]=0;xp[3]=0;
	yp[0]=-d;yp[1]=-d;yp[2]=d;yp[3]=d;
	zp[0]=-h;zp[1]=h;zp[2]=h;zp[3]=-h;

	//theta
	twist (yp,zp,theta);
	twist (yp+1,zp+1,theta);
	twist (yp+2,zp+2,theta);
	twist (yp+3,zp+3,theta);

	//delta
	twist (xp,zp,delta);
	twist (xp+1,zp+1,delta);
	twist (xp+2,zp+2,delta);
	twist (xp+3,zp+3,delta);

	//phi
	twist (xp,yp,phi);
	twist (xp+1,yp+1,phi);
	twist (xp+2,yp+2,phi);
	twist (xp+3,yp+3,phi);

	// x y z FRONT plate 0
	rp->ri[0]=x+px+xp[0];
	rp->ri[1]=y+py+yp[0];
	rp->ri[2]=z+pz+zp[0];
	//rp->ri[2]=z-h;
	// x y z FRONT plate 1
	rp->rj[0]=x+px+xp[1];
	rp->rj[1]=y+py+yp[1];
	rp->rj[2]=z+pz+zp[1];
	//rp->rj[2]=z+h;
	// x y z FRONT plate 2
	rp->rk[0]=x+px+xp[2];
	rp->rk[1]=y+py+yp[2];
	rp->rk[2]=z+pz+zp[2];
	// x y z FRONT plate 3
	rp->rl[0]=x+px+xp[3];
	rp->rl[1]=y+py+yp[3];
	rp->rl[2]=z+pz+zp[3];
	// normals
	rp->n[0]=cos(delta)*cos(phi);
	rp->n[1]=cos(theta)*sin(phi);
	rp->n[2]=cos(theta)*sin(delta);
	//rp->n[0]=px/w;
	//rp->n[1]=py/w;
	//rp->n[2]=pz/w;

	// x y z Back plate 4
	(rp+1)->ri[0]=x-px+xp[0];
	(rp+1)->ri[1]=y-py+yp[0];
	(rp+1)->ri[2]=z-pz+zp[0];
	// x y z Back plate 5
	(rp+1)->rj[0]=x-px+xp[1];
	(rp+1)->rj[1]=y-py+yp[1];
	(rp+1)->rj[2]=z-pz+zp[1];
	// x y z Back plate 6
	(rp+1)->rk[0]=x-px+xp[2];
	(rp+1)->rk[1]=y-py+yp[2];
	(rp+1)->rk[2]=z-pz+zp[2];
	// x y z Back plate 7
	(rp+1)->rl[0]=x-px+xp[3];
	(rp+1)->rl[1]=y-py+yp[3];
	(rp+1)->rl[2]=z-pz+zp[3];
	// normals
	(rp+1)->n[0]=-rp->n[0];
	(rp+1)->n[1]=-rp->n[1];
	(rp+1)->n[2]=-rp->n[2];

	//1/5/6/2 TOP plate 2
	(rp+2)->ri[0]=rp->rj[0]; (rp+2)->ri[1]=rp->rj[1]; (rp+2)->ri[2]=rp->rj[2];
	(rp+2)->rj[0]=(rp+1)->rj[0]; (rp+2)->rj[1]=(rp+1)->rj[1]; (rp+2)->rj[2]=(rp+1)->rj[2];
	(rp+2)->rk[0]=(rp+1)->rk[0]; (rp+2)->rk[1]=(rp+1)->rk[1]; (rp+2)->rk[2]=(rp+1)->rk[2];
	(rp+2)->rl[0]=rp->rk[0]; (rp+2)->rl[1]=rp->rk[1]; (rp+2)->rl[2]=rp->rk[2];
	// normals
	(rp+2)->n[0]=-sin(delta)*cos(phi);
	(rp+2)->n[1]=-sin(theta)*cos(phi);
	(rp+2)->n[2]=cos(theta)*cos(delta);

	//0/4/7/3 Bottom plate 3
	(rp+3)->ri[0]=rp->ri[0]; (rp+3)->ri[1]=rp->ri[1]; (rp+3)->ri[2]=rp->ri[2];
	(rp+3)->rj[0]=(rp+1)->ri[0]; (rp+3)->rj[1]=(rp+1)->ri[1]; (rp+3)->rj[2]=(rp+1)->ri[2];
	(rp+3)->rk[0]=(rp+1)->rl[0]; (rp+3)->rk[1]=(rp+1)->rl[1]; (rp+3)->rk[2]=(rp+1)->rl[2];
	(rp+3)->rl[0]=rp->rl[0]; (rp+3)->rl[1]=rp->rl[1]; (rp+3)->rl[2]=rp->rl[2];
	// normals
	(rp+3)->n[0]=-(rp+2)->n[0];
	(rp+3)->n[1]=-(rp+2)->n[1];
	(rp+3)->n[2]=-(rp+2)->n[2];

	//0/1/5/4 left plate 4
	(rp+4)->ri[0]=rp->ri[0]; (rp+4)->ri[1]=rp->ri[1]; (rp+4)->ri[2]=rp->ri[2];
	(rp+4)->rj[0]=rp->rj[0]; (rp+4)->rj[1]=rp->rj[1]; (rp+4)->rj[2]=rp->rj[2];
	(rp+4)->rk[0]=(rp+1)->rj[0]; (rp+4)->rk[1]=(rp+1)->rj[1]; (rp+4)->rk[2]=(rp+1)->rj[2];
	(rp+4)->rl[0]=(rp+1)->ri[0]; (rp+4)->rl[1]=(rp+1)->ri[1]; (rp+4)->rl[2]=(rp+1)->ri[2];
	// normals
	(rp+4)->n[0]=-cos(delta)*sin(phi);
	(rp+4)->n[1]=cos(theta)*cos(phi);
	(rp+4)->n[2]=sin(theta)*cos(delta);

	//3/2/6/7 right plate 5
	(rp+5)->ri[0]=rp->rl[0]; (rp+5)->ri[1]=rp->rl[1]; (rp+5)->ri[2]=rp->rl[2];
	(rp+5)->rj[0]=rp->rk[0]; (rp+5)->rj[1]=rp->rk[1]; (rp+5)->rj[2]=rp->rk[2];
	(rp+5)->rk[0]=(rp+1)->rk[0]; (rp+5)->rk[1]=(rp+1)->rk[1]; (rp+5)->rk[2]=(rp+1)->rk[2];
	(rp+5)->rl[0]=(rp+1)->rl[0]; (rp+5)->rl[1]=(rp+1)->rl[1]; (rp+5)->rl[2]=(rp+1)->rl[2];
	// normals
	(rp+5)->n[0]=-(rp+4)->n[0];
	(rp+5)->n[1]=-(rp+4)->n[1];
	(rp+5)->n[2]=-(rp+4)->n[2];

	int i;
	for (i=0;i<6;i++)
	{
		(rp+i)->r=r;
		(rp+i)->g=g;
		(rp+i)->b=b;
	}

	free (xp);
	free (yp);
	free (zp);

}



int t_plate(struct real_plate *rp, struct plate *p, struct view *v) 
{
	float in[3];
	float out[2];
	in[0]=rp->ri[0]; in[1]=rp->ri[1]; in[2]=rp->ri[2];
	all_one(in,out, v);
	p->pi[0]=out[0]; p->pi[1]=out[1];

	in[0]=rp->rj[0]; in[1]=rp->rj[1]; in[2]=rp->rj[2];
	all_one(in,out, v);
	p->pj[0]=out[0]; p->pj[1]=out[1];


	in[0]=rp->rk[0]; in[1]=rp->rk[1]; in[2]=rp->rk[2];
	all_one(in,out, v);
	p->pk[0]=out[0]; p->pk[1]=out[1];

	in[0]=rp->rl[0]; in[1]=rp->rl[1]; in[2]=rp->rl[2];
	all_one(in,out, v);
	p->pl[0]=out[0]; p->pl[1]=out[1];
}


void matrix (struct view *p)
{
	p->cp=cos(p->pan);
	p->sp=sin(p->pan);
	p->ct=cos(p->tilt);
	p->st=sin(p->tilt);
}


int all_one(float *in,float *out, struct view *v)
{

	float dx;
	float dy; 
	float dz; 
	float tt,X,Y,Z;

	//nx=-ct*cp;
	//ny=ct*sp;
	//nz=-st;

	dx=(in[0]-v->x);
	dy=(in[1]-v->y);
	dz=(in[2]-v->z);

	//printf("nx %f ny  %f nz %f\n",nx,ny,nz);
	//nl=sqrt((nx*nx)+(ny*ny)+(nz*nz));

	//xp=xv+(l*(ct*cp));
	///yp=yv-(l*(ct*sp));
	//zp=zv+(l*st);

	//tt=((nx*(xp-xv))+(ny*(yp-yv))+(nz*(zp-zv)))/(((nx*(x1-xv))+(ny*(y1-yv))+(nz*(z1-zv))));
	tt=((((v->ct*v->cp*v->ct*v->cp)))+(v->ct*v->sp*v->ct*v->sp)+(v->st*v->st))/(-(v->ct*v->cp*dx)-(v->ct*v->sp*dy)+(v->st*dz));
	if (tt<0 || tt>v->l){ return 0;}

	X=v->l*((tt*dx)-(v->ct*v->cp));
	Y=v->l*((v->ct*v->sp)-(tt*dy));
	Z=v->l*((tt*dz)+v->st);


	out[0]=(v->sp*X)+(v->cp*Y)+(X_SIZE/2);
	out[1]=(Y_SIZE/2)-(((v->cp*X)-(v->sp*Y))*v->st)-(Z*v->ct);


	if (out[0]<-100 || out[0]>X_SIZE+100 || out[1]< -100 || out[1]>Y_SIZE+100){ return 0;}else{return 1;}
}


void plane_intersect(float *in,float *out,float *view, float pan,float tilt, float l)
{
	// Compute the plane.
	//
	//point on the plane.
	//l in the x/y plane 
	//

	float view_x,view_y,view_z;
	float in_x,in_y,in_z,t;
	float x0_plane,y0_plane,z0_plane;
	float nx_plane,ny_plane,nz_plane;

	view_x=view[0]; view_y=view[1]; view_z=view[2];
	in_x=in[0]; in_y=in[1]; in_z=in[2];

	//point on plane
	x0_plane=view_x+(l*cos(tilt)*cos(pan));
	y0_plane=view_y+(l*cos(tilt)*sin(pan));
	z0_plane=view_z+(l*sin(tilt));


	//normal
	nx_plane=-cos(tilt)*cos(pan);
	ny_plane=-cos(tilt)*sin(pan);
	nz_plane=-sin(tilt);


	//plane equation
	//nx_plane(x-x0_plane)+ny_plane(y-y0_plane)+nz_plane(z-z0_plane)=0;

	//line equation
	//x=view_x+t*(in_x-view_x);
	//y=view_y+t*(in_y-view_y);
	//z=view_z+t*(in_z-view_z);

	// t
	//nx_plane(view_x+t*(in_x-view_x)-x0_plane)+ny_plane(view_x+t*(in_x-view_x)-y0_plane)+nz_plane(view_x+t*(in_x-view_x)-z0_plane)=0;
	// just x the rest the same
	//nx_plane*(view_x+(t*(in_x-view_x)))-(nx_plane*x0_plane)
	//(nx_plane*(view_x-x0_plane))+(nx_plane*t)*(in_x-view_x)
	//(nx_plane*t)*(in_x-view_x)=nx_plane*(x0_plane-view_x)
	t=((nx_plane*(x0_plane-view_x))+(ny_plane*(y0_plane-view_y))+(nz_plane*(z0_plane-view_z)))/
		((nx_plane*(in_x-view_x))+(ny_plane*(in_y-view_y))+(nz_plane*(in_z-view_z)));

	
	printf ("n len %f t %f\n",(nx_plane*nx_plane)+(ny_plane*ny_plane)+(nz_plane*nz_plane),t);

	out[0]=view_x+(t*(in_x-view_x))-x0_plane;
	out[1]=view_y+(t*(in_y-view_y))-y0_plane;
	out[2]=view_z+(t*(in_z-view_z))-z0_plane;


	//out[0]=view_x+(t*(in_x-view_x));
	//out[1]=view_y+(t*(in_y-view_y));
	//out[2]=view_z+(t*(in_z-view_z));

}

void three_two(float *in,float *out,float pan,float tilt)
{
	float  a,b,c,u,v,w,r,s,t,n,m,o;
	a=in[0];
	b=in[1];
	c=in[2];

	// flattened values.
	//tilt 0 = (x , y) panned
	//tilt 90 = z,(x,y panned)
	//

	// normal	
	n=-cos(tilt)*cos(pan);
	m=-cos(tilt)*sin(pan);
	o=-sin(tilt);

	//printf ("N %f\n",(n*n)+(m*m)+(o*o));

	// X vector
	u=-sin(pan);
	v=cos(pan);
	w=0;

	//printf ("X %f\n",(u*u)+(v*v)+(w*w));

	// Y vector
	r=-sin(tilt)*cos(pan);
	s=-sin(tilt)*sin(pan);
	t=-cos(tilt);

	//printf ("Y %f\n",(r*r)+(s*s)+(t*t));

	float theta, phi, x,y,delta,eta;

	//theta=((n*w)-(o*u));
	//phi=((n*v)-(m*u));

	//delta=((n*t)-(o*r));
	//eta=((n*s)-(m*r));

	//printf("theta phi delta eta %f %f %f %f \n",theta,phi,delta,eta);

	//printf ("Y2 %f\n",((b*u)-(a*v))/((u*s)-(v*r)));
	//printf ("X2 %f\n",((b*r)-(a*s))/((v*r)-(s*u)));

	//out[0]=((theta*((n*b)-(m*a)))+(phi*((o*a)-(c*n))))/((theta*((n*s)-(m*r)))+(phi*((o*r)-(n*t))));
	//out[1]=((eta*((n*c)-(o*a)))+(delta*((m*a)-(n*b))))/((eta*((n*w)-(o*u)))+(delta*((m*u)-(n*v))));
	out[0]=(((b*u)-(a*v))/((u*s)-(v*r)));
	out[1]=(((b*r)-(a*s))/((v*r)-(s*u)));
}


void draw_plate (unsigned char *painted, struct plate *p)
{
	struct triangle *t;
	float in[2];

	t=(struct triangle *)malloc(sizeof(struct triangle));
	t->pi=(float *)malloc(sizeof(float)*3);
	t->pj=(float *)malloc(sizeof(float)*3);
	t->pk=(float *)malloc(sizeof(float)*3);

	t->pi[0]=p->pi[0]; t->pi[1]=p->pi[1];
	t->pj[0]=p->pj[0]; t->pj[1]=p->pj[1];
	t->pk[0]=p->pk[0]; t->pk[1]=p->pk[1];
//void draw_biod (unsigned char *painted, float *in,float rad,int r,int g, int b, float n)
	draw_triangle (painted,t,p->r,p->g,p->b);
	//in[0]=t->pi[0]; in[1]=t->pi[1]; draw_biod (painted,in,10,0,0,255,1);
	//in[0]=t->pj[0]; in[1]=t->pj[1]; draw_biod (painted,in,10,0,255,0,1);
	t->pi[0]=p->pi[0]; t->pi[1]=p->pi[1];
	t->pj[0]=p->pk[0]; t->pj[1]=p->pk[1];
	t->pk[0]=p->pl[0]; t->pk[1]=p->pl[1];
	draw_triangle (painted,t,p->r,p->g,p->b);

	free (t->pi);
	free (t->pj);
	free (t->pk);
	free (t);
}	

void draw_triangle (unsigned char *painted, struct triangle *flat,int r,int g, int b)
{
	float lk,lj,dxj,dyj,dxk,dyk,along;
	float xi,xj,xk,yi,yj,yk,inc;

	xi=flat->pi[0]; yi=flat->pi[1];
	xj=flat->pj[0]; yj=flat->pj[1];
	xk=flat->pk[0]; yk=flat->pk[1];

	if (xi>X_SIZE || xj>X_SIZE || xk>X_SIZE){return;}
	if (yi>Y_SIZE || yj>Y_SIZE || yk>Y_SIZE){return;}
	if (xi<0 || xj<0 || xk<0 ){return;}
	if (yi<0 || yj<0 ||  yk<0 ){return;}

	//printf ("%f %f %f %f %f %f\n",xi,yi,xj,yj,xk,yk);

	dxj=xi-xj; dyj=yi-yj;
	dxk=xi-xk; dyk=yi-yk;

	lk=sqrt((dxk*dxk)+(dyk*dyk));
	lj=sqrt((dxj*dxj)+(dyj*dyj));

	if (lk>lj){ inc=lk; }else{inc=lj;}

	for (along=0;along<inc;along+=0.8)
	{
		float xb,yb,xe,ye,dx,dy,len,fill;
		xb=xi-((along*dxj)/inc);
		yb=yi-((along*dyj)/inc);
		xe=xi-((along*dxk)/inc);
		ye=yi-((along*dyk)/inc);
		dx=xb-xe;dy=yb-ye;
		len=sqrt((dx*dx)+(dy*dy));
		for (fill=0;fill<len;fill+=0.8)
		{
			float x,y;
			x=xb-(fill*dx/len);
			y=yb-(fill*dy/len);
			plott(painted,x,y,r,g,b);			
		}
	}	
}


void draw_biod (unsigned char *painted, float *in,float rad,int r,int g, int b, float n)
{
	float radius;
	float phi,x,y;
	int rm,gm,bm;

	rm=r*n;
	gm=g*n;
	bm=b*n;

	if (rad==1){
			plott(painted,in[0],in[1],rm,gm,bm);			
			plott(painted,in[0]-1,in[1]+1,rm,gm,bm);			
			plott(painted,in[0]-2,in[1]+2,rm,gm,bm);			
			plott(painted,in[0]+1,in[1]+1,rm,gm,bm);			
			plott(painted,in[0]+2,in[1]+2,rm,gm,bm);			
		}else{


	for (radius=0.4;radius<rad;radius+=0.4)
	{
		for (phi=0;phi<2*M_PI;phi+=M_PI/(4*radius))
		{
			x=in[0]+(radius*sin(phi));
			y=in[1]+(radius*cos(phi));
			plott(painted,x,y,rm,gm,bm);			
		}
	}} 
}


void edge_detect (unsigned char *painted, struct ap *a)
{
        long *escore;
        int *counter;

        escore=(long *)malloc(sizeof(long)*X_SIZE*Y_SIZE);
        counter=(int *)malloc(sizeof(int)*3*3*255);

        int x,y,total;
        long score;


        for (x=0;x<(3*3*255);x++){ counter[x]=0;}

        for (y=2;y<Y_SIZE-1;y++)
        {
                int yp,ym;
                yp=y*3*X_SIZE;
                ym=(y-1)*3*X_SIZE;
                for (x=2;x<X_SIZE-1;x++)
                {
                        int pix,pp,pm;
                        score=0;
                        pp=yp+(x*3);
                        pm=ym+(x*3);
                        for (pix=0;pix<3;pix++)
                        {
                                int b;
                                b=painted[pp+pix]-painted[pp+pix-3]; //left
                                score+=(b*b);
                                b=painted[pp+pix]-painted[pm+pix]; //up
                                score+=(b*b);
                                b=painted[pp+pix]-painted[pm+pix-3]; //up left
                                score+=(b*b);
                        }
                        score=sqrt(score);
                        escore[x+(y*X_SIZE)]=score;
                        counter[score]++;
                }
        }
        total=0;
        x=(3*3*255)-1;
        while (total<AP && x>=0)
        {
                total+=counter[x];
                score=x;
                x--;
        }
        score ++;
        printf ("paint %d %d %d\n",x,total,AP);

	int ap_count;
	ap_count=0;

        for (y=2;y<Y_SIZE-1;y++)
        {
                int yp,ep,xb,xe;
                int nr,ng,nb,along;
                yp=y*3*X_SIZE;
                for (x=2;x<X_SIZE-1;x++)
                {
                        ep=x+(y*X_SIZE);
                        if (escore[ep]>=score && ap_count<AP)
			{
				int t;
                            t=painted[yp+(x*3)]+painted[yp+(x*3)+1]+painted[yp+(x*3)+2];
			    a->x[ap_count]=x;
			    a->y[ap_count]=y;
			    a->z[ap_count]=t*Y_SIZE/(255*3);
                            //painted[yp+(x*3)]=255;
                            //painted[yp+(x*3)+1]=255;
                           // painted[yp+(x*3)+2]=255;
			    ap_count++;
                        } else {
                            //painted[yp+(x*3)]=0;
                            //painted[yp+(x*3)+1]=0;
                            //painted[yp+(x*3)+2]=0;
			}
                }
        }
        free (escore);
        free (counter);
}

void load_image( unsigned char *image, char *kk, int r)
{
        int x_size,y_size,x,y,pims[3];
        unsigned char *buff;
        buff=(unsigned char *)malloc(sizeof (char)*20*4096*4096); // image loader  buffer
        read_JPEG_file (kk, buff, pims);

        //streach it to fit both ways.
        float x_scale,y_scale,scale;

        x_size=pims[0];
        y_size=pims[1];

        x_scale=(float)x_size/(float)X_SIZE;
        y_scale=(float)y_size/(float)Y_SIZE;

        if (x_scale<=y_scale){ scale=x_scale;}else{scale=y_scale;}

        for (y=0;y<Y_SIZE;y++)
        {
                for (x=0;x<X_SIZE;x++)
                {
                        int b,i;
                        i=(x*3)+(y*X_SIZE*3);
                        b=((int)(x*scale)*3)+((int)(y*scale)*x_size*3);
                        image[i]=buff[b];
                        image[i+1]=buff[b+1];
                        image[i+2]=buff[b+2];
                }
        }
        free (buff);
}


void clear (unsigned char *painted)
{
	int x;
	for (x=0;x<X_SIZE*Y_SIZE*3;x++)
	{
		painted[x]=255;
	}
}

void load_point (struct ap *a, int gen)
{
	int m;
	if (gen%2==0)
	{
		for (m=0;m<AP;m++)
		{
			a->x[m]=rand()%X_SIZE;
			a->y[m]=rand()%Y_SIZE;
			a->z[m]=rand()%X_SIZE;
		}
	}else {
		for (m=0;m<AP;m++)
		{
			a->x[m]=X_SIZE/2+((X_SIZE/2)*sin(2*m*(float)(m)/AP));
			a->y[m]=Y_SIZE/2+((Y_SIZE/2)*cos((2+(float)gen/100)*m*(float)(m)/AP));
			a->z[m]=X_SIZE*m/AP;
		}
	}
	
}	


void calc_point (unsigned char *image, struct floc *f, struct ap *a, int *point)
{
	int m;
	int *visit;
	visit=(int *)malloc(sizeof(int)*AP);
	for (m=0;m<AP;m++){ visit[m]=0;}
	for (m=0;m<MEM;m++)
	{
		float min,dx,dy,dz;
		min=X_SIZE*X_SIZE*10;
		int fp,got;
		for (fp=0;fp<AP;fp++)
		{
			if (visit[fp]){continue;}
			float r;
			dx=f->x[m]-a->x[fp]; dy=f->y[m]-a->y[fp]; dz=f->z[m]-a->z[fp];
			r=(dx*dx)+(dy*dy)+(dz*dz);
			if (r<min){min=r;got=fp;}
		}
		point[m]=got;
		visit[got]=1;
	}
	free(visit);
}

void show_flock (unsigned char *image, struct floc *f, int w)
{
	int m,oo[2];
	int *ind;
	float pah;

	ind=(int *)malloc(sizeof(int)*MEM);
	for (m=0;m<MEM;m++)
	{
		ind[m]=m;
	}


	pah=M_PI*(float)w*100;


	// first thing lets bubble sort the z values.
	int sorting;
	sorting=1;
	while (sorting)
	{
		sorting=0;
		for (m=0;m<MEM-1;m++)
		{
			int temp;
			if (f->z[ind[m]] > f->z[ind[m+1]])
			{
				sorting=1;
				temp=ind[m];
				ind[m]=ind[m+1];
				ind[m+1]=temp;
			}
		}
	}
	//printf ("min %f max %f \n",f->z[ind[0]],f->z[ind[MEM-1]]);

	

	for (m=0;m<MEM;m++)
	{
		float d,phi,ddy,pha,x,y,z,dx,dy,v;
		x=f->x[ind[m]];
		y=f->y[ind[m]];
		z=f->z[ind[m]];
		v=f->v[ind[m]];
		dx=f->dx[ind[m]];
		dy=f->dy[ind[m]];

		d=3+(9*(z/X_SIZE));
		oo[0]=x; oo[1]=y;
		ddy=dy;
		if (ddy<0){ddy=-ddy;}
		pha=pah*v;
		phi=atan(dx/ddy);
		phi=-phi*180/M_PI;
		linet( image, oo, z*255/X_SIZE, 255, 255-(255*z/X_SIZE), d, phi+20, d/10, 0.5, d/2, pha);
		oo[0]=x; oo[1]=y;
		linet( image, oo, z*255/X_SIZE, 255*y/Y_SIZE, 255-(255*z/X_SIZE), d, phi-20, d/10, 0.5, d/2, -M_PI-pha);
	}

	free (ind);
}


int move_flock (unsigned char *image, struct floc *f, struct ap *a, int *ap)
{
	int m,n,inc;
	float xw,yw,zw,xfac,yfac,zfac;
	inc=0;

	for (m=0;m<MEM;m++)
	{
		xw=a->x[ap[m]];
		yw=a->y[ap[m]];
		zw=a->z[ap[m]]; 
		float min;
		float v,dx,dy,dz,t,rx,ry,rz,r;
		v=f->v[m];
		dx=xw-f->x[m];
		dy=yw-f->y[m];
		dz=zw-f->z[m];
		f->dx[m]=dx;f->dy[m]=dy;
		r=sqrt((dx*dx)+(dy*dy)+(dz*dz));
		if (r<10){ inc++;}
		xfac=1; yfac=1; zfac=1;
		min=100;

		for (n=0;n<MEM;n++)
		{
			float tdis,ddx,ddy,ddz;
			if (m==n){ continue;}
			ddx=f->x[m]-f->x[n];
			ddy=f->y[m]-f->y[n];
			ddz=f->z[m]-f->z[n];
			tdis=(ddx*ddx)+(ddy*ddy)+(ddz*ddz);
			if (tdis<min ){ 
				float hisr;
				ddx=xw-f->x[n];
				ddy=yw-f->y[n];
				ddz=zw-f->z[n];
				hisr=sqrt((ddx*ddx)+(ddy*ddy)+(ddz*ddz));
				if (hisr<r)
				{
					min=tdis;
					/*xfac=1; yfac=1; zfac=1;
					if (dx>dy && dx>dz){ xfac=tdis/100;}
					if (dy>dx && dy>dz){ yfac=tdis/100;}
					if (dz>dy && dz>dx){ zfac=tdis/100;} */
					xfac=tdis/100; yfac=tdis/100; zfac=tdis/100;
				}
			}
			
		}

		f->x[m]+=xfac*v*dx/r;
		f->y[m]+=yfac*v*dy/r;
		f->z[m]+=zfac*v*dz/r;
	}
	//if (inc==1 ){ ap[0]++; printf("Moved %d\n",ap[0]);}
	return inc;
}



void ripple (unsigned char *painted, int depth, int freq)
{
        unsigned char *buff;
        buff=(unsigned char *)malloc(sizeof (char)*3*X_SIZE*Y_SIZE); // buffer
	int x,y;

	for (y=0;y<Y_SIZE;y++)
	{
		int pp,bp;
		pp=y*3*X_SIZE;
		for (x=0;x<X_SIZE;x++)
		{
			int yn;
			yn=y+(depth*sin(2*M_PI*(float)freq*x/(X_SIZE*100)));
			if (yn>0 && yn<Y_SIZE)
			{
				buff[pp+(x*3)]=painted[(yn*X_SIZE*3)+(x*3)];
				buff[pp+(x*3)+1]=painted[(yn*X_SIZE*3)+(x*3)+1];
				buff[pp+(x*3)+2]=painted[(yn*X_SIZE*3)+(x*3)+2];
			}else{
				buff[pp+(x*3)]=255;
				buff[pp+(x*3)+1]=255;
				buff[pp+(x*3)+2]=255;
			}

		}
	}
	memcpy(painted,buff,X_SIZE*Y_SIZE*3);
	free( buff);
}



void gaus(unsigned char *painted, int depth)
{
	float phase,freq;
	int x,y;


	for (y=0;y<Y_SIZE;y++)
	{
		float rbf,rbp;
		float gbf,gbp;
		float bbf,bbp;
		float rmin; rmin=1000000000;
		float gmin; gmin=1000000000;
		float bmin; bmin=1000000000;
		rmin=0;
		gmin=0;
		bmin=0;
		int yp;
		yp=y*3*X_SIZE;
		for (freq=M_PI*2/(X_SIZE);freq<M_PI*8/(X_SIZE);freq+=(M_PI/(depth*X_SIZE)))
		{
			for (phase=0;phase<2*M_PI;phase+=M_PI/(float)depth)
			{
				float rtot; rtot=0;
				float gtot; gtot=0;
				float btot; btot=0;
				int p;
				for (x=0;x<X_SIZE;x++) {
					float guess,col;
					p=yp+(x*3);
					/*
					guess=127*(1+(sin(phase+((float)x*freq))));
					col=guess-painted[p]; rtot+=(col*col);
					col=guess-painted[p+1]; gtot+=(col*col);
					col=guess-painted[p+2]; btot+=(col*col);
					if (rtot>rmin && gtot>gmin && btot>bmin){ break;}
					*/
					guess=(sin(phase+((float)x*freq)));
					rtot+=(float)painted[p]*guess;
					gtot+=(float)painted[p+1]*guess;
					btot+=(float)painted[p+2]*guess;
					if (rtot>rmin && gtot>gmin && btot>bmin){ break;}
				}
				if (rtot>rmin){ rmin=rtot; rbf=freq;rbp=phase;}
				if (gtot>gmin){ gmin=gtot; gbf=freq;gbp=phase;}
				if (btot>bmin){ bmin=btot; bbf=freq;bbp=phase;}
			}
		}
		//printf ("R Phase %f Freq %f y %d \n",rbp,rbf,y);
		//printf ("G Phase %f Freq %f y %d \n",gbp,gbf,y);
		printf ("B Phase %f Freq %f y %d de %d\n",bbp,bbf*X_SIZE/(2*M_PI),y,depth);
		for (x=0;x<X_SIZE;x++) {
			painted[yp+(x*3)]=127*(1+(sin(rbp+((float)x*rbf))));
			painted[yp+(x*3)+1]=127*(1+(sin(gbp+((float)x*gbf))));
			painted[yp+(x*3)+2]=127*(1+(sin(bbp+((float)x*bbf))));
		}
	}
}	



void pline(unsigned char *painted, int angle)
{
	int x,y,along;
	int *xc,*yc;
	xc=(int *)malloc(sizeof(int)*X_SIZE*2);
	yc=(int *)malloc(sizeof(int)*X_SIZE*2);
	for (x=0;x<X_SIZE;x++)
	{
		int len;
		len=gline(xc,yc,x,0,angle);
		for (along=0;along<len;along++) { plott(painted,xc[along],yc[along],along,angle,x); }
		len=gline(xc,yc,x,Y_SIZE-1,angle);
		for (along=0;along<len;along++) { plott(painted,xc[along],yc[along],along,angle,x); }
	}
	for (y=0;y<Y_SIZE;y++)
	{
		int len;
		len=gline(xc,yc,0,y,angle);
		for (along=0;along<len;along++) { plott(painted,xc[along],yc[along],along,angle,y); }
		len=gline(xc,yc,X_SIZE-1,y,angle);
		for (along=0;along<len;along++) { plott(painted,xc[along],yc[along],along,angle,y); }
	}
	free (xc);
	free (yc);
}


int gline(int *xc, int *yc, int x, int y, int angle)
{
	float dx,dy,along;
	dx=sin((float)angle*2*M_PI/360);
	dy=cos((float)angle*2*M_PI/360);
	for (along=0;along<2*X_SIZE;along++)
	{
		int xp,yp;
		float wx,wy;
		xp=x+(along*dx); yp=y+(along*dy);
		if (xp>=X_SIZE || yp>=Y_SIZE || xp<0 || yp<0){break;}
		xc[(int)along]=xp;
		yc[(int)along]=yp;
	}
	return along; 
}


void apaint(unsigned char *painted, int angle)
{
	unsigned char *buff;
	int *escore,*counter;
        int *xc,*yc;
        xc=(int *)malloc(sizeof(int)*X_SIZE*2);
        yc=(int *)malloc(sizeof(int)*X_SIZE*2);
	int depth; depth=5;

	buff=(unsigned char *)malloc(sizeof(unsigned char)*X_SIZE*Y_SIZE*3);
	escore=(int *)malloc(sizeof(int)*X_SIZE*Y_SIZE);
	counter=(int *)malloc(sizeof(int)*3*3*255);

	int x,y,total,score,along;

	if (depth<1){depth=1;}

	for (x=0;x<(3*3*255);x++){ counter[x]=0;}

	for (y=1;y<Y_SIZE;y++)
	{
		int yp,ym;
		yp=y*3*X_SIZE;
		ym=(y-1)*3*X_SIZE;
		for (x=1;x<X_SIZE;x++)
		{
			int pix,pp,pm;
			score=0;
			pp=yp+(x*3);
			pm=ym+(x*3);
			for (pix=0;pix<3;pix++)
			{
				int a;
				a=painted[pp+pix]-painted[pp+pix-3]; //left
				score+=(a*a);
				a=painted[pp+pix]-painted[pm+pix]; //up
				score+=(a*a);
				a=painted[pp+pix]-painted[pm+pix-3]; //up left
				score+=(a*a);
			}
			score=sqrt(score);
			escore[x+(y*X_SIZE)]=score;
			counter[score]++;
		}
	}
	total=0;
	x=(3*3*255)-1;
	while (total<(X_SIZE*Y_SIZE*depth/100) && x>=0)
	{
		total+=counter[x];
		score=x;
		x--;	
	}
	score ++;
	printf ("paint %d %d %d\n",x,total,X_SIZE*Y_SIZE/20);

	int xo;
        for (xo=1;xo<X_SIZE-1;xo++)
        {
                int len;
                len=gline(xc,yc,xo,1,angle);
                for (along=0;along<len;along++)
		{
	       		x=xc[along]; y=yc[along];
			//printf ("%d %d \n",x,y);
			int yp,ep,nr,ng,nb;
			yp=y*3*X_SIZE;
			ep=x+(y*X_SIZE);
			if (escore[ep]>=score || along==0){ 
		 		nr=painted[yp+(x*3)];ng=painted[yp+1+(x*3)];nb=painted[(x*3)+yp+2];
			}
		 	buff[yp+(x*3)]=nr;buff[yp+1+(x*3)]=ng;buff[(x*3)+yp+2]=nb;
		}
	}

        for (xo=1;xo<X_SIZE-1;xo++)
        {
                int len;
                len=gline(xc,yc,xo,Y_SIZE-2,angle);
                for (along=0;along<len;along++)
                {
                        x=xc[along]; y=yc[along];
                        //printf ("%d %d \n",x,y);
                        int yp,ep,nr,ng,nb;
                        yp=y*3*X_SIZE;
                        ep=x+(y*X_SIZE);
                        if (escore[ep]>=score || along==0){
                                nr=painted[yp+(x*3)];ng=painted[yp+1+(x*3)];nb=painted[(x*3)+yp+2];
                        }
                        buff[yp+(x*3)]=nr;buff[yp+1+(x*3)]=ng;buff[(x*3)+yp+2]=nb;
                }
        }

        for (xo=1;xo<Y_SIZE-1;xo++)
        {
                int len;
                len=gline(xc,yc,1,xo,angle);
                for (along=0;along<len;along++)
		{
	       		x=xc[along]; y=yc[along];
			int yp,ep,nr,ng,nb;
			yp=y*3*X_SIZE;
			ep=x+(y*X_SIZE);
			if (escore[ep]>=score || along==0){ 
		 		nr=painted[yp+(x*3)];ng=painted[yp+1+(x*3)];nb=painted[(x*3)+yp+2];
			}
		 	buff[yp+(x*3)]=nr;buff[yp+1+(x*3)]=ng;buff[(x*3)+yp+2]=nb;
		}
	}


        for (xo=1;xo<Y_SIZE-1;xo++)
        {
                int len;
                len=gline(xc,yc,X_SIZE-2,xo,angle);
                for (along=0;along<len;along++)
                {
                        x=xc[along]; y=yc[along];
                        int yp,ep,nr,ng,nb;
                        yp=y*3*X_SIZE;
                        ep=x+(y*X_SIZE);
                        if (escore[ep]>=score || along==0){
                                nr=painted[yp+(x*3)];ng=painted[yp+1+(x*3)];nb=painted[(x*3)+yp+2];
                        }
                        buff[yp+(x*3)]=nr;buff[yp+1+(x*3)]=ng;buff[(x*3)+yp+2]=nb;
                }
        }


	memcpy(painted,buff,X_SIZE*Y_SIZE*3);
	free (escore);
	free (counter);
	free (buff);
        free (xc);
        free (yc);
}

void paint(unsigned char *painted, int depth)
{
	unsigned char *buffer;
	long *escore;
	int *counter;

	buffer=(unsigned char *)malloc(sizeof(unsigned char)*X_SIZE*Y_SIZE*3);
	escore=(long *)malloc(sizeof(long)*X_SIZE*Y_SIZE);
	counter=(int *)malloc(sizeof(int)*3*3*255);

	int x,y,total;
	long score;

	if (depth<1){depth=1;}

	for (x=0;x<(3*3*255);x++){ counter[x]=0;}

	for (y=1;y<Y_SIZE;y++)
	{
		int yp,ym;
		yp=y*3*X_SIZE;
		ym=(y-1)*3*X_SIZE;
		for (x=1;x<X_SIZE;x++)
		{
			int pix,pp,pm;
			score=0;
			pp=yp+(x*3);
			pm=ym+(x*3);
			for (pix=0;pix<3;pix++)
			{
				int a;
				a=painted[pp+pix]-painted[pp+pix-3]; //left
				score+=(a*a);
				a=painted[pp+pix]-painted[pm+pix]; //up
				score+=(a*a);
				a=painted[pp+pix]-painted[pm+pix-3]; //up left
				score+=(a*a);
			}
			score=sqrt(score);
			escore[x+(y*X_SIZE)]=score;
			counter[score]++;
		}
	}
	total=0;
	x=(3*3*255)-1;
	while (total<(X_SIZE*Y_SIZE*depth/100) && x>=0)
	{
		total+=counter[x];
		score=x;
		x--;	
	}
	score ++;
	printf ("paint %d %d %d\n",x,total,X_SIZE*Y_SIZE/20);

	for (y=1;y<Y_SIZE;y++)
	{
		int yp,ep,xb,xe;
		int edgy,nr,ng,nb,along;
		yp=y*3*X_SIZE;
		xb=0;xe=0;
		edgy=1;
		for (x=1;x<X_SIZE;x++)
		{
			ep=x+(y*X_SIZE);
			if (escore[ep]>=score ){ 
				xb=xe;xe=x;
				nr=0;ng=0;nb=0;
				for (along=xb;along<xe;along++)
				{
		 			nr+=painted[yp+(along*3)];ng+=painted[(along*3)+yp+1];nb+=painted[(along*3)+yp+2];
				}	
				nr=nr/(xe-xb); ng=ng/(xe-xb); nb=nb/(xe-xb);
				for (along=xb;along<xe;along++)
				{
		  			buffer[yp+(along*3)]=nr ;buffer[(along*3)+yp+1]=ng;buffer[(along*3)+yp+2]=nb;
				}	
		  			buffer[yp+(xb*3)]=nr/2 ;buffer[(xb*3)+yp+1]=ng/2;buffer[(xb*3)+yp+2]=nb/2;
			}
		}
	}
	memcpy(painted,buffer,X_SIZE*Y_SIZE*3);
/*
        for (x=0;x<(3*3*255);x++){ counter[x]=0;}
        for (x=1;x<X_SIZE;x++)
        {
                for (y=1;y<Y_SIZE;y++)
                {
                	int yp,ym;
                	yp=y*3*X_SIZE;
                	ym=(y-1)*3*X_SIZE;
                        int pix,pp,pm;
                        score=0;
                        pp=yp+(x*3);
                        pm=ym+(x*3);
                        for (pix=0;pix<3;pix++)
                        {
                                int a;
                                //a=painted[pp+pix]-painted[pp+pix-3]; //left
                                //score+=(a*a);
                                a=painted[pp+pix]-painted[pm+pix]; //up
                                score+=(a*a);
                                //a=painted[pp+pix]-painted[pm+pix-3]; //up left
                                //score+=(a*a);
                        }
                        score=sqrt(score);
                        escore[x+(y*X_SIZE)]=score;
                        counter[score]++;
                }
        }

        total=0;
        x=(3*3*255)-1;
        while (total<(X_SIZE*Y_SIZE*depth/100) && x>=0)
        {
                total+=counter[x];
                score=x;
                x--;
        }
        score ++;
        printf ("paint %d %d %d\n",x,total,X_SIZE*Y_SIZE/30);


        for (x=1;x<X_SIZE;x++)
        {
                int yp,ep,nr,ng,nb;
		int xb,xe,along;
		xb=0;xe=0;
        	for (y=1;y<Y_SIZE;y++)
                {
                        ep=x+(y*X_SIZE);
                        if (escore[ep]>=score ){
                                xb=xe;xe=y;
                                nr=0;ng=0;nb=0;
                                for (along=xb;along<xe;along++)
                                {
                			yp=along*3*X_SIZE;
                                        nr+=painted[yp+(x*3)];ng+=painted[(x*3)+yp+1];nb+=painted[(x*3)+yp+2];
                                }
                                nr=nr/(xe-xb); ng=ng/(xe-xb); nb=nb/(xe-xb);
                                for (along=xb;along<xe;along++)
                                {
                			yp=along*3*X_SIZE;
                                        buffer[yp+(x*3)]=(buffer[yp+(x*3)]+nr)/2;
                                        buffer[yp+(x*3)+1]=(buffer[yp+(x*3)+1]+ng)/2;
                                        buffer[yp+(x*3)+2]=(buffer[yp+(x*3)+2]+nb)/2;

                                }
                        }

                }
        }
	memcpy(painted,buffer,X_SIZE*Y_SIZE*3);*/

	free (escore);
	free (counter);
	free (buffer);
}


void blurt (unsigned char *image2, int depth)
{
	unsigned char *image1;
        image1=(unsigned char *)malloc(sizeof (char)*3*X_SIZE*Y_SIZE); // blur buffer
	int x,y,d,pix,ly,xe;
	ly=STRIDE*(Y_SIZE-2);
	xe=STRIDE-3;

	for (d=0;d<depth;d++)
	{

		// corners
		for (pix=0;pix<3;pix++){ 
			image1[pix]=(image2[pix]+image2[pix+3]+image2[pix+3+STRIDE]+image2[pix+STRIDE])/4;
			image1[xe+pix]=(image2[xe+pix]+image2[xe-3+pix]+image2[xe+pix-3+STRIDE]+image2[xe+pix+STRIDE])/4;
			image1[pix+ly]=(image2[pix+ly]+image2[ly+pix+3]+image2[ly+pix+3-STRIDE]+image2[ly+pix-STRIDE])/4;
			image1[xe+pix+ly]=(image2[xe+pix+ly]+image2[xe-3+pix+ly]+image2[xe+pix-3+STRIDE+ly]+image2[xe+pix+STRIDE+ly])/4;
		}
		// top/bottom
		for (x=3;x<STRIDE-4;x+=3)
		{
			for (pix=0;pix<3;pix++){ 
				image1[pix+x]=(48+image2[pix+x]+image2[pix+3+x]+image2[pix+x-3]+
						image2[pix+x+STRIDE]+image2[pix+3+x+STRIDE]+image2[pix+x-3+STRIDE])/6;
				image1[pix+x+ly]=(48+image2[pix+x+ly]+image2[pix+3+x+ly]+image2[pix+x-3+ly]+
						image2[pix+x+ly-STRIDE]+image2[pix+3+x+ly-STRIDE]+image2[pix+x-3+ly-STRIDE])/6;
			}
		}
		// left/right
		for (y=1;y<Y_SIZE-1;y++)
		{
			int yp,yr;
			yp=y*STRIDE;
			yr=yp+STRIDE-3;
			for (pix=0;pix<3;pix++){ 
				image1[yp+pix]=(48+image2[yp+pix]+image2[yp+pix-STRIDE]+image2[yp+pix+STRIDE]+
					image2[yp+pix+3]+image2[yp+pix-STRIDE+3]+image2[yp+pix+STRIDE+3])/6;
				image1[yr+pix]=(48+image2[yr+pix]+image2[yr+pix-STRIDE]+image2[yr+pix+STRIDE]+
					image2[yr+pix-3]+image2[yr+pix-STRIDE-3]+image2[yr+pix+STRIDE-3])/6;
			}
		}

		for (y=1;y<Y_SIZE-1;y++)
		{
			for (x=3;x<STRIDE-4;x+=3)
			{
				int u,p,d;
				p=(x)+(y*STRIDE);	
				u=p+STRIDE;	
				d=p-STRIDE;	
				for (pix=0;pix<3;pix++)
				{
					//int n,m;
					//if (rand()%2000<500){ n=rand()%3;if (pix+n>2) {n-=2;}}
					//n=rand()%3;if (pix+n>2) {n-=2;}
					//m=rand()%3;if (pix+m>2) {m-=2;}
					image1[p+pix]=(image2[p-3+pix]+image2[p+pix]+image2[p+3+pix]+
						image2[u-3+pix]+image2[u+pix]+image2[u+3+pix]+
						image2[d-3+pix]+image2[d+pix]+image2[d+3+pix]
					  )/9;
				}
			}
		}

		// corners
		for (pix=0;pix<3;pix++){ 
			image2[pix]=(image1[pix]+image1[pix+3]+image1[pix+3+STRIDE]+image1[pix+STRIDE])/4;
			image2[xe+pix]=(image1[xe+pix]+image1[xe-3+pix]+image1[xe+pix-3+STRIDE]+image1[xe+pix+STRIDE])/4;
			image2[pix+ly]=(image1[pix+ly]+image1[ly+pix+3]+image1[ly+pix+3-STRIDE]+image1[ly+pix-STRIDE])/4;
			image2[xe+pix+ly]=(image1[xe+pix+ly]+image1[xe-3+pix+ly]+image1[xe+pix-3+STRIDE+ly]+image1[xe+pix+STRIDE+ly])/4;
		}
		// top/bottom
		for (x=3;x<STRIDE-12;x+=3)
		{
			for (pix=0;pix<3;pix++){ 
				image2[pix+x]=(image1[pix+x]+image1[pix+3+x]+image1[pix+x-3]+
						image1[pix+x+STRIDE]+image1[pix+3+x+STRIDE]+image1[pix+x-3+STRIDE])/6;
				image2[pix+x+ly]=(image1[pix+x+ly]+image1[pix+3+x+ly]+image1[pix+x-3+ly]+
						image1[pix+x+ly-STRIDE]+image1[pix+3+x+ly-STRIDE]+image1[pix+x-3+ly-STRIDE])/6;
			}
		}

		// left/right
		for (y=1;y<Y_SIZE-1;y++)
		{
			int yp,yr;
			yp=y*STRIDE;
			yr=yp+STRIDE-3;
			for (pix=0;pix<3;pix++){ 
				image2[yp+pix]=(image1[yp+pix]+image1[yp+pix-STRIDE]+image1[yp+pix+STRIDE]+
					image1[yp+pix+3]+image1[yp+pix-STRIDE+3]+image1[yp+pix+STRIDE+3])/6;
				image2[yr+pix]=(image1[yr+pix]+image1[yr+pix-STRIDE]+image1[yr+pix+STRIDE]+
					image1[yr+pix-3]+image1[yr+pix-STRIDE-3]+image1[yr+pix+STRIDE-3])/6;
			}
		}

                for (y=1;y<Y_SIZE-1;y++)
                {
                        for (x=3;x<STRIDE-4;x+=3)
                        {
                                int u,p,d;
                                p=(x)+(y*STRIDE);
                                u=p+STRIDE;
                                d=p-STRIDE;
                                for (pix=0;pix<3;pix++)
                                {
                                        image2[p+pix]=(image1[p-3+pix]+image1[p+pix]+image1[p+3+pix]+
                                                image1[u-3+pix]+image1[u+pix]+image1[u+3+pix]+
                                                image1[d-3+pix]+image1[d+pix]+image1[d+3+pix]
                                          )/9;
                                }
                        }
                }

	}
	free (image1);
}


void eartht(unsigned char *image2, int stage)
{
        int x,y;
	int oo[2];
	int atten;

	oo[0]=0;
	oo[1]=0;


        //square (0,0,90,90,255,X_SIZE,Y_SIZ2); //sky
        //square (0,Y_SIZ2,120,80,10,X_SIZE,Y_SIZ2); //earth
	//sky

	srand(0);
	float mod,moe;
	if (stage<120){ mod=(float)stage/120;}else{mod=1;}
	if (stage<240) { for (x=0;x<X_SIZE*Y_SIZE*3;x++){image2[x]=(rand()%256);} blurt(image2,10);} // background
	srand(0);
	for (y=0;y<mod*Y_SIZ2;y+=4)
	{
		oo[0]=0;
		oo[1]=y;
		atten=100+((155*(Y_SIZ2-y))/Y_SIZE);
		linet(image2,oo,atten*(110+(mod*(rand()%30)))/255,atten*(130+(mod*(rand()%50)))/255,atten,X_SIZE,90,10*mod,mod*((float)(rand()%30)),3,(float)(rand()%100)*mod/10);
	}		
	
	srand(0);
	if (stage>=120) {
		if (stage<240){ moe=((float)stage-120)/120;}else{moe=1;}
		for (y=Y_SIZ2;y<Y_SIZ2+(moe*Y_SIZ2);y+=4)
		{
			oo[0]=0;
			oo[1]=y;
			atten=55+((200*(y+Y_SIZ2))/Y_SIZE);
			linet(image2,oo,
					atten, // red
					atten*(80+(moe*(rand()%(50))))/255, //green
					atten*((moe*(rand()%50)))/255,X_SIZE, //blue
					90, //ang
					1+(moe*((y-Y_SIZ2))/10), // thick
					20+((float)(Y_SIZE-y)/20), // wiggles
					moe*8*((float)y-Y_SIZ2)/Y_SIZ2,  // wiggle depth
					2*M_PI*(float)(moe*(rand()%Y_SIZE))/Y_SIZE  );
		}
	}
        if (stage>240 && stage<300){blurt(image2,(stage-240)/4); }
	if (stage>=300){blurt(image2,15);}
}


void linet( unsigned char *image2, int *origin, int r, int g, int b, float length, int angle, float thick, float wig, float depth, float pha)
{
        float phi,dx,dy,along,tx,ty,wide;
        phi=(2*M_PI*(float)angle/360);

        dx=(float)length*sin(phi);
        dy=-(float)length*cos(phi);
        tx=(float)thick*sin(phi+(M_PI/2));
        ty=-(float)thick*cos(phi+(M_PI/2));

        for (along=0;along<length;along+=0.6)
        {
		float off;
		off=depth*(sin(pha+(along*wig*M_PI/length)))/thick;
                for (wide=-thick/2;wide<thick/2;wide+=0.6)
                {
                        float xp,yp;
                        xp=((along*dx)/length)+((tx*wide)/thick)+(off*tx);
                        yp=((along*dy)/length)+((ty*wide)/thick)+(off*ty);
                        plott(image2,origin[0]+xp,origin[1]+yp,r,g,b);
                }
        }
        origin[0]+=dx;
        origin[1]+=dy;
}

void lineu( unsigned char *image2, float *beg, float *end ,int r, int g, int b, float thick)
{
        float dx,dy,along,tx,ty,wide,len;

        dx=end[0]-beg[0];
        dy=end[1]-beg[1];
	len=sqrt((dx*dx)+(dy*dy));
        tx=-dy/len;
        ty=dx/len;

	printf("%f %f %f\n",thick,tx,ty);

        for (along=0;along<len;along+=0.6)
        {
                for (wide=-thick/2;wide<thick/2;wide+=0.6)
                {
                        float xp,yp;
                        xp=((along*dx)/len)+((tx*wide));
                        yp=((along*dy)/len)+((ty*wide));
                        plott(image2,beg[0]+xp,beg[1]+yp,r,g,b);
                }
        }
}

// plot into the final
void plott (unsigned char *image2, int x,int y,int r,int g,int b)
{
        if (x>=X_SIZE){ return;}
        if (x<0){ return;}
        if (y>=Y_SIZE){ return;}
        if (y<0){ return;}
        int xpoint;
        xpoint=(y*X_SIZE*3)+(x*3);
        image2[xpoint]=r;
        image2[xpoint+1]=g;
        image2[xpoint+2]=b;
}

void merge(unsigned char *out, unsigned char *in, int amount)
{
        if (amount<0){ amount=0;}
        if (amount>255){ amount=255;}
        int prime,x;
        prime=255-amount;
        for (x=0;x<X_SIZE*Y_SIZE*3;x++)
        {
                out[x]=((amount*out[x])+(prime*in[x]))/256;
        }
}

