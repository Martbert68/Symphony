#include <time.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <string.h>
#include <jerror.h>
#include <jpeglib.h>
#include <setjmp.h>
#include "martin.h"

/* here are our X variables */
Display *dis;
int screen;
Window win;
GC gc;
XImage *x_image;
unsigned char *x_buffer;

/* here are our X routines declared! */
void init_x();
void close_x();
void redraw();


void disp (unsigned char *,int,int);


/*void alloc_try (struct try *tri,int count)
{
	tri=(struct try *)malloc(sizeof (struct try)*count);
	int loop;
	for (loop=0;loop<count;loop++)
	{
		tri->pi=(float *)malloc(sizeof(float)*3);
		tri->pj=(float *)malloc(sizeof(float)*3);
		tri->pk=(float *)malloc(sizeof(float)*3);
	}
}*/

void alloc_try (struct try *tri,struct trf *flt,int count)
{
	int loop;
	for (loop=0;loop<count;loop++)
	{
		tri[loop].pi=(float *)malloc(sizeof(float)*3);
		tri[loop].pj=(float *)malloc(sizeof(float)*3);
		tri[loop].pk=(float *)malloc(sizeof(float)*3);
		flt[loop].pi=(float *)malloc(sizeof(float)*2);
		flt[loop].pj=(float *)malloc(sizeof(float)*2);
		flt[loop].pk=(float *)malloc(sizeof(float)*2);
	}
}

void usage ()
{
	printf("usage: font filename threshold [20-40 ish] star,framestcode [65A 97a]\n");
	exit (1);
}


int main(int argc,char *argv[])
{
	unsigned char *image2,*image9;
	struct floc *f;
	struct ap *a;
	struct try *t; 
	struct trf *flat; 
	float vel;
	int along;
	int jay,*order;
	int *point;
	FILE *stl;
        image2=(unsigned char *)malloc(sizeof (char)*3*X_SIZE*Y_SIZE); // disp buffer
        image9=(unsigned char *)malloc(sizeof (char)*3*X_SIZE*Y_SIZE); // disp buffer
        f=(struct floc *)malloc(sizeof (struct floc)); // flock 
        a=(struct ap *)malloc(sizeof (struct ap)); // attractors 

        point=(int *)malloc(sizeof (int)*MEM); // xpos 
        f->x=(float *)malloc(sizeof (float)*MEM); // xpos 
        f->dx=(float *)malloc(sizeof (float)*MEM); // xpos 
        f->y=(float *)malloc(sizeof (float)*MEM); // xpos 
        f->dy=(float *)malloc(sizeof (float)*MEM); // xpos 
        f->z=(float *)malloc(sizeof (float)*MEM); // xpos 
        f->v=(float *)malloc(sizeof (float)*MEM); // xpos 

        a->x=(float *)malloc(sizeof (float)*AP); // xpos 
        a->y=(float *)malloc(sizeof (float)*AP); // xpos 
        a->z=(float *)malloc(sizeof (float)*AP); // xpos 

	if (argv[1]){jay=1;}else{jay=0;}
	init_x();

	//stl=fopen("./Pig.stl","r");
	stl=fopen("./jaguar.stl","r");
	if (stl==NULL){exit(1);}

	char dummy[80];
	unsigned int tri_count;
	int read,loop;
	read=fread(dummy,sizeof(char),80,stl); printf("read %d\n",read);
	read=fread(&tri_count,4,1,stl); printf("read %d\n",read);
	printf ("I reckon there are %d trianges\n",tri_count);

	//allocate space for the triangles.
	t=(struct try *)malloc(sizeof (struct try)*tri_count);
	flat=(struct trf *)malloc(sizeof (struct trf)*tri_count);
	order=(int *)malloc(sizeof (int)*tri_count);
	for (along=0;along<tri_count;along++){ order[along]=along;}
	alloc_try(t,flat,tri_count);

	short attrib;
	float normal[3];
	int ll;
	loop=0;
	while ( fread(normal,sizeof(float),3,stl))
	{
		ll=loop/40;
		read=fread(t[ll].pi,sizeof(float),3,stl);
		read=fread(t[ll].pj,sizeof(float),3,stl);
		read=fread(t[ll].pk,sizeof(float),3,stl);
		read=fread(&attrib,sizeof(short),1,stl);
		loop++;
	}
	loop=ll;
	tri_count/=40;

	printf ("I reckon there are %d trianges and I read %d\n",tri_count,loop);
	fclose(stl);	

	// We have loaded the file!

	float xmax,ymax,zmax; float xmin,ymin,zmin;
	xmax=-1000000;ymax=-1000000;zmax=-1000000; xmin=1000000;ymin=1000000;zmin=1000000;
	for (along=0;along<tri_count;along++)
	{
		if (t[along].pi[0]>xmax){xmax=t[along].pi[0];} if (t[along].pj[0]>xmax){xmax=t[along].pj[0];} if (t[along].pk[0]>xmax){xmax=t[along].pk[0];}
		if (t[along].pi[0]<xmin){xmin=t[along].pi[0];} if (t[along].pj[0]<xmin){xmin=t[along].pj[0];} if (t[along].pk[0]<xmin){xmin=t[along].pk[0];}
		if (t[along].pi[1]>ymax){ymax=t[along].pi[1];} if (t[along].pj[1]>ymax){ymax=t[along].pj[1];} if (t[along].pk[1]>ymax){ymax=t[along].pk[1];}
		if (t[along].pi[1]<ymin){ymin=t[along].pi[1];} if (t[along].pj[1]<ymin){ymin=t[along].pj[1];} if (t[along].pk[1]<ymin){ymin=t[along].pk[1];}
		if (t[along].pi[2]>zmax){zmax=t[along].pi[2];} if (t[along].pj[2]>zmax){zmax=t[along].pj[2];} if (t[along].pk[2]>zmax){zmax=t[along].pk[2];}
		if (t[along].pi[2]<zmin){zmin=t[along].pi[2];} if (t[along].pj[2]<zmin){zmin=t[along].pj[2];} if (t[along].pk[2]<zmin){zmin=t[along].pk[2];}
	}
	printf ("xmax %f xmin %f ymax %f ymin %f zmax %f zmin %f \n",xmax,xmin,ymax,ymin,zmax,zmin);

	// as a guess lets scale xmax to 200.
	float scale;
	scale=200/xmax;
	for (along=0;along<tri_count;along++)
	{
		int i;
		for (i=0;i<3;i++)
		{
			t[along].pi[i]*=scale; t[along].pj[i]*=scale; t[along].pk[i]*=scale;
		}
	}

	// let try
	float in[3];
	float out[3];
	float view[3];
	float cm[3];
	float xy[2];

	view[0]=1550;view[1]=-1520;view[2]=-300;

	float pan,tilt,tlt,l;
	tilt=-M_PI/1.5;l=-3000 ; pan=M_PI;
	int frame,red,green;
	float p;
	float cp,sp,ct,st;


	// calculate pan and tilt.
	// pan is based on x/y
	//
	cm[0]=0;cm[1]=-200;cm[2]=0;
	
	//pan=atan((view[1]-cm[1])/(view[0]-cm[0]));
	//tilt=asin((view[2]-cm[2])/sqrt(((view[0]-cm[0])*(view[0]-cm[0]))+((view[1]-cm[1])*(view[1]-cm[1]))+((view[2]-cm[2])*(view[2]-cm[2]))));


	for (p =0;p<2000;p++)
	{
	//view[2]++;
	view[0]-=4;
	view[1]+=4;
	view[2]+=9;
	pan=atan((view[1]-cm[1])/(view[0]-cm[0]));
	tilt=asin((view[2]-cm[2])/sqrt(((view[0]-cm[0])*(view[0]-cm[0]))+((view[1]-cm[1])*(view[1]-cm[1]))+((view[2]-cm[2])*(view[2]-cm[2]))));
	clear(image2);
	matrix(pan,tilt,&cp,&sp,&ct,&st);

	dist_tri(t,view,tri_count);
	sort_tri(t,order,tri_count);

	for (along=0;along<tri_count;along++)
	{
		if (rend_tri(t+order[along], flat,view, cp,sp, ct, st, l))
		{
			draw_triangle (image2, flat,255*along/tri_count,255*along/tri_count,255*along/tri_count);
		}
		//if (along%1000==0){disp(image2,frame,0);}
	}
	disp(image2,p,1);
	}
	char c;
	scanf("%c",&c);
	close_x();


	exit(0);

}	

void disp (unsigned char *image2,int fram,int ab)
{
	int x,y;
	char input[100];


       	for (y=0;y<Y_SIZE;y++)
       	{
               	int p=y*STRIDE;
               	int XYP=X_SIZE*4*y;
               	for (x=0;x<X_SIZE;x++)
               	{
			int xpoint;
			int X_POINT;
			X_POINT=XYP+(4*x);
			xpoint=(x*3)+(p);

			x_buffer[X_POINT+2]=image2[xpoint];
			x_buffer[X_POINT+1]=image2[xpoint+1];
			x_buffer[X_POINT]=image2[xpoint+2];
                }
        }
	XPutImage(dis, win, gc, x_image, 0, 0, 0, 0, X_SIZE, Y_SIZE);
	sprintf(input,"./jpegs/jmage%04d.jpg",fram);
	if (ab){jayit(image2,X_SIZE, Y_SIZE, input);}
}


struct my_error_mgr {
  struct jpeg_error_mgr pub;	/* "public" fields */

  jmp_buf setjmp_buffer;	/* for return to caller */
};

typedef struct my_error_mgr * my_error_ptr;

/*
 * Here's the routine that will replace the standard error_exit method:
 */

METHODDEF(void)
my_error_exit (j_common_ptr cinfo)
{
  /* cinfo->err really points to a my_error_mgr struct, so coerce pointer */
  my_error_ptr myerr = (my_error_ptr) cinfo->err;

  /* Always display the message. */
  /* We could postpone this until after returning, if we chose. */
  (*cinfo->err->output_message) (cinfo);

  /* Return control to the setjmp point */
  longjmp(myerr->setjmp_buffer, 1);
}

GLOBAL(int)
read_JPEG_file (char * filename, unsigned char * dots, int * params)
{
  /* This struct contains the JPEG decompression parameters and pointers to
   * working space (which is allocated as needed by the JPEG library).
   */
  struct jpeg_decompress_struct cinfo;
  /* We use our private extension JPEG error handler.
   * Note that this struct must live as long as the main JPEG parameter
   * struct, to avoid dangling-pointer problems.
   */
  struct my_error_mgr jerr;
  /* More stuff */
  FILE * infile;		/* source file */
  JSAMPARRAY buffer;		/* Output row buffer */
  int row_stride;		/* physical row width in output buffer */

  if ((infile = fopen(filename, "rb")) == NULL) {
    fprintf(stderr, "can't open %s\n", filename);
    return 0;
  }

  /* Step 1: allocate and initialize JPEG decompression object */

  /* We set up the normal JPEG error routines, then override error_exit. */
  cinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = my_error_exit;
  /* Establish the setjmp return context for my_error_exit to use. */
  if (setjmp(jerr.setjmp_buffer)) {
    /* If we get here, the JPEG code has signaled an error.
     * We need to clean up the JPEG object, close the input file, and return.
     */
    jpeg_destroy_decompress(&cinfo);
    fclose(infile);
    return 0;
  }
  /* Now we can initialize the JPEG decompression object. */
  jpeg_create_decompress(&cinfo);

  /* Step 2: specify data source (eg, a file) */

  jpeg_stdio_src(&cinfo, infile);

  /* Step 3: read file parameters with jpeg_read_header() */

  (void) jpeg_read_header(&cinfo, TRUE);
  /* We can ignore the return value from jpeg_read_header since
   *   (a) suspension is not possible with the stdio data source, and
   *   (b) we passed TRUE to reject a tables-only JPEG file as an error.
   * See libjpeg.txt for more info.
   */

  /* Step 5: Start decompressor */

  (void) jpeg_start_decompress(&cinfo);
  /* We can ignore the return value since suspension is not possible
   * with the stdio data source.
   */

  /* We may need to do some setup of our own at this point before reading
   * the data.  After jpeg_start_decompress() we have the correct scaled
   * output image dimensions available, as well as the output colormap
   * if we asked for color quantization.
   * In this example, we need to make an output work buffer of the right size.
   */ 
  /* JSAMPLEs per row in output buffer */
  row_stride = cinfo.output_width * cinfo.output_components;
  /* Make a one-row-high sample array that will go away when done with image */
  buffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);


  /* Step 6: while (scan lines remain to be read) */
  /*           jpeg_read_scanlines(...); */

  /* Here we use the library's state variable cinfo.output_scanline as the
   * loop counter, so that we don't have to keep track ourselves.
   */

  while (cinfo.output_scanline < cinfo.output_height) {
    /* jpeg_read_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could ask for
     * more than one scanline at a time if that's more convenient.
     */
    (void) jpeg_read_scanlines(&cinfo, buffer, 1);
    memcpy (dots+(row_stride*cinfo.output_scanline),buffer[0],row_stride);
    /* Assume put_scanline_someplace wants a pointer and sample count. */
    /* put_scanline_someplace(buffer[0], row_stride); */

  }
  /* Step 7: Finish decompression */
  params[0]=cinfo.output_width;
  params[1]=cinfo.output_height;
  params[2]=cinfo.output_components;

  (void) jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);
  fclose(infile);

  /* And we're done! */
  return 1;
}

int jayit(unsigned char *screen,int image_width, int image_height, char *name)
{

int row_stride,ex,why,cmp,div,set;
unsigned char *image,**row_pointer,*cr,*cg,*cb;
row_pointer=(unsigned char **)malloc(1);

struct jpeg_compress_struct cinfo;
struct jpeg_error_mgr jerr;
FILE * outfile;		/* target file */
cinfo.err = jpeg_std_error(&jerr);
jpeg_create_compress(&cinfo);
if ((outfile = fopen(name, "wb")) == NULL) { 
	fprintf(stderr, "can't open file\n");
	exit(1);
}
jpeg_stdio_dest(&cinfo, outfile);
cinfo.image_width = image_width; 	/* image width and height, in pixels */
cinfo.image_height = image_height;
cinfo.input_components = 3;		/* # of color components per pixel */
cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
jpeg_set_defaults(&cinfo);
jpeg_set_quality(&cinfo,100,TRUE); /* limit to baseline-JPEG values */
jpeg_start_compress(&cinfo, TRUE);

  row_stride = image_width * 3;	/* JSAMPLEs per row in image_buffer */

  while (cinfo.next_scanline < cinfo.image_height) {
    /* jpeg_write_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could pass
     * more than one scanline at a time if that's more convenient.
     */
    row_pointer[0] = & screen[cinfo.next_scanline * row_stride];
    (void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
  }
jpeg_finish_compress(&cinfo);
fclose(outfile);
jpeg_destroy_compress(&cinfo);
}

void init_x()
{
/* get the colors black and white (see section for details) */
        unsigned long black,white;

        x_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        //y_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        //z_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        dis=XOpenDisplay((char *)0);
        screen=DefaultScreen(dis);
        black=BlackPixel(dis,screen),
        white=WhitePixel(dis,screen);
        win=XCreateSimpleWindow(dis,DefaultRootWindow(dis),0,0,
                X_SIZE, Y_SIZE, 5, white,black);
        XSetStandardProperties(dis,win,"image","images",None,NULL,0,NULL);
        gc=XCreateGC(dis, win, 0,0);
        XSetBackground(dis,gc,black); XSetForeground(dis,gc,white);
        XClearWindow(dis, win);
        XMapRaised(dis, win);
        //XMoveWindow(dis, win,window_x,100);
        Visual *visual=DefaultVisual(dis, 0);
        x_image=XCreateImage(dis, visual, DefaultDepth(dis,DefaultScreen(dis)), ZPixmap, 0, x_buffer, X_SIZE, Y_SIZE, 32, 0);
};

void close_x() {
        XFreeGC(dis, gc);
        XDestroyWindow(dis,win);
        XCloseDisplay(dis);
        exit(1);
};

void redraw() {
        XClearWindow(dis, win);
};

