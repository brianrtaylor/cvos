#include "qx_basic.h"
#include "qx_ppm.h"
#include "qx_recursive_bilateral_filter.h"
int example_filtering(int argc,char*argv[]);
int example_detail_enhancement(int argc,char*argv[]);
int main(int argc,char*argv[])
{
	if(argc!=6)
	{
		printf("Usage:\n");
		printf("--------------------------------------------------------------------\n\n");
		printf("*.exe: filename_out filename_in (only ppm image) \n");
		printf("       sigma_spatial(e.g., 0.03) sigma_range(e.g., 0.1)\n");
		printf("       application_id (0 for GDBF; 1 for RBF and 2 for detail enhancement)\n\n");
		printf("--------------------------------------------------------------------\n");
		return(-1);
	}
	else
	{
		if(atoi(argv[5])==0||atoi(argv[5])==1) return(example_filtering(argc,argv));
		else if(atoi(argv[5])==2) return(example_detail_enhancement(argc,argv));
	}
}
int example_filtering(int argc,char*argv[])
{
	char*filename_out=argv[1];
	char*filename_in=argv[2];
	float sigma_spatial=atof(argv[3]);
	float sigma_range=atof(argv[4]);
	int filter_type=atoi(argv[5]);

	int h,w;//height and width of the input image
	qx_image_size(filename_in,h,w);//obtain the height and width of the input image
	unsigned char***texture=qx_allocu_3(h,w,3);//allocate memory
	double***image=qx_allocd_3(h,w,3);
	double***image_filtered=qx_allocd_3(h,w,3);
	double***temp=qx_allocd_3(h,w,3);
	double***temp_2=qx_allocd_3(2,w,3);

	qx_loadimage(filename_in,texture[0][0],h,w);//read an color image
	for(int y=0;y<h;y++) for(int x=0;x<w;x++) for(int c=0;c<3;c++) image[y][x][c]=texture[y][x][c];//initialize the original color image

	qx_timer timer;//run filter: qx_gradient_domain_recursive_bilateral_filter
	timer.start();
	int nr_iteration=10;
	if(filter_type==0)//GDBF
	{
		for(int i=0;i<nr_iteration;i++) qx_gradient_domain_recursive_bilateral_filter(image_filtered,image,texture,sigma_spatial,sigma_range,h,w,temp,temp_2);
	}
	else//RBF
	{
		double**temp_factor=qx_allocd(h*2+2,w);
		for(int i=0;i<nr_iteration;i++) qx_recursive_bilateral_filter(image_filtered,image,texture,sigma_spatial,sigma_range,h,w,temp,temp_2,temp_factor,&(temp_factor[h]),&(temp_factor[h+h]));
		qx_freed(temp_factor); temp_factor=NULL;
	}
	timer.fps_display("recursive bilateral filter",nr_iteration);

	for(int y=0;y<h;y++) for(int x=0;x<w;x++) for(int c=0;c<3;c++) texture[y][x][c]=image_filtered[y][x][c]; //write the filtered image to the hard drive
	qx_saveimage(filename_out,texture[0][0],h,w,3);
	
	qx_freeu_3(texture);
	qx_freed_3(image);
	qx_freed_3(image_filtered);
	qx_freed_3(temp);
	qx_freed_3(temp_2);
	return(0);
}
int example_detail_enhancement(int argc,char*argv[])
{
	char*filename_out=argv[1];
	char*filename_in=argv[2];
	float sigma_spatial=atof(argv[3]);
	float sigma_range=atof(argv[4]);

	int h,w;//height and width of the input image
	qx_image_size(filename_in,h,w);//obtain the height and width of the input image
	unsigned char***texture=qx_allocu_3(h,w,3);//allocate memory
	double***image=qx_allocd_3(h,w,3);
	double***image_filtered=qx_allocd_3(h,w,3);
	double***temp=qx_allocd_3(h,w,3);
	double***temp_2=qx_allocd_3(2,w,3);

	qx_loadimage(filename_in,texture[0][0],h,w);//read an color image
	for(int y=0;y<h;y++) for(int x=0;x<w;x++) for(int c=0;c<3;c++) image[y][x][c]=texture[y][x][c];//initialize the original color image

	qx_gradient_domain_recursive_bilateral_filter(image_filtered,image,texture,sigma_spatial,sigma_range,h,w,temp,temp_2);//filtering

	for(int y=0;y<h;y++) for(int x=0;x<w;x++) for(int c=0;c<3;c++) //detail enhancement
	{
		texture[y][x][c]=(unsigned char)min(255.0,max(0.0,image_filtered[y][x][c]+5.0*(image[y][x][c]-image_filtered[y][x][c])));
	}
	qx_saveimage(filename_out,texture[0][0],h,w,3);//write the filtered image to the hard drive

	qx_freeu_3(texture);
	qx_freed_3(image);
	qx_freed_3(image_filtered);
	qx_freed_3(temp);
	qx_freed_3(temp_2);
	return(0);
}