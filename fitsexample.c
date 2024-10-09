#include "fitsio.h" /* required by every program that uses CFITSIO */
#include "stdio.h"

int main()
{
  fitsfile *fptr, *gptr; /* pointer to the FITS file; defined in fitsio.h */
  int status, ii, jj;
  long fpixel=1, naxis = 2, nelements, exposure;
  long naxes[2] = { 300, 200 }; /* image is 300 pixels wide by 200 rows */
  short array[200][300];

  status = 0; /* initialize status before calling fitsio routines */
  fits_open_file(&fptr, "image_out_inc52.9_PA_226.2.fits", READONLY, &status);
  /* Create the primary array image (16-bit short integer pixels */
  //fits_create_img(fptr, SHORT_IMG, naxis, naxes, &status);
  
  fits_get_img_size(fptr, 2, naxes, &status);
  double *myimage = (double *) malloc(naxes[0]*naxes[1]*sizeof(double));
  int nval[1024*1024];
  
  printf("size = (%ld, %ld)\n", naxes[1], naxes[0]);

  fits_read_img(fptr, TDOUBLE, fpixel, naxes[0]*naxes[1], NULL, myimage, NULL, &status);
  printf("%.10e \n", myimage[1024*512+128]);
  
  int i, j;
  for (i=0; i < naxes[0]*naxes[1]; i++)
    if (myimage[i] <= 4e-8)
      myimage[i] = 0;
      
  
  
  fits_create_file(&gptr, "image_out_inc52.9_PA_226.2_clipped.fits", &status); /* create new file */
  fits_copy_header(fptr, gptr, &status);
  
  
  fits_write_img(gptr, TDOUBLE, fpixel, naxes[0]*naxes[1], myimage, &status);
  fits_close_file(gptr, &status);
  
  exit(0);
  
  return( status );
}
