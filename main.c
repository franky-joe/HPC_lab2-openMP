#include "fitsio.h" /* Requerido para usar CFITSIO */
#include "stdio.h"
#include "stdlib.h" /* Para malloc y free */


/******************************************************************************
 *                                Utilidades                                  *
 ******************************************************************************/


/* Función para leer la imagen desde un archivo FITS y devolver un puntero a la lista de píxeles */
double* read_fits_image(const char* filepath, long* naxes, int* status) {
    fitsfile *fptr;
    long fpixel = 1;
    
    /* Abrir el archivo FITS en modo de solo lectura */
    fits_open_file(&fptr, filepath, READONLY, status);
    if (*status) {
        printf("Error al abrir el archivo: %s\n", filepath);
        return NULL;
    }
    
    /* Obtener las dimensiones de la imagen */
    fits_get_img_size(fptr, 2, naxes, status);
    
    /* Reservar memoria para la imagen */
    double *myimage = (double *) malloc(naxes[0] * naxes[1] * sizeof(double));
    if (myimage == NULL) {
        printf("Error al asignar memoria\n");
        return NULL;
    }
    
    /* Leer la imagen completa desde el archivo FITS */
    fits_read_img(fptr, TDOUBLE, fpixel, naxes[0] * naxes[1], NULL, myimage, NULL, status);
    if (*status) {
        printf("Error al leer la imagen\n");
        free(myimage);
        return NULL;
    }

    /* Cerrar el archivo FITS */
    fits_close_file(fptr, status);

    return myimage;
}


/* Función para escribir una imagen en un archivo FITS */
void write_fits_image(const char* filepath, long* naxes, double* myimage, int* status) {
    fitsfile *gptr;
    long fpixel = 1;
    
    /* Crear un nuevo archivo FITS con la opción de sobrescribir si existe */
    fits_create_file(&gptr, filepath, status);
    if (*status) {
        printf("Error al crear el archivo: %s\n", filepath);
        return;
    }

    /* Crear la imagen primaria en el nuevo archivo (double) */
    fits_create_img(gptr, DOUBLE_IMG, 2, naxes, status);
    
    /* Escribir los datos de la imagen en el archivo FITS */
    fits_write_img(gptr, TDOUBLE, fpixel, naxes[0] * naxes[1], myimage, status);
    if (*status) {
        printf("Error al escribir la imagen\n");
        return;
    }

    /* Cerrar el archivo FITS */
    fits_close_file(gptr, status);
}


int main() {
    long naxes[2]; /* Para almacenar las dimensiones de la imagen */
    int status = 0; /* Inicializar estado */
    
    /* Leer la imagen desde un archivo FITS */
    double* myimage = read_fits_image("imagein1.fits", naxes, &status);
    if (myimage == NULL) {
        printf("Error al leer la imagen\n");
        return status;
    }

    /* Imprimir las dimensiones de la imagen */
    printf("Tamaño de la imagen: (%ld, %ld)\n", naxes[1], naxes[0]);

    /* Imprimir los píxeles de la imagen */
    for (int ii = 0; ii < naxes[1]; ii++) {
        for (int jj = 0; jj < naxes[0]; jj++) {
            printf("%.10e ", myimage[ii * naxes[0] + jj]);
        }
        printf("\n"); /* Imprimir una nueva línea al final de cada fila */
    }

    /* Realizar alguna modificación en la imagen (opcional) */
    for (int i = 0; i < naxes[0] * naxes[1]; i++) {
        if (myimage[i] <= 4e-8)
            myimage[i] = 0;
    }

    /* Escribir la imagen modificada en un nuevo archivo FITS */
    write_fits_image("salida2.fits", naxes, myimage, &status);
    if (status) {
        printf("Error al escribir la imagen\n");
    }

    /* Liberar la memoria de la imagen */
    free(myimage);

    return status;
}
