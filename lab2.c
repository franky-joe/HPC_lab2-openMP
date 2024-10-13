/* Laboratorio N°2 High Performance Computing
   Autores: Marcelo Alvarez, Franco Cornejo */

#include "fitsio.h" /* Requerido para usar CFITSIO */
#include "stdio.h"
#include "stdlib.h" /* Para malloc y free */
#include "math.h"   /* Para sqrt y log */
#include "time.h"   /* Para inicializar la semilla de rand */
#include "unistd.h" /* Para getopt() */
#include <omp.h>
#include <string.h>

/******************************************************************************
 *                                Utilidades                                  *
 ******************************************************************************/

/* Generar un número aleatorio con distribución normal usando el método de Box-Muller */
double randn(double mu, double sigma) {
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    return z0 * sigma + mu;
}

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
    
    /* Crear un nuevo archivo FITS con la opción de sobrescribir si existe. */
    char overwrite_filepath[256];
    snprintf(overwrite_filepath, sizeof(overwrite_filepath), "!%s", filepath);
    fits_create_file(&gptr, overwrite_filepath, status);

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

/******************************************************************************
 *                              Primer filtro                                 *
 ******************************************************************************/

/* Función para crear el kernel Gaussiano */
void create_gaussian_kernel(double* kernel, int kernel_size, double sigma, double theta) {
    int mid = kernel_size / 2;
    double sigma_sq = sigma * sigma;
    double cos_theta = cos(theta);
    double sin_theta = sin(theta);

    /* Variables para la rotación del kernel */
    double a = cos_theta * cos_theta / (2.0 * sigma_sq) + sin_theta * sin_theta / (2.0 * sigma_sq);
    double b = sin(2.0 * theta) / (4.0 * sigma_sq);
    double c = sin_theta * sin_theta / (2.0 * sigma_sq) + cos_theta * cos_theta / (2.0 * sigma_sq);

    /* Inicializar la suma */
    double sum = 0.0;

    /* Calcular los valores del kernel */
    #pragma omp parallel
    {
        /* Cada hebra tendrá su propia suma parcial */
        double sum_local = 0.0;

        #pragma omp for collapse(2) schedule(static)
        for (int x = -mid; x <= mid; x++) {
            for (int y = -mid; y <= mid; y++) {
                int idx = (x + mid) * kernel_size + (y + mid);
                kernel[idx] = exp(-(a * x * x + 2 * b * x * y + c * y * y));
                sum_local += kernel[idx];  /* Acumular en la suma local */
            }
        }

        /* Reducir la suma parcial de cada hebra en la variable global */
        #pragma omp atomic
        sum += sum_local;
    }

    /* Normalizar el kernel en paralelo */
    #pragma omp parallel for
    for (int i = 0; i < kernel_size * kernel_size; i++) {
        kernel[i] /= sum;
    }
}

/* Función para aplicar la convolución Gaussiana */
void gaussian_convolution(double* image, double* output, long* naxes, double* kernel, int kernel_size) {
    int mid = kernel_size / 2;
    long width = naxes[0];
    long height = naxes[1];

    /* Convolución de la imagen con el kernel Gaussiano */
    #pragma omp parallel for collapse(2)
    for (long x = 0; x < width; x++) {
        for (long y = 0; y < height; y++) {
            double sum = 0.0;
            for (int i = -mid; i <= mid; i++) {
                for (int j = -mid; j <= mid; j++) {
                    long xi = x + i;
                    long yj = y + j;

                    /* Considerar el borde de la imagen */
                    if (xi >= 0 && xi < width && yj >= 0 && yj < height) {
                        sum += image[xi * width + yj] * kernel[(i + mid) * kernel_size + (j + mid)];
                    }
                }
            }
            output[x * width + y] = sum;
        }
    }
}

/******************************************************************************
 *                              Segundo filtro                                *
 ******************************************************************************/

/* Filtro de ruido gaussiano */
void add_gaussian_noise(double* image, long* naxes, double sigma) {
    long num_pixels = naxes[0] * naxes[1];

    /* Añadir ruido gaussiano de forma paralela */
    #pragma omp parallel for
    for (long i = 0; i < num_pixels; i++) {
        double noise = randn(0, sigma); /* Ruido gaussiano con media 0 y desviación sigma */
        image[i] += noise;
    }
}

/******************************************************************************
 *                              Tercer filtro                                 *
 ******************************************************************************/

/* Filtro de normalización */
void normalize_image(double* image, long* naxes, double Jmax) {
    long num_pixels = naxes[0] * naxes[1];
    double min_val = image[0];
    double max_val = image[0];

    /* Encontrar los valores mínimo y máximo de la imagen */
    #pragma omp parallel for reduction(min:min_val) reduction(max:max_val)
    for (long i = 0; i < num_pixels; i++) {
        if (image[i] < min_val) {
            min_val = image[i];
        }
        if (image[i] > max_val) {
            max_val = image[i];
        }
    }

    /* Normalizar la imagen */
    #pragma omp parallel for
    for (long i = 0; i < num_pixels; i++) {
        image[i] = ((image[i] - min_val) / (max_val - min_val)) * Jmax;
    }
}

int main(int argc, char *argv[]) {
    long naxes[2]; /* Para almacenar las dimensiones de la imagen */
    int status = 0; /* Inicializar estado */
    char *input_filename = NULL;  /* Nombre del archivo de entrada */
    char *output_filename = NULL; /* Nombre del archivo de salida */
    double sigma = 0.0; /* Desviación estándar del ruido */
    double fraciontheta = 0.0; /* Fracción theta */
    double Jmax = 0.0; /* Intensidad máxima */
    int num_threads = 1; /* Número de hebras */

    /* Procesar argumentos con getopt() */
    int opt;
    while ((opt = getopt(argc, argv, "i:o:s:t:j:T:")) != -1) {
        switch (opt) {
            case 'i':
                input_filename = optarg;
                break;
            case 'o':
                output_filename = optarg;
                break;
            case 's':
                sigma = atof(optarg);
                break;
            case 't':
                fraciontheta = atof(optarg);
                break;
            case 'j':
                Jmax = atof(optarg);
                break;
            case 'T':
                num_threads = atoi(optarg);
                break;
            default:
                fprintf(stderr, "Uso: %s -i archivo_entrada -o archivo_salida -s sigma -t fraciontheta -j intensidad -T numerohebras\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    if (input_filename == NULL || output_filename == NULL || sigma == 0.0 || Jmax == 0.0 || fraciontheta == 0.0) {
        fprintf(stderr, "Todos los parámetros son obligatorios.\n");
        exit(EXIT_FAILURE);
    }

    /* Configurar el número de hebras con OpenMP */
    omp_set_num_threads(num_threads);

    /* Inicializar la semilla para generar números aleatorios */
    srand(time(NULL));

    /* Leer la imagen desde un archivo FITS */
    double* myimage = read_fits_image(input_filename, naxes, &status);
    if (myimage == NULL) {
        printf("Error al leer la imagen\n");
        return status;
    }

    /* Imprimir las dimensiones de la imagen */
    printf("Tamaño de la imagen: (%ld, %ld)\n", naxes[1], naxes[0]);

    /* Parámetros del kernel Gaussiano */
    int kernel_size = 27;
    double sigma_kernel = 3.5;
    double theta_kernel = fraciontheta * M_PI;

    /* Crear el kernel Gaussiano */
    double* kernel = (double*) malloc(kernel_size * kernel_size * sizeof(double));
    create_gaussian_kernel(kernel, kernel_size, sigma_kernel, theta_kernel);

    /* Crear la imagen de salida para la convolución */
    double* smoothed_image = (double*) malloc(naxes[0] * naxes[1] * sizeof(double));

    /* Aplicar la convolución Gaussiana (paralelizada con OpenMP) */
    gaussian_convolution(myimage, smoothed_image, naxes, kernel, kernel_size);

    /* Liberar la memoria del kernel */
    free(kernel);

    /* Reemplazar la imagen original con la suavizada */
    memcpy(myimage, smoothed_image, naxes[0] * naxes[1] * sizeof(double));
    free(smoothed_image);

    /* Aplicar el filtro de ruido gaussiano */
    add_gaussian_noise(myimage, naxes, sigma);

    /* Aplicar la normalización */
    normalize_image(myimage, naxes, Jmax);

    /* Escribir la imagen modificada en un nuevo archivo FITS */
    write_fits_image(output_filename, naxes, myimage, &status);
    if (status) {
        printf("Error al escribir la imagen\n");
    }

    /* Liberar la memoria de la imagen */
    free(myimage);

    return status;
}
