/*
Diego-Santiago-Gutierrez
Histogram equalization.c
*/

//Librerias que vamos a utilizar 
#include <stdio.h>
#include<stdlib.h>
#include<omp.h>
#include<math.h>

//Stn master nos ayudará a manejar las imagenes 
#define STB_IMAGE_IMPLEMENTATION
#include "stb-master/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb-master/stb_image_write.h"

//Caracteristicas  de una imagen 
#define COLOR_VALUE 256
#define QUALITY_IMAGE 100

//FUNCIONES DE CONTROL SECUENCIAL-PARALELO
double equalize_image_sequential(unsigned char *image, int width, int height, int channels, char *name_UNsuffix);
double equalize_image_parallel(unsigned char *image, int width, int height, int number_of_channels, char *name_UNsuffix);

//-----------FUNCIONES COMPARTIDAS------------
//############################################
//ARREGLOS VACIOS
long *empty_array_LONG(long size);
unsigned char *empty_array_UC(long size);
//RECONOCIMIENTO DE ARCHIVOS
char *filename_UNsuffix(char *path);
char *new_file(char *original, char *suffix);
//VALOR MINIMO DEL CDF
long min_cdf(long *original_histogram);

//----------------SECUENCIAL---------------// 
//###########################################
//GENERAR ARCHIVO CSV
void csv_secuencial(long *original_histogram, long *equalized_histogram);
//PROCESOS SECUENCIALES PARA PROCESARLA IMAGEN 
void cdf_SECUENCIAL(long *cdf, long *original_histogram);
void equalized_cdf_SECUENCIAL(long *cdf, long minimum, long size);
void image_SEQUENTIAL(unsigned char *original_image, unsigned char *new_image, int channels, long *equalized_cdf, long size);
void histogram_SEQUENTIAL(long *original_histogram, unsigned char *image, long size, int channels);

//----------------PARALELO---------------// 
//#########################################
//GENERA ARCHIVO CSV
void csv_parallel(long *original_histogram, long *equalized_histogram);
//PROCESOS PARALELOS PARA PROCESARLA IMAGEN 
void equalized_cdf_PARALLEL(long *cdf, long minimum, long size);
void image_PARALLEL(unsigned char *original_image, unsigned char *new_image, int number_of_channels, long *equalized_cdf, long size);
void histogram_PARALLEL(long *original_histogram, unsigned char *image, long size, int number_of_channels);

///////////////////////////////////////////
int main(int argc, char *argv[])
{
    if (argc < 2){
        printf("ERROR AL PASAR ARGUMENTO");
        return -1;
    }

    char *image_source = argv[1];
    int width, height, channels;

    double t1 = omp_get_wtime();
    unsigned char *image = stbi_load(image_source, &width, &height, &channels, 0);
    double t2 = omp_get_wtime();


    printf("-------IMAGEN LEIDA----------");
    char *name_UNsuffix = filename_UNsuffix(image_source);
    printf("Nombre de la imagen: %s\n", name_UNsuffix);
    printf("Tiempo de carga: %f\n", t2 - t1);
    printf("ANCHO: %d\t ALTO: %d\n", width, height);
    printf("Numero de canales: %d\n", channels);
    printf("Tamano: %d\n", (width * height)/1000);

    double sequential_time = equalize_image_sequential(image, width, height, channels, name_UNsuffix);
    double parallel_time = equalize_image_parallel(image, width, height, channels, name_UNsuffix);

    int total_processors = omp_get_num_procs();
    double speedup = sequential_time / parallel_time;
    double efficiency = speedup / total_processors;
    double overhead = parallel_time - sequential_time / total_processors;

    printf("\nNumero de procesadores: %d\n", total_processors);
    printf("SECUENCIAL: %f\t PARALELO: %f\n", sequential_time, parallel_time);
    printf("Speedup: %f\n", speedup);
    printf("Eficiencia: %f\n", efficiency);
    printf("Overhead: %f\n", overhead);

    stbi_image_free(image);
    free(name_UNsuffix);
    return 0;
}


unsigned char *empty_array_UC(long size){
    unsigned char *empty_array_UC = malloc(sizeof(unsigned char) * size);
    for (int i = 0; i < size; ++i){
        empty_array_UC[i] = 0;
    }
    return empty_array_UC;
}
long *empty_array_LONG(long size){
    long *empty_array_LONG = malloc(sizeof(long) * size);
    for (int i = 0; i < size; ++i){
        empty_array_LONG[i] = 0;
    }
    return empty_array_LONG;
}

long min_cdf(long *histogram){
    for (int i = 0; i < COLOR_VALUE; ++i){
        if (histogram[i] != 0) return histogram[i];
    }
    return -1;
}

char *filename_UNsuffix(char *path){
    char *filename_UNsuffix = strrchr(path, '/');
    if (filename_UNsuffix == NULL) filename_UNsuffix = path;
    else filename_UNsuffix++;
    return strtok(filename_UNsuffix, ".");
}

char *new_file(char *original, char *suffix)
{
    char *new_file = malloc(sizeof(char) * 1000);
    snprintf(new_file, sizeof(char) * 1000, "%s%s", original, suffix);
    return new_file;
}

//---------------SECUENCIAL-----------------

void csv_secuencial(long *original_histogram, long *equalized_histogram){   
    FILE *csv = fopen("histograma_secuencial.csv", "w+");
    for (int i = 0; i < COLOR_VALUE; i++){
        fprintf(csv, "%d,%ld,%ld\n", i, original_histogram[i], equalized_histogram[i]);
    }
    fclose(csv);
}

void cdf_SECUENCIAL(long *cdf, long *original_histogram){
    cdf[0] = original_histogram[0];
    for (int i = 1; i < COLOR_VALUE; ++i)
    {
        cdf[i] = cdf[i - 1] + original_histogram[i];
    }
}

void equalized_cdf_SECUENCIAL(long *cdf, long minimum, long size)
{
    double numerator = (double) COLOR_VALUE - 2;
    double denominator = (double) size - minimum;

    for (int i = 0; i < COLOR_VALUE; ++i)
    {
        cdf[i] = (long) round( (double)(cdf[i] - minimum) * (numerator / denominator))+ 1;
    }
}

double equalize_image_sequential(unsigned char *image, int width, int height, int channels, char *name_UNsuffix){
    long size = width * height;
    long *original_histogram = empty_array_LONG(COLOR_VALUE);
    long *new_histogram = empty_array_LONG(COLOR_VALUE);
    long *cdf = empty_array_LONG(COLOR_VALUE);

    unsigned char *new_image = empty_array_UC(size * channels);

    double t1 = omp_get_wtime();
    histogram_SEQUENTIAL(original_histogram, image, size, channels);
    cdf_SECUENCIAL(cdf, original_histogram);
    long min_cdf_value = min_cdf(cdf);
    equalized_cdf_SECUENCIAL(cdf, min_cdf_value, size);
    image_SEQUENTIAL(image, new_image, channels, cdf, size);
    histogram_SEQUENTIAL(new_histogram, new_image, size, channels);
    double t2 = omp_get_wtime();

    double t1_csv = omp_get_wtime();
    csv_secuencial(original_histogram, new_histogram);
    double t2_csv = omp_get_wtime();

    char *new_sequential_name = new_file(name_UNsuffix, "_eq_sequential.jpg");

    double t1_image = omp_get_wtime();
    stbi_write_jpg(new_sequential_name, width, height, channels, new_image, QUALITY_IMAGE);
    double t2_image = omp_get_wtime();

    printf("\nGENERACION CSV (secuencial): %f\n", t2_csv - t1_csv);
    printf("GENERACION IMAGEN: (secuencial): %f\n", t2_image - t1_image);

    free(original_histogram);
    free(new_histogram);
    free(cdf);
    free(new_image);

    return t2 - t1;
}

void image_SEQUENTIAL(unsigned char *original_image, unsigned char *new_image, int channels, long *equalized_cdf, long size){
    for (int i = 0; i < size * channels; i += channels){
        unsigned char original_value = original_image[i];
        for (int j = i; j < channels + i; ++j){
            new_image[j] = (unsigned char) equalized_cdf[original_value];
        }
    }
}

void histogram_SEQUENTIAL(long *original_histogram, unsigned char *image, long size, int channels){
    for (int i = 0; i < size * channels; i += channels){
        original_histogram[image[i]]++;
    }
}

//-----------------------PARALELO-------------------------
void csv_parallel(long *original_histogram, long *equalized_histogram){   
    FILE *csv_parallel = fopen("histograma_parallel.csv", "w+");
    for (int i = 0; i < COLOR_VALUE; i++){
        fprintf(csv_parallel, "%d,%ld,%ld\n", i, original_histogram[i], equalized_histogram[i]);
    }
    fclose(csv_parallel);
}

void equalized_cdf_PARALLEL(long *cdf, long minimum, long size){
    double numerator = (double) COLOR_VALUE - 2;
    double denominator = (double) size - minimum;
    #pragma omp parallel for
    for (int i = 0; i < COLOR_VALUE; ++i)
    {
        cdf[i] = (long) round( (double)(cdf[i] - minimum) * (numerator / denominator)) + 1;
    }
}

double equalize_image_parallel(unsigned char *image, int width, int height, int number_of_channels, char *name_UNsuffix)
{
    long size = width * height;
    long *original_histogram = empty_array_LONG(COLOR_VALUE);
    long *new_histogram = empty_array_LONG(COLOR_VALUE);
    long *cdf = empty_array_LONG(COLOR_VALUE);
    unsigned char *new_image = empty_array_UC(size * number_of_channels);

    int total_processors = omp_get_num_procs();
    omp_set_num_threads(total_processors);

    double t1 = omp_get_wtime();
    histogram_PARALLEL(original_histogram, image, size, number_of_channels);
    cdf_SECUENCIAL(cdf, original_histogram);
    long min_cdf_value = min_cdf(cdf);
    equalized_cdf_PARALLEL(cdf, min_cdf_value, size);
    image_PARALLEL(image, new_image, number_of_channels, cdf, size);
    histogram_PARALLEL(new_histogram, new_image, size, number_of_channels);
    double t2 = omp_get_wtime();

    double t1_csv = omp_get_wtime();
    csv_parallel(original_histogram, new_histogram);
    double t2_csv = omp_get_wtime();

    char *new_parallel_name = new_file(name_UNsuffix, "_eq_parallel.jpg");

    double t1_image = omp_get_wtime();
    stbi_write_jpg(new_parallel_name, width, height, number_of_channels, new_image, QUALITY_IMAGE);
    double t2_image = omp_get_wtime();

    printf("\nPARALELO CSV: %f\n", t2_csv - t1_csv);
    printf("IMAGEN PARALALELA: %f\n", t2_image - t1_image);

    free(original_histogram);
    free(new_histogram);
    free(cdf);
    free(new_image);

    return t2 - t1;
}

void image_PARALLEL(unsigned char *original_image, unsigned char *new_image, int number_of_channels, long *equalized_cdf, long size){
    #pragma omp parallel for
    for (int i = 0; i < size * number_of_channels; i += number_of_channels){
        unsigned char original_value = original_image[i];
        for (int j = i; j < number_of_channels + i; ++j)
        {
            new_image[j] = (unsigned char) equalized_cdf[original_value];
        }
    }
}

void histogram_PARALLEL(long *original_histogram, unsigned char *image, long size, int number_of_channels){
    #pragma omp parallel for reduction(+:original_histogram[:COLOR_VALUE])
    for (int i = 0; i < size * number_of_channels; i += number_of_channels){
        original_histogram[image[i]]++;
    }
}
