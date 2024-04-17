#include "dbg.h"
#include "datafile.h"
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>

int main(){

    if (init_logging("debug.log")!= 0){
        return 1;
    }

    clock_t start_load, start_process,  end_load, end_process;
    start_load = clock();
    double cpu_time_used;

    int M = 2;
    int N = 2;
    int O = 0;
    int P = 0;

    char file1[256];
    char file2[256];

    printf("Please specify the name of the csv for the first matrix. \n");
    scanf("%255s", file1);
    getchar();
    printf("Please specify the name the csv file for the second matrix. \n");
    scanf("%255s", file2);
    getchar();

    // File loading of Matrix Values.
    count_rows_cols(file1, &M, &N);
    count_rows_cols(file2, &O, &P);

    gsl_matrix *A = gsl_matrix_alloc(M,N);
    gsl_matrix *B = gsl_matrix_alloc(O,P);

    if( M == N && N == O && N == P){
        load_data(file1, A);
        load_data(file2, B);
    } else {
        printf("You must provide 2 square matrixes for this program as of now.\n");
    }
    end_load = clock();
    cpu_time_used = ((double) (end_load - start_process)) / CLOCKS_PER_SEC;
    process_time(cpu_time_used, "End of the processing time. ");

    start_process = clock ();

    gsl_matrix_mul_elements(A,B);
    double max = gsl_matrix_max(A);
    gsl_matrix_transpose(B);

    for (int i = 0; i < M; i++ ){
        printf("( ");
        for (int j = 0; j < N; j++){
            printf(" %g ", gsl_matrix_get(A,i,j));
        }
        printf(" ) \n");
    }
    printf(" maximum value %f \n", max);
    printf(" The transposed matrix B is: \n");

        for (int i = 0; i < M; i++ ){
        printf("( ");
        for (int j = 0; j < N; j++){
            printf(" %g ", gsl_matrix_get(B,i,j));
        }
        printf(" ) \n");
    }

    // Find the determinant.
    gsl_permutation *p = gsl_permutation_alloc(N);
    int s;
    gsl_linalg_LU_decomp(B, p, &s);
    double det = gsl_linalg_LU_det(B, s);

    printf( "The determinant is: %g\n", det);

    // Find the inverse matrix of B
    if(det == 0){
        printf("The matrix is singular and cannot be inverted. \n");
    } else {
        gsl_matrix *inverse = gsl_matrix_alloc(N, N);
        gsl_linalg_LU_invert(B, p, inverse);

        printf("The inverse matrix is: \n");
        for( int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                printf("%g\t", gsl_matrix_get(inverse, i, j));
            }
            printf("\n");
        }
        gsl_matrix_free(inverse);
    }

    // Find the eigenvalues and eigenvectors of B
    if(M == N && N == O && N == P ){
        gsl_vector *eval = gsl_vector_alloc(N);
        gsl_matrix *evec = gsl_matrix_alloc(N, N);

        gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(N);
        gsl_eigen_symmv(B, eval, evec, w);

        gsl_eigen_symmv_free(w);

    // Output of eigenvalues and eigenvectors
        for(int i = 0; i < N; i++){
            double eigenvalue = gsl_vector_get(eval, i);
            gsl_vector_view eigenvector = gsl_matrix_column(evec, i);
            printf("Eigenvalue: %g\n", eigenvalue);
            printf("Eigenvector: \n");
            printf(" ( ");
            for (int j = 0; j < N; j++){
                if(j < N-1){
                printf("%g, ", gsl_vector_get(&eigenvector.vector, j));
                } else {
                printf("%g ", gsl_vector_get(&eigenvector.vector, j));
                }
            }
            printf(")\n");
        }
    } else {
        printf("The matrix is not symmetric and new functionality must be added to handle non-symmetric matrix. \n");
    }

    gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_permutation_free(p);

    end_process = clock();
    cpu_time_used = ((double) (end_process - start_process)) / CLOCKS_PER_SEC;
    process_time(cpu_time_used, "End of the processing time. ");
    cpu_time_used = ((double) (end_process - start_load)) / CLOCKS_PER_SEC;
    process_time(cpu_time_used, "Total runtime. ");


    return 0;

}