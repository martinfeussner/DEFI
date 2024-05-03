#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include "parameters.h"
#include "rng_functions.h"
#include "keccak.h"

/**
 * This tables stores all permutations for s=3. It is used when trying to permute rows/columns of matrices.
 */
const int P[6][3] = {{0,1,2}, {0,2,1}, {1,0,2}, {1,2,0}, {2,0,1}, {2,1,0}};

/**
 * This function allocates memory for a 2D vector where each entry is to store a polynomial of degree < M.
 *
 * @param n   The length of the vector.
 * @return A pointer to the first element of the newly allocated 2D vector.
 *         If the allocation fails at any point, any memory that was 
 *         successfully allocated is freed and NULL is returned.
 */
int64_t** allocate_ring_vector(int n)
{
    int64_t** A = malloc(n*sizeof(int64_t*));
    
    if(A==NULL)
    {
    	printf("Memory Allocation Failure!");
        exit(2);
	}
    
    for(int i=0; i<n; i++)
    {
        A[i] = calloc(M, sizeof(int64_t));
        
        if(A[i]==NULL)
        {
            for(int k=0; k<i; k++)
            	free(A[k]);
            
            free(A);
            printf("Memory Allocation Failure!");
            exit(2);
        }
    }
    
    return A;
}

/**
 * This function allocates memory for a 3D matrix where each entry is to store a polynomial of degree < M.
 *
 * @param m   The number of rows in the matrix.
 * @param n   The number of columns in the matrix.
 * @return A pointer to the first element of the newly allocated 3D matrix.
 *         If the allocation fails at any point, any memory that was 
 *         successfully allocated is freed and NULL is returned.
 */
int64_t*** allocate_ring_matrix(int m, int n)
{
    int64_t*** A = malloc(m*sizeof(int64_t**));
    
    if(A==NULL)
    {
    	printf("Memory Allocation Failure!");
        exit(2);
	}
    
    for(int i=0; i<m; i++)
	{
        A[i] = malloc(n*sizeof(int64_t*));
        
        if(A[i]==NULL)
		{
            for(int k=0; k<i; k++)
            	free(A[k]);
                
            free(A);
            printf("Memory Allocation Failure!");
            exit(2);
        }
        
        for(int j=0; j<n; j++)
		{
            A[i][j] = calloc(M, sizeof(int64_t));
            
            if(A[i][j]==NULL)
			{
                for(int k=0; k<j; k++)
                	free(A[i][k]);
                    
                free(A[i]);
                
                for(int k=0; k<i; k++)
                {
                	for(int l=0; l<n; l++)
                    	free(A[k][l]);
                    
                	free(A[k]);
				}
                    
               	free(A);
	        	printf("Memory Allocation Failure!");
	            exit(2);
            }
        }
    }
    
    return A;
}

/**
 * This function frees the memory for a 2D vector.
 *
 * @param m   The length of the vector.
 * @param A   A pointer to the first element of the matrix to be freed.
 */
void free_ring_vector(int n, int64_t** A)
{
    for(int i=0; i<n; i++)
    	free(A[i]);
    
    free(A);
}

/**
 * This function frees memory for a 3D matrix.
 *
 * @param m   The number of rows in the matrix.
 * @param n   The number of columns in the matrix.
 * @param A   A pointer to the first element of the 3D matrix to be deallocated.
 */
void free_ring_matrix(int m, int n, int64_t*** A)
{
    for(int i=0; i<m; i++)
    {
        for(int j=0; j<n; j++)
        	free(A[i][j]);
        
        free(A[i]);
    }
    
    free(A);
}

/**
 * This function sets all elements of a vector to zero.
 *
 * @param n   The length of the vector.
 * @param A   A pointer to the first element of the vector.
 */
void zero_vector(int n, int64_t* A)
{
	for(int i=0; i<n; i++)
		A[i] = 0;
}

/**
 * This function sets all elements of a 2D vector to zero.
 *
 * @param n   The length of the vector.
 * @param A   A pointer to the first element of the 2D vector to be zeroed.
 */
void zero_ring_vector(int n, int64_t** A)
{
	for(int i=0; i<n; i++)
		for(int j=0; j<M; j++)
			A[i][j] = 0;
}

/**
 * This function sets all elements of a 3D matrix to zero.
 *
 * @param m   The number of rows in the matrix.
 * @param n   The number of columns in the matrix.
 * @param A   A pointer to the first element of the 3D matrix to be zeroed.
 */
void zero_ring_matrix(int m, int n, int64_t*** A)
{
	for(int i=0; i<m; i++)
		for(int j=0; j<n; j++)
			for(int k=0; k<M; k++)
				A[i][j][k] = 0;
}

/**
 * This function sets a 3D square matrix to a 3D identity matrix.
 *
 * @param n   The number of rows/coloumns of the matrix.
 * @param A   A pointer to the first element of the 3D square matrix.
 */
void identity_ring_matrix(int n, int64_t*** A)
{
	zero_ring_matrix(n, n, A);		
	for(int i=0; i<n; i++)
		A[i][i][0] = 1;
}

/**
 * This function copies a 3D matrix to another 3D matrix.
 *
 * @param m   The number of rows in the matrix.
 * @param n   The number of columns in the matrix.
 * @param A   A pointer to the first element of the 3D matrix to be copied.
 * @param B   A pointer to the first element of the 3D matrix storing the copy.
 */
void copy_ring_matrix(int m, int n, int64_t*** A, int64_t*** B)
{
	for(int i=0; i<m; i++)
		for(int j=0; j<n; j++)
			for(int k=0; k<M; k++)
				B[i][j][k] = A[i][j][k];
}

/**
 * This function initializes an vector with consecutive integers starting from 0 up to n-1.
 *
 * @param n   The number of elements in the vector.
 * @param A   A pointer to the first element of the vector to be initialized.
 */
void range_vector(int n, int* A)
{
	for(int i=0; i<n; i++)
		A[i] = i;		
}

/**
 * This function permutes the first m elements of a vector such that after
 * permutation, these elements contain m random selections from the total n
 * elements of the vector. The permutation is done in-place, using a swapping
 * mechanism that ensures each of the first m elements is randomly exchanged
 * with another element from the remaining part of the vector, or itself, ensuring
 * a random sample of m elements from the range of n without replacement.
 *
 * @param m   The number of elements in the vector to be randomly selected and permuted.
 * @param n   The total number of elements in the vector, which defines the range limit for the random index.
 * @param A   A pointer to the first element of the vector that will be permuted.
 */
void permute_vector(int m, int n, int* A)
{
	for(int i=0; i<m; i++)
	{
		int j = i + rng(n-i);
		int t = A[j];
		A[j] = A[i];
		A[i] = t;
	}
}

/**
 * This function permutes the rows of a 3D matrix based on a predefined permutation lookup table.
 * It reorders the rows of the matrix E according to the permutation defined by P[idx], and stores
 * the result in matrix EP. The permutation affects all layers of the matrix uniformly.
 *
 * @param idx   The index into the permutation table P, specifying which permutation to apply.
 * @param E     A pointer to the first element of the 3D matrix whose rows are to be permuted.
 * @param EP    A pointer to the first element of the 3D matrix that will store the permuted result.
 */
void row_permute(int idx, int64_t*** E, int64_t*** EP)
{
	for(int i=0; i<S; i++)
		for(int j=0; j<S; j++)
			for(int k=0; k<M; k++)
				EP[i][j][k] = E[P[idx][i]][j][k];
}

/**
 * This function permutes the columns of a 3D matrix based on a predefined permutation lookup table.
 * It reorders the columns of the matrix E according to the permutation defined by P[idx], and stores
 * the result in matrix EP. The permutation affects all layers of the matrix uniformly.
 *
 * @param idx   The index into the permutation table P, specifying which permutation to apply.
 * @param E     A pointer to the first element of the 3D matrix whose columns are to be permuted.
 * @param EP    A pointer to the first element of the 3D matrix that will store the permuted result.
 */
void col_permute(int idx, int64_t*** E, int64_t*** EP)
{
	for(int i=0; i<S; i++)
		for(int j=0; j<S; j++)
			for(int k=0; k<M; k++)
				EP[i][j][k] = E[i][P[idx][j]][k];
}

/**
 * This function performs polynomial multiplication over the ring defined by x^M + 1. It multiplies
 * two polynomials, poly1 and poly2, each of degree less than M, and stores the result in result_poly.
 * If the sum of the degrees of two terms exceeds or equals M, the result is adjusted according to
 * the specific ring rules. The function also includes an option to zero out the result_poly
 * before performing the multiplication if the overwrite flag is set to true.
 *
 * @param poly1        A pointer to the first coefficient of the first polynomial.
 * @param poly2        A pointer to the first coefficient of the second polynomial.
 * @param result_poly  A pointer to the first coefficient of the resulting polynomial.
 * @param overwrite    A boolean flag indicating whether to zero out result_poly before computation.
 */
void product_in_ring(int64_t* poly1, int64_t* poly2, int64_t* result_poly, bool overwrite)
{
	if(overwrite==true)
		zero_vector(M, result_poly);
	
	for(int i=0; i<M; i++)
	{
		if(poly1[i]!=0)
		{
			for(int j=0; j<M; j++)
			{
				if(poly2[j]!=0)
				{
					int degree = (i+j)%M;
					
					if(i+j>=M)
						result_poly[degree] -= poly1[i]*poly2[j];
					else
						result_poly[degree] += poly1[i]*poly2[j];
				}	
			}	
		}		
	}
}

/**
 * This function performs matrix-vector multiplication in the context of a polynomial ring defined by x^M + 1.
 * The matrix A and vector b are multiplied, and the result is stored in vector c. Before starting the multiplication,
 * the vector c is zeroed out to ensure clean storage for the result. The multiplication is carried out by taking each
 * row of matrix A and multiplying it with vector b, storing the result in the corresponding entry of vector c using
 * the polynomial ring multiplication rules defined in product_in_ring.
 *
 * @param m    The number of rows in matrix A and the size of vector c.
 * @param l    The number of columns in matrix A and the size of vector b.
 * @param A    A pointer to the first element of the 3D matrix A, which holds polynomial coefficients.
 * @param b    A pointer to the first element of the vector b, which holds polynomial coefficients.
 * @param c    A pointer to the first element of the result vector c, which holds polynomial coefficients.
 */
void rmv_multiply(int m, int l, int64_t*** A, int64_t** b, int64_t** c)
{
	zero_ring_vector(m, c);
	
	for(int i=0; i<m; i++)
			for(int j=0; j<l; j++)
				product_in_ring(A[i][j], b[j], c[i], false);
}

/**
 * This function performs matrix-matrix multiplication in the context of a polynomial ring defined by x^M + 1.
 * The matrices A and B are multiplied, and the result is stored in matrix C. Before starting the multiplication,
 * the matrix C is zeroed out to ensure clean storage for the result. The multiplication is carried out using the
 * standard matrix multiplication algorithm, but each element multiplication is done according to the polynomial
 * ring rules defined in product_in_ring, allowing for polynomial coefficients to be combined appropriately.
 *
 * @param m    The number of rows in matrix A and the result matrix C.
 * @param l    The number of columns in matrix A and the number of rows in matrix B.
 * @param n    The number of columns in matrix B and the result matrix C.
 * @param A    A pointer to the first element of the 3D matrix A, which holds polynomial coefficients.
 * @param B    A pointer to the first element of the 3D matrix B, which holds polynomial coefficients.
 * @param C    A pointer to the first element of the 3D matrix C, which holds the result polynomial coefficients.
 */
void rmm_multiply(int m, int l, int n, int64_t*** A, int64_t*** B, int64_t*** C)
{
	zero_ring_matrix(m, n, C);
	
	for(int i=0; i<m; i++)
		for(int j=0; j<n; j++)
			for(int k=0; k<l; k++)
				product_in_ring(A[i][k], B[k][j], C[i][j], false);
}

 /**
 * This function hashes a message using the SHAKE256 hash function and formats the output into a polynomial representation.
 * The output of the hash is such to fit a specific number of coefficients in polynomial format, where each coefficient
 * is centered about zero. The hash bytes are first calculated to be twice the security level in bits (due to the birthday
 * paradox considerations), then each byte is parsed to extract multiple smaller coefficients based on the bits per entry
 * calculated from the security level and the polynomial degree M. Each extracted coefficient is then adjusted by
 * subtracting half the maximum value of the coefficient size to center the values around zero.
 *
 * @param m     A pointer to the message to be hashed.
 * @param mlen  The length of the message.
 * @param h     A pointer to the output vector storing the appropriate hash, where each entry is a polynomial in the ring.
 */
void hash_of_message(const unsigned char* m, unsigned long long mlen, int64_t* h)
{
	int required_length = (SECURITY*2)/8;
	int bits_per_entry = (SECURITY*2)/M;
    unsigned char required_hash[required_length];
	
    FIPS202_SHAKE256(m, mlen, required_hash, required_length);
	zero_vector(M, h);
	
	int index = 0;
	
	for(int i=0; i<required_length; i++)
    {
    	unsigned char val = required_hash[i];
    	unsigned char mask = (1 << bits_per_entry) - 1;
    	
    	for(int j=0; j<8; j+=bits_per_entry)
    		h[index++] = (val >> j) & mask;
    }
    
    int sub_val = 1<<(bits_per_entry-1);

	for(int k=0; k<M; k++)
		h[k] -= sub_val;
}

