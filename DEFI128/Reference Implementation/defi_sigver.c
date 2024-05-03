#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>
#include <math.h>
#include "parameters.h"
#include "general_functions.h"

/**
 * This function decodes a byte array pk into a 3D matrix C using predefined bit constraints for each section of the matrix.
 * The byte array is unpacked into individual bits, which are then assembled into integer coefficients using bit masking
 * and shifting operations. These coefficients correspond to the entries of the matrix C and are adjusted by subtracting
 * respective bounds (C1_BOUND, C2_BOUND, C3_BOUND) based on their positions within the matrix:
 * 1. The top-left element uses C1_BITS per coefficient.
 * 2. The first row (excluding the top-left element) uses C2_BITS per coefficient.
 * 3. The upper triangle (excluding the first row) uses C3_BITS per coefficient.
 * These sections are decoded sequentially, and symmetry is utilized to fill the lower triangle of the matrix.
 *
 * @param pk  A pointer to the public key byte array containing the encoded data.
 * @param C   A pointer to the 3D matrix C where the decoded data will be stored.
 */
void pk_to_C(const unsigned char* pk, int64_t*** C)
{
	int total_bits = C1_BITS*M + C2_BITS*M*S + C3_BITS*M*N*S/2;
    bool C_bits[total_bits];
    
    // pk to C_bits
    int pkc = 0;
    for(int i=0; i<total_bits; i+=8)
    {
    	for(int j=0; j<8; j++)
    	   	C_bits[i+j] = (pk[pkc] >> j) & 1;
    	   	
    	pkc++;
	}
	
	int bitc = 0;
    	
    // C1
	for(int k=0; k<M; k++)
	{
		int64_t val = 0;
		
		for(int bit=0; bit<C1_BITS; bit++)
        	val |= (int64_t)C_bits[bitc++] << bit;
        
        C[0][0][k] = val - C1_BOUND;
	}
	
	// C2
	for(int j=1; j<N; j++)
	{
		for(int k=0; k<M; k++)
		{
			int64_t val = 0;
		
			for(int bit=0; bit<C2_BITS; bit++)
	        	val |= (int64_t)C_bits[bitc++] << bit;
	        
	        C[0][j][k] = val - C2_BOUND;
		}
	}
	
	//print_ring_matrix(N, N, C);
	
	// C3
	for(int i=1; i<N; i++)
	{
		for(int j=i; j<N; j++)
		{
			for(int k=0; k<M; k++)
			{
				int64_t val = 0;
		
				for(int bit=0; bit<C3_BITS; bit++)
			    	val |= (int64_t)C_bits[bitc++] << bit;
			    
			    C[i][j][k] = val - C3_BOUND;
			}
		}
	}
	
	//print_ring_matrix(N, N, C);
	
	for(int i=0; i<N; i++)
        for(int j=0; j<i; j++)
			for(int k=0; k<M; k++)
				C[i][j][k] = C[j][i][k];
}

/**
 * This function decodes a combined array sm, containing a signature encoded as a series of bits packed into bytes
 * followed by an original message, into the signature matrix y and extracts the original message m. The function
 * first unpacks the bits from the beginning of sm to fill the matrix y, adjusting each coefficient by subtracting
 * the bound Y_BOUND. It then calculates the remaining length of the array, which corresponds to the original message,
 * and copies this part into the output array m.
 *
 * @param sm      The signed message.
 * @param smlen   The length of the signed message.
 * @param m       The output byte array representing the message.
 * @param mlen    A pointer to the length of message.
 * @param y       A pointer to the 2D vector y where the decoded signature coefficients will be stored.
 */
void sm_to_my(const unsigned char* sm, unsigned long long smlen, unsigned char* m, unsigned long long* mlen, int64_t** y)
{   
	int total_bits = Y_BITS*M*S;
	bool bits[total_bits];
    
    // sm to bits
    int smi = 0;
    for(int i=0; i<total_bits; i+=8)
    {
    	for(int j=0; j<8; j++)
    	   	bits[i+j] = (sm[smi] >> j) & 1;
    	   	
    	smi++;
	}
	
	
	int biti = 0;
    	
	for(int i=0; i<S; i++)
	{
		for(int k=0; k<M; k++)
		{
			int64_t val = 0;
		
			for(int bit=0; bit<Y_BITS; bit++)
	        	val |= (int64_t)bits[biti++] << bit;
	        
	        y[i][k] = val - Y_BOUND;
		}
	}
    
    *mlen = smlen - smi;
    
    for(int i=0; i<*mlen; i++)
		m[i] = sm[smi+i];
}

/**
 * Performs signature verification.
 *
 * @param m       The input byte array representing the message.
 * @param mlen    A pointer to the length of the message.
 * @param sm      The input byte array representing the signed message.
 * @param smlen   The length of the signed message.
 * @param pk      The input byte array representing the public key.
 * @return 0 if the verification is successful, -1 if unsuccessful.
 */
int sig_ver(unsigned char *m, unsigned long long *mlen, const unsigned char *sm, unsigned long long smlen, const unsigned char *pk)
{
	int64_t*** C = allocate_ring_matrix(N, N);
	pk_to_C(pk, C);
	
	int64_t** y = allocate_ring_vector(S);
	sm_to_my(sm, smlen, m, mlen, y);
	
	for(int i=0; i<S; i++)
	{
		for(int k=0; k<M; k++)
		{
			if(abs(y[i][k]) >= Y_BOUND)
			{
				free_ring_matrix(N, N, C); C = NULL;
				free_ring_vector(S, y); y = NULL;
				return -1; // Verification Unsuccessfull
			}
		}
	}		
	
	int64_t h[M];
	hash_of_message(m, *mlen, h);
	
	int64_t** z = allocate_ring_vector(N);	
	
	for(int k=0; k<M; k++)
		z[0][k] = h[k];
			
	for(int i=0; i<S; i++)
		for(int k=0; k<M; k++)
			z[R+i][k] = y[i][k];
	
	int64_t** Cz = allocate_ring_vector(N);
	rmv_multiply(N, N, C, z, Cz);

	int64_t result[M];
	zero_vector(M, result);
	
	for(int i=0; i<N; i++)
		product_in_ring(z[i], Cz[i], result, false);
			
	free_ring_matrix(N, N, C); C = NULL;
	free_ring_vector(S, y); y = NULL;
	free_ring_vector(N, z); z = NULL;
	free_ring_vector(N, Cz); Cz = NULL;
	
	for(int k=0; k<M; k++)
		if(result[k]!=0)
			return -1; // Verification Unsuccessfull
	
	return 0; // Verification Successfull
}


