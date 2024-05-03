#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>
#include <math.h>
#include "parameters.h"
#include "rng.h"
#include "rng_functions.h"
#include "general_functions.h"

/**
 * LOG2FACT is an array where each element at index i contains the floor of the binary logarithm of i!.
 */
const int LOG2FACT[65] = {0, 0, 0, 2, 4, 6, 9, 12, 15, 18, 21, 25, 28, 32, 36, 40, 44, 48, 52, 56, 61, 65, 69, 74, 79, 83, 88, 93, 97, 102, 107, 112, 117, 122, 127, 132, 138, 143, 148, 153, 159, 164, 169, 175, 180, 186, 191, 197, 202, 208, 214, 219, 225, 231, 237, 242, 248, 254, 260, 266, 272, 278, 284, 289, 295};

 /**
 * This function generates the matrix B_22 based on the construction described in the paper. It involves:
 * 1. Generating random elementary matrices and permuting them.
 * 2. Multiplying these matrices iteratively to form B_22.
 * 3. Finally, multiplying it with a hard to guess matrix to finalize B_22.
 *
 * @param B22   A pointer to the 3D matrix that will store the resulting B_22 matrix.
 */
void generate_B22(int64_t*** B22)
{
	int64_t*** E = allocate_ring_matrix(S, S);
	int64_t*** PE = allocate_ring_matrix(S, S);
	int64_t*** T = allocate_ring_matrix(S, S);
	identity_ring_matrix(S, T);
	
	for(int r=0; r<K; r++)
	{
		int i = rng_byte()%SF;
		int k = rng_byte()%M;
		
		identity_ring_matrix(S, E);
		E[P[i][0]][P[i][1]][k] = rng2();
		
		i = rng_byte()%SF;
		row_permute(i, E, PE);
		
		rmm_multiply(S, S, S, T, PE, B22);
		copy_ring_matrix(S, S, B22, T);
	}
	
	int idx[M];
 	range_vector(M, idx);
 	
 	int64_t x[M];
 	int64_t y[M];
 	zero_vector(M, x);
 	zero_vector(M, y);
 	
 	permute_vector(LF, M, idx);
 	for(int k=0; k<LF; k++)
		x[idx[k]] = rng2();
		
	permute_vector(LF, M, idx);
 	for(int k=0; k<LF; k++)
		y[idx[k]] = rng2();
	
	identity_ring_matrix(S, E);
	
	for(int k=0; k<M; k++)
    {
    	E[0][1][k] = -y[k];
    	E[1][0][k] = x[k];
    	E[1][2][k] = y[k];
    	E[2][1][k] = x[k];
	}	

	rmm_multiply(S, S, S, T, E, B22);
	
	free_ring_matrix(S, S, E); E = NULL;
	free_ring_matrix(S, S, PE); PE = NULL;
    free_ring_matrix(S, S, T); T = NULL;
}

/**
 * The conservative security metric descripted in the paper. This metric provides an estimate of the complexity
 * involved in guessing a polynomial under conservative assumptions.
 *
 * @param u   A pointer to the polynomial coefficients, whose guessing security to is be estimated.
 * @return    The computed base 2 conservative security metric as an integer.
 */
int security_metric(int64_t* u)
{
	int bits_security = LOG2FACT[M];
	
	int64_t minv = u[0];
	for(int k=1; k<M; k++)
		if(u[k]<minv)
			minv = u[k];
			
	int64_t maxv = u[0];
	for(int k=1; k<M; k++)
		if(u[k]>maxv)
			maxv = u[k];
	
	for(int i=minv; i<=maxv; i++)
	{
		int count = 0;
		for(int k=0; k<M; k++)
			if(u[k]==i)
				count++;
		
		bits_security -= LOG2FACT[count];
	}
	
	return bits_security;
}

/**
 * This function validates the matrix B_22 based on the guessing complexity of each entry.
 * It checks each polynomial stored in the matrix B_22 to ensure that the guessing complexity of
 * each polynomial meets or exceeds a defined threshold G. If any element's guessing complexity
 * is below the threshold G, the matrix is deemed invalid for the scheme.
 *
 * @param B22   A pointer to the 3D matrix B_22 that is to be validated.
 * @return      A boolean value; true if all polynomials in B_22 meet the security threshold, false otherwise.
 */
bool valid_B22(int64_t*** B22)
{
	for(int i=0; i<S; i++)
		for(int j=0; j<S; j++)
			if(security_metric(B22[i][j])<G)
				return false;
	
	return true;
}

/**
 * This function polpulates B with B_22 and also generates the submatrix B_21 that is part of B.
 *
 * @param B22   A pointer to the 3D matrix B_22 that is to be copied to B.
 * @param B     A pointer to the 3D matrix B that gets populated.
 */
void generate_B(int64_t*** B22, int64_t*** B)
{
 	zero_ring_matrix(N, N, B);
 	B[0][0][0] = 1;
 	
 	int idx[M];
 	range_vector(M, idx);
 	
    for(int i=1; i<N; i++)
    {
    	for(int j=0; j<R; j++)
    	{
    		permute_vector(LB, M, idx);
    		
    		for(int k=0; k<LB; k++)
		    	B[i][j][idx[k]] = rng4();
		}	
	}
    	
    for(int i=0; i<S; i++)
    	for(int j=0; j<S; j++)
    		for(int k=0; k<M; k++)
    			B[R+i][R+j][k] = B22[i][j][k];
}

/**
 * This function uses B and J = {1, 1, -1, -1} to compute C as described in the paper.
 *
 * @param B     A pointer to the 3D matrix B used in the computation.
 * @param C     A pointer to the 3D matrix C that gets computed.
 */
void compute_C(int64_t*** B, int64_t*** C)
{
	int64_t*** BTJ = allocate_ring_matrix(N, N);
	
	for(int i=0; i<N; i++)
		for(int j=0; j<N; j++)
			for(int k=0; k<M; k++)
				BTJ[j][i][k] = B[i][j][k];
	
	for(int i=0; i<N; i++)
		for(int j=S-R; j<N; j++)
			for(int k=0; k<M; k++)
				BTJ[i][j][k] *= -1;
	
	rmm_multiply(N, N, N, BTJ, B, C);
	
	free_ring_matrix(N, N, BTJ); BTJ = NULL;
}

/**
 * This function checks if the matrix C is valid by ensuring that all its elements are within predefined bounds.
 * The matrix C is divided conceptually into different regions, each with a specific bound:
 * 1. The top-left element (at index [0][0]) must have all its polynomial coefficients' absolute values less than C1_BOUND.
 * 2. The first row (excluding the top-left element) must have polynomial coefficients' absolute values less than C2_BOUND.
 * 3. The elements in the upper triangle of the matrix (excluding the first row) must have polynomial coefficients' absolute values less than C3_BOUND.
 *
 * @param C   A pointer to the 3D matrix C that needs to be validated.
 * @return    A boolean value; true if all parts of the matrix C meet their respective bounds, false otherwise.
 */
bool valid_C(int64_t*** C)
{
	// C1_BOUND
	for(int k=0; k<M; k++)
		if(abs(C[0][0][k]) >= C1_BOUND)
	        return false;

	
	// C2_BOUND
	for(int j=1; j<N; j++)
		for(int k=0; k<M; k++)
			if(abs(C[0][j][k]) >= C2_BOUND)
				return false;
									
	// C3_BOUND
	for(int i=1; i<N; i++)
		for(int j=i; j<N; j++)
			for(int k=0; k<M; k++)
				if(abs(C[i][j][k]) >= C3_BOUND)
					return false;
		
	return true;	
}

/**
 * This function encodes the matrix C into a compact byte array pk (public key). The encoding is done by first offsetting
 * the elements by their respective bounds and then representing them as binary bits of specified lengths.
 * The matrix is divided into different sections, each treated with specific bit constraints:
 * 1. The top-left element of C uses C1_BITS per coefficient.
 * 2. The first row (excluding the top-left element) uses C2_BITS per coefficient.
 * 3. The upper triangle of the matrix (excluding the first row) uses C3_BITS per coefficient.
 * The function compacts these bits into bytes, storing the result in the byte array pk to minimize storage requirements.
 *
 * @param C    A pointer to the 3D matrix C containing the coefficients to be encoded.
 * @param pk   A pointer to the public key byte array where the encoded matrix C will be stored.
 */
void C_to_pk(int64_t*** C, unsigned char* pk)
{	
	int total_bits = C1_BITS*M + C2_BITS*M*S + C3_BITS*M*N*S/2;
	bool C_bits[total_bits];
	
	int bitc = 0;
	
	// C1 bits
	for(int k=0; k<M; k++)
	{
		int64_t val = C[0][0][k] + C1_BOUND;
		
		for(int bit=0; bit<C1_BITS; bit++)
        	C_bits[bitc++] = (val >> bit) & 1;
	}
	
	// C2 bits
	for(int j=1; j<N; j++)
	{
		for(int k=0; k<M; k++)
		{
			int64_t val = C[0][j][k] + C2_BOUND;
		
			for(int bit=0; bit<C2_BITS; bit++)
        		C_bits[bitc++] = (val >> bit) & 1;
		}
	}
	
	// C3 bits
	for(int i=1; i<N; i++)
	{
		for(int j=i; j<N; j++)
		{
			for(int k=0; k<M; k++)
			{
				int64_t val = C[i][j][k] + C3_BOUND;
		
				for(int bit=0; bit<C3_BITS; bit++)
        			C_bits[bitc++] = (val >> bit) & 1;
			}
		}
	}
	
	// C_bits to pk
	int pkc = 0;
	for(int i=0; i<total_bits; i+=8)
	{
		unsigned char val = 0;
		
		for(int j=0; j<8; j++)
			val += C_bits[i+j] << j;
				
		pk[pkc++] = val;
	}
}

/**
 * This function generates a pair of public and private keys for DEFI
 *
 * @param pk   A pointer to the public key which stores matrix C efficiently.
 * @param sk   A pointer to the private key which stores the seed to reproduce the secret matrix B.
 * @return 0 if the function was successful
 */
int key_gen(unsigned char *pk, unsigned char *sk)
{
	int64_t*** B22 = allocate_ring_matrix(S, S);
	int64_t*** B = allocate_ring_matrix(N, N);
	int64_t*** C = allocate_ring_matrix(N, N);
	
	do
	{
		do
		{
			randombytes(sk, 48);
    		initialize_rng(sk, 48);
    		generate_B22(B22);
		}
		while(valid_B22(B22)==false);
		
		generate_B(B22, B);
    	compute_C(B, C);
	}
	while(valid_C(C)==false);
	
	clear_rng();
	
	C_to_pk(C, pk);
	
	free_ring_matrix(S, S, B22); B22 = NULL;
	free_ring_matrix(N, N, B); B = NULL;
	free_ring_matrix(N, N, C); C = NULL;
	
    return 0;
}


