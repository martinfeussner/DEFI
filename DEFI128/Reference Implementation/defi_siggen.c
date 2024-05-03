#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>
#include <math.h>
#include "parameters.h"
#include "keccak.h"
#include "rng_functions.h"
#include "general_functions.h"

 /**
 * This function generates the matrix B_{22}^-1 based on the construction described in the paper. It involves:
 * 1. Generating random elementary matrices and permuting them.
 * 2. Multiplying these matrices iteratively to form B_{22}^-1.
 * 3. Finally, multiplying it with a hard to guess matrix to finalize B_{22}^-1.
 *
 * @param B22inv   A pointer to the 3D matrix that will store the resulting B_{22}^-1 matrix.
 */
void generate_B22inv(int64_t*** B22inv)
{
	int64_t*** E = allocate_ring_matrix(S, S);
	int64_t*** EP = allocate_ring_matrix(S, S);
	int64_t*** T = allocate_ring_matrix(S, S);
	identity_ring_matrix(S, T);
	
	for(int r=0; r<K; r++)
	{
		int i = rng_byte()%SF;
		int k = rng_byte()%M;
		
		identity_ring_matrix(S, E);
		E[P[i][0]][P[i][1]][k] = -rng2();
		
		i = rng_byte()%SF;
		col_permute(i, E, EP);
		
		rmm_multiply(S, S, S, EP, T, B22inv);
		copy_ring_matrix(S, S, B22inv, T);	
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
	
	int64_t x2[M];
 	int64_t y2[M];	
	int64_t xy[M];
	product_in_ring(x, x, x2, true);
	product_in_ring(y, y, y2, true);
	product_in_ring(x, y, xy, true);
	
	identity_ring_matrix(S, E);
	
	for(int k=0; k<M; k++)
    {
    	E[0][0][k] -= xy[k];
    	E[0][1][k] += y[k];
    	E[0][2][k] -= y2[k];
    	
    	E[1][0][k] -= x[k];
    	E[1][2][k] -= y[k];
    	
    	E[2][0][k] += x2[k];
    	E[2][1][k] -= x[k];
    	E[2][2][k] += xy[k];
	}	
	
	rmm_multiply(S, S, S, E, T, B22inv);
	
	free_ring_matrix(S, S, E); E = NULL;
	free_ring_matrix(S, S, EP); EP = NULL;
    free_ring_matrix(S, S, T); T = NULL;
}

 /**
 * This function generates the submatrix B_21 of matrix B based on the construction described in the paper.
 *
 * @param B21   A pointer to the 3D matrix that will store the matrix B_21.
 */
void generate_B21(int64_t** B21)
{
	zero_ring_vector(S, B21);
	
	int idx[M];
 	range_vector(M, idx);
	
	for(int i=0; i<S; i++)
	{
		permute_vector(LB, M, idx);
		
		for(int k=0; k<LB; k++)
	    	B21[i][idx[k]] = rng4();
	}
}

/**
 * This function verifies that all coefficients in the signature vector y are within the allowable range specified by Y_BOUND.
 * If any coefficient exceeds this bound, the signature is deemed invalid.
 *
 * @param y    A pointer to the 2D matrix containing the coefficients to be validated.
 * @return     A boolean value; true if all coefficients in y are within the specified bound, false if any exceed Y_BOUND.
 */
bool valid_y(int64_t** y)
{
	for(int i=0; i<S; i++)
		for(int k=0; k<M; k++)
			if(abs(y[i][k]) >= Y_BOUND)
				return false;

	return true;
}

 /**
 * This function encodes the signature vector y into a compact binary format and concatenates it with the original message m.
 * The coefficients of the matrix y are offset by Y_BOUND and then encoded as binary bits with a specified number of bits per coefficient (Y_BITS).
 * These bits are packed into bytes, and these bytes are stored in the array sm. The original message m is then appended to this encoded data,
 * forming a single combined array. The function calculates the total length of this combined array and stores it in smlen.
 *
 * @param m       A pointer to the original message array.
 * @param mlen    The length of the original message.
 * @param y       A pointer to the vector y containing the signature coefficients.
 * @param sm      A pointer to the signed message.
 * @param smlen   A pointer to the length of the signed message.
 */
void my_to_sm(const unsigned char* m, unsigned long long mlen, int64_t** y, unsigned char* sm, unsigned long long* smlen)
{
	int total_bits = Y_BITS*M*S;
	bool bits[total_bits];
	int biti = 0;
	
	for(int i=0; i<S; i++)
	{
		for(int k=0; k<M; k++)
		{
			int64_t val = y[i][k] + Y_BOUND;
		
			for(int bit=0; bit<Y_BITS; bit++)
        		bits[biti++] = (val >> bit) & 1;
		}
	}
	
	int smi = 0;
	for(int i=0; i<total_bits; i+=8)
	{
		unsigned char val = 0;
		
		for(int j=0; j<8; j++)
			val += bits[i+j] << j;
				
		sm[smi++] = val;
	}
	
	*smlen = mlen + smi;
	
	for(int i=0; i<mlen; i++)
		sm[smi+i] = m[i];
}

/**
 * This function generates the DEFI signature for a given message.
 * The signature is encoded in the 'sm' array, and the length of the signature is stored in 'smlen'.
 *
 * @param sm      A pointer to the signed message.
 * @param smlen   A pointer to the length of the signed message.
 * @param m       A pointer to the message.
 * @param mlen    The length of the message.
 * @param sk      A pointer to the secret key.
 * @return 0 for successful execution
 */
int sig_gen(unsigned char *sm, unsigned long long *smlen, const unsigned char *m, unsigned long long mlen, const unsigned char *sk)
{
	initialize_rng((unsigned char*)sk, 48);
	
	int64_t*** B22inv = allocate_ring_matrix(S, S);
	generate_B22inv(B22inv);
	
	int64_t** B21 = allocate_ring_vector(S);
	generate_B21(B21);
	
	int64_t h[M]; // is also Z1
	hash_of_message(m, mlen, h);
	
	int64_t** B21h = allocate_ring_vector(S);
	for(int i=0; i<S; i++)
		product_in_ring(B21[i], h, B21h[i], true);
	
	
	// This is to provide a deterministic version. Ideally, the rng is initialised using some entropy source.
	//.....................................//
	unsigned char new_seed[48];
	FIPS202_SHAKE256(m, mlen, new_seed, 48);
	for(int i=0; i<48; i++)
		new_seed[i] += sk[i];
	initialize_rng(new_seed, 48);
	//.....................................//
	
	
	int64_t u1[M];
	int64_t u2[M];
	int64_t v[M];
	int64_t v2[M];
	int64_t u1u1[M];
	int64_t t[M];
	int64_t u2u1u1[M];
	int64_t Z1u1[M];
	int64_t u1u2[M];
	
	int idx[M];
 	range_vector(M, idx);
	
	int64_t** y = allocate_ring_vector(S);
	
	do
	{
		zero_vector(M, u1);
		zero_vector(M, u2);
		zero_vector(M, v2);
		
		permute_vector(LU, M, idx);
		for(int k=0; k<LU; k++)
			u1[idx[k]] = rng4();
		
		permute_vector(LU, M, idx);
		for(int k=0; k<LU; k++)
		{
			v2[idx[k]] = rng4();
			u2[idx[k]] = 2*v2[idx[k]];
		}
		
		product_in_ring(u1, u1, u1u1, true);
	
		zero_vector(M, t);
		t[0] = 1;
		for(int k=0; k<M; k++)
			t[k] -= u1u1[k];
			
		product_in_ring(v2, t, v, true);
		product_in_ring(u2, u1u1, u2u1u1, true);
		product_in_ring(h, u1, Z1u1, true);
		product_in_ring(u1, u2, u1u2, true);
		
		for(int k=0; k<M; k++) // Z" - B21h >>> stored in B21
		{
			B21[0][k] = v[k] + u2u1u1[k] - Z1u1[k] - B21h[0][k];
			B21[1][k] = v[k] + Z1u1[k]- B21h[1][k];
			B21[2][k] = u1u2[k] - h[k] - B21h[2][k];
		}
		
		rmv_multiply(S, S, B22inv, B21, y);
	}
	while(valid_y(y)==false);
	
	clear_rng();
	
	my_to_sm(m, mlen, y, sm, smlen);

	free_ring_matrix(S, S, B22inv); B22inv = NULL;
	free_ring_vector(S, B21); B21 = NULL;
	free_ring_vector(S, B21h); B21h = NULL;
	free_ring_vector(S, y); y = NULL;

	return 0;
}



