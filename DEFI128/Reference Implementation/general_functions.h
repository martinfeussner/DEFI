#ifndef general_functions_h
#define general_functions_h

#include <stdbool.h>

extern const int P[6][3];

int64_t** allocate_ring_vector(int n);
int64_t*** allocate_ring_matrix(int m, int n);
void free_ring_vector(int n, int64_t** A);
void free_ring_matrix(int m, int n, int64_t*** A);
void zero_vector(int n, int64_t* A);
void zero_ring_vector(int n, int64_t** A);
void zero_ring_matrix(int m, int n, int64_t*** A);
void identity_ring_matrix(int n, int64_t*** A);
void copy_ring_matrix(int m, int n, int64_t*** A, int64_t*** B);
void range_vector(int n, int* A);
void permute_vector(int m, int n, int* A);
void row_permute(int idx, int64_t*** E, int64_t*** EP);
void col_permute(int idx, int64_t*** E, int64_t*** EP);
void product_in_ring(int64_t* poly1, int64_t* poly2, int64_t* result_poly, bool overwrite);
void rmv_multiply(int m, int l, int64_t*** A, int64_t** b, int64_t** c);
void rmm_multiply(int m, int l, int n, int64_t*** A, int64_t*** B, int64_t*** C);
void hash_of_message(const unsigned char* m, unsigned long long mlen, int64_t* h);

#endif
