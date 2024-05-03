#ifndef rng_functions_h
#define rng_functions_h

#include <stdbool.h>

void refill_rng_buffer();
bool rng_bit();
unsigned char rng_byte();
int rng(int x);
int rng2();
int rng4();
void initialize_rng(unsigned char* seed, int seed_size);
void clear_rng();

#endif
