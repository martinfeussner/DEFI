#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include "rng.h"

// Defines the buffer size and buffer holding the random values for the RNG.
#define RNG_BUFFER_SIZE 1024 
unsigned char rng_buffer[RNG_BUFFER_SIZE];
int rng_buffer_idx = RNG_BUFFER_SIZE;  // Trigger refill on first use.
int current_bit_idx = 8;  // Start at 8 to trigger byte fetch on first use.
unsigned char current_byte;  // Holds the byte from which bits are being extracted.

void refill_rng_buffer()
{
    randombytes(rng_buffer, RNG_BUFFER_SIZE);
    rng_buffer_idx = 0;
}

bool rng_bit()
{
    if(current_bit_idx == 8)  // All bits in the current byte are used
    {
        if(rng_buffer_idx == RNG_BUFFER_SIZE)  // If buffer is empty, refill
        	refill_rng_buffer();
        	
        current_byte = rng_buffer[rng_buffer_idx];
        current_bit_idx = 0;  // Reset bit index for the new byte
        rng_buffer_idx++;     // Move to the next byte
    }

    bool bit = (current_byte & (1 << current_bit_idx)) != 0;
    current_bit_idx++;

    return bit;
}

unsigned char rng_byte()
{
    if(rng_buffer_idx == RNG_BUFFER_SIZE)  // If buffer is empty, refill
    	refill_rng_buffer();
    	
    return rng_buffer[rng_buffer_idx++];
}

int rng(int x)
{
    if(rng_buffer_idx + sizeof(uint16_t) > RNG_BUFFER_SIZE)
    	refill_rng_buffer();

    uint16_t random_value;
    memcpy(&random_value, &rng_buffer[rng_buffer_idx], sizeof(uint16_t));
    rng_buffer_idx += sizeof(uint16_t);
    
    return random_value % x;
}

int rng2()
{
	return (rng_bit()==0)?1:-1;
}

int rng4()
{
	int a = (rng_bit()==0)?1:-1;
	int b = (rng_bit()==0)?1:2;
	
	return a*b;
}

void initialize_rng(unsigned char* seed, int seed_size)
{
	randombytes_init(seed, NULL, 8*seed_size);
    rng_buffer_idx = RNG_BUFFER_SIZE;
	current_bit_idx = 8;
}

void clear_rng()
{
	for(int i=0; i<RNG_BUFFER_SIZE; i++)
		rng_buffer[i] = 0;
	
	current_byte = 0;
}



