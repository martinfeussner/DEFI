#include <stdio.h>
#include <stdlib.h>

#include "defi_keygen.h"
#include "defi_siggen.h"
#include "defi_sigver.h"

#define CRYPTO_SECRETKEYBYTES 48
#define CRYPTO_PUBLICKEYBYTES 800
#define CRYPTO_BYTES 432

int crypto_sign_keypair(unsigned char *pk, unsigned char *sk)
{
	return key_gen(pk, sk);
}

int crypto_sign(unsigned char *sm, unsigned long long *smlen, const unsigned char *m, unsigned long long mlen, const unsigned char *sk)
{
	return sig_gen(sm, smlen, m, mlen, sk);
}

int crypto_sign_open(unsigned char *m, unsigned long long *mlen, const unsigned char *sm, unsigned long long smlen, const unsigned char *pk)
{
	return sig_ver(m, mlen, sm, smlen, pk);
}

