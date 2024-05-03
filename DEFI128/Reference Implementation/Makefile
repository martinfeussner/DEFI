CC = /usr/bin/gcc
CFLAGS = -g -O3 -std=c99
LDFLAGS = -static-libgcc -lssl -lcrypto

SOURCES = sign.c defi_keygen.c defi_siggen.c defi_sigver.c keccak.c rng.c rng_functions.c general_functions.c PQCgenKAT_sign.c
HEADERS = api.h parameters.h defi_keygen.h defi_siggen.h defi_sigver.h keccak.h rng.h rng_functions.h general_functions.h

PQCgenKAT_sign: $(HEADERS) $(SOURCES)
	$(CC) $(CFLAGS) -o $@ $(SOURCES) $(LDFLAGS)

.PHONY: clean

clean:
	-rm PQCgenKAT_sign