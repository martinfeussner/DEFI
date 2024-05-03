#ifndef parameters_h
#define parameters_h

// Parameters from the description

#define SECURITY 128 // security in bits

#define N 4 // $n$
#define S 3 // $s$
#define R 1 // $r$
#define M 64 // $m$
#define K 14 // $k$
#define LB 24 // $\lambda_{B_{21}}$
#define LE 1 // $\lambda_{E}$ --- NOT USED
#define LF 17 // $\lambda_{F}$
#define LU 35 // $\lambda_{F}$
#define SF 6 // Computed: $s!$
#define G 106 // Computed: Computed using $m$ and $\lambda_{B_{21}}$ - Base 2 of the guessing complexity for entries of $B_{21}$ and $B_{22}$

#define C1_BOUND 64 // $\delta_{C_1}$
#define C2_BOUND 256 // $\delta_{C_2}$
#define C3_BOUND 1024 // $\delta_{C_3}$
#define Y_BOUND 131072 // $\delta_{C_4}$
#define C1_BITS 7 // Computed: $log_2(2*\delta_{C_1})$
#define C2_BITS 9 // Computed: $log_2(2*\delta_{C_2})$
#define C3_BITS 11 // Computed: $log_2(2*\delta_{C_3})$
#define Y_BITS 18 // Computed: $log_2(2*\delta_{y})$

#endif
