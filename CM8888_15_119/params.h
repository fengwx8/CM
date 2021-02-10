#ifndef PARAMS_H
#define PARAMS_H

#include <stdint.h>

// in this test, GOPPA_N is in range [8192,16384)
// GOPPA_N also must be a multiple of 8
#define GOPPA_N 8888 
#define GOPPA_T 119
#define M 15
#define Q (1 << M)
#define MU 16
#define NU 32

#define sigma1 16
#define sigma2 32

#define POLY_BYTE 2*GOPPA_T
#define CONTROLBYTES (2*M - 1)*(1 << (M - 4))

#define PK_ROW (M * GOPPA_T)
#define PK_COLUMN (GOPPA_N - PK_ROW)
#define PK_ROW_BYTE ((PK_COLUMN + 7) / 8)

#define PK_SIZE (PK_ROW * PK_ROW_BYTE)
#define SK_SIZE (32 + 4 + 2*GOPPA_T + CONTROLBYTES + GOPPA_N/8)
#define HASHBYTE 32
#define CIPHERBYTE ((PK_ROW + 7) / 8)
#define CIPHERTEXTBYTE (32 + CIPHERBYTE)

typedef uint16_t fq;

#endif
