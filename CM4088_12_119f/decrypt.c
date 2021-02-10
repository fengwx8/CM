#include <stdio.h>

#include "decrypt.h"

#include "bm.h"
#include "evaluation.h"
#include "support.h"
#include "syndrome.h"
#include "little_endian.h"

/* input sk: secret key                */
/*       rword: C0                     */
/* output vector: the error e          */
/* return 0 for success, 1 for failure */
int decrypt(unsigned char* vector, unsigned char* sk, unsigned char* rword)
{
    unsigned char v[ GOPPA_N/8 ];
    
    fq g[ GOPPA_T+1 ]; // polynomial
    fq S[ GOPPA_N ]; // support
    fq synd[ 2*GOPPA_T ]; // syndrome
    fq locator[ GOPPA_T+1 ]; // error-locator polynomial
    fq root[ GOPPA_N ]; // the root of error-locator polynomial
    //// root[i] = 0 means e[i] = 1;

    fq synd_e[ 2*GOPPA_T ]; // syndrome for check
    int i = 0, weight = 0;

    // v = C0 append k zeros, k = n - m*t
    for (; i < GOPPA_N/8; ++i)
        v[i] = ((i < CIPHERBYTE) ? rword[i] : 0);

    // a monic and irreducible polynomial of degree t
    // extract from sk
    g[ GOPPA_T ] = 1;
    for (i = 0; i < GOPPA_T; ++i){
        g[i] = load2(sk);
        sk += 2;
    }

    // extract n field elements
    support(S, sk);
    // compute syndrome
    syndrome(synd, g, S, v);
    // get the error-locator polynomial
    berlekamp_massey(locator, synd);

    for (i = 0; i < GOPPA_N; ++i)
        root[i] = evaluation(locator, S[i]);

    for (i = 0; i < GOPPA_N/8; ++i)
        vector[i] = 0;

    printf("error e in decrypt: ");

    for (i = 0; i < GOPPA_N; ++i) {
        if (root[i] == 0)
        {
            vector[i >> 3] |= (1 << (i & 7));
            printf("%x ", i);
            weight++;
        }
    }
    printf("\n");

    // if weight != t, failed!
    if (weight ^ GOPPA_T) return 1;

    syndrome(synd_e, g, S, vector);

    // verify Hv = He
    for (i = 0; i < 2*GOPPA_T; ++i) 
        if (synd_e[i] ^ synd[i]) 
            return 1; 
    
    return 0;
}