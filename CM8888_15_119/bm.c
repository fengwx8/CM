/* 
    This file refers to the following paper to implements the Berlekamp-Massey algorithm
    "Implementation of Berlekamp-Massey algorithm without inversion" Xu Youzhi 
*/
#include "bm.h"
#include <stdbool.h>

/* input list: a list of field elements of size 2*GOPPA_T 
   output poly: */
void berlekamp_massey(fq *poly, fq *list){
    int i;

    uint16_t k = 0;
    uint16_t lk = 0;

    // dk = 0 iff dk_EQU_0 is true
    // k < 2lk iff k_LSS_2lk is true
    bool dk_EQU_0, k_LSS_2lk;
    /* cond is true iff dk_EQU_0 is true 
        or k_LSS_2lk is true, otherwise cond is false*/ 
    // in other words, cond = dk_EQU_0 | k_LSS_2lk
    bool cond;

    fq beta[ GOPPA_T+1 ], sigma[ GOPPA_T+1 ], temp[ GOPPA_T+1 ];

    fq delta = 1, dk;

    // initialize β(x) = x, σ(x) = 1;
    for (i = 0; i <= GOPPA_T ; ++i)
        beta[i] = sigma[i] = 0;
    
    beta[1] = sigma[0] = 1;

    for (k = 0; k < 2 * GOPPA_T; ++k)
    {
        dk = 0;
        //d = ∑ σi * S[k-i] , i=[0,min(k,t)]
        int bounds = (k < GOPPA_T ? k : GOPPA_T);
        for (i = 0; i <= bounds; ++i)
            dk ^= fq_mul(sigma[i], list[k-i]);

        dk_EQU_0 = (dk ? false : true);
        k_LSS_2lk = (k < 2*lk ?  true: false);

        cond = dk_EQU_0 | k_LSS_2lk;

        for (i = 0; i <= GOPPA_T; ++i)
            temp[i] = beta[i];

        // update β(x)
        for (i = 0; i <= GOPPA_T; ++i)
            beta[i] = cond ? beta[i] : sigma[i];
        // β(x) = xβ(x)
        for (i = GOPPA_T; i >= 1; --i)
            beta[i] = beta[i-1];
        beta[0] = 0;

        for (i = 0; i <= GOPPA_T; ++i)
		{
			sigma[i] = fq_mul(delta, sigma[i]);
			sigma[i] ^= fq_mul(dk, temp[i]);
		} 

        // update lk
        lk = cond ? lk : k+1-lk;

        // update delta
        delta = cond ? delta : dk;
    }

    for (i = 0; i <= GOPPA_T; ++i)
        poly[i] = sigma[ GOPPA_T-i ];
}
