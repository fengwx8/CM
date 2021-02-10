#ifndef SUPPORT_H
#define SUPPORT_H

#include "permutation.h"
#include "fq.h"

/* input cbits: controlbits */
/* output s: the first n elements after permutation and bit reverse */
void support(fq *s, unsigned char* cbits){
    uint16_t pi[Q];
    int i;
    for (i = 0;i < Q; ++i)
        pi[i] = elementrev((fq)i);
    
    permutation(pi,cbits,M);

    for (i = 0; i < GOPPA_N; ++i)
        s[i] = pi[i];
}

#endif