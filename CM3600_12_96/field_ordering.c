#include "field_ordering.h"
#include "fq.h"

// min(a,b) = (a^(b^a)&(-(a<b)))
// max(a,b) = (b^(b^a)&(-(a<b)))
// bit shift is needed
#define MINMAX64(a,b) \
do{ \
    uint64_t c = b - a; \
    c >>= 63; \
    c = -c; \
    c &= a ^ b; \
    a ^= c; \
    b ^= c; \
}while(0)

void sort64(uint64_t *e,int num)
{
    int bound, tp = 1;
    int i,j,m,n;
    do{
        bound = tp;
        tp <<= 1;
    } while (tp < num);

    for (i = bound; i > 0; i >>= 1)
    {
        for (j = 0; j < num - i; ++j)
            if(!(j & i))
                MINMAX64(e[j], e[i+j]);
        j = 0;

        for (m = bound; m > i; m >>= 1)
        {
            for (; j < num - m; ++j)
            {
                if(!(j & i)){
                    uint64_t temp = e[i + j];
                    for (n = m; n > i; n >>= 1)
                        MINMAX64(temp, e[j+n]);
                    e[i + j] = temp;
                }
            }
        }
    }
}


/* input a: q 4-byte integers      */
/* output pi: permutation          */
int field_ordering(uint16_t *pi, uint32_t *a)
{
    int i;
    uint64_t sort[Q];
    uint16_t bit_and = Q - 1;

    // step 2, if not distinct, return ⊥
    for(i = 1; i < Q; ++i)
        if(a[i] == a[i-1])
            return -1;

    // step 1, reverse the bit 
    for(i = 0; i < Q; ++i)
    {
        sort[i] = (uint64_t)a[i];
        sort[i] <<= 31;
        sort[i] |= i;
    }

    // step 3, sorting
    sort64(sort, Q);

    // step 4, generate α1...αq
    for(i = 0; i < Q; ++i)
        pi[i] = sort[i] & bit_and;
        
        
    return 0;
}