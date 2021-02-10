#include "permutation.h"

void permutation(uint16_t *pi, unsigned char* cbits, int m)
{   
    int i,j,gap,pos;
    int n = 1 << m;
    int index = 0;
    unsigned char bit_and[8] = {1,2,4,8,16,32,64,128};
    uint16_t swap;
    
    // for (i = 0; i < n; ++i)
    //     pi[i] = i;
    
    for (i = 0; i < 2*m-1; ++i)
    {
        pos = i < 2*m-2-i ? i : 2*m-2-i;
        gap = 1 << pos;
        for (j = 0; j < n/2; ++j)
        {
            if(cbits[index >> 3] & bit_and[index & 7])
            {
                pos = (j%gap) + 2*gap*(j/gap);
                swap = pi[pos] ^ pi[pos + gap];
                pi[pos] ^= swap;
                pi[pos + gap] ^= swap;
            }
            ++index;
        }
    }
}