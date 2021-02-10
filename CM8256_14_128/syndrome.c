#include "syndrome.h"
#include "evaluation.h"
#include "params.h"

/* input poly: the monic and irreducible polynomial */
/*       list: n elements in Fq                     */
/*       word: C0 append k zeros(v)                 */
/* output synd: syndrome of length 2t */
void syndrome(fq *synd, fq *poly, fq *list, unsigned char *word)
{
    int i,j,k = 0;
    fq poly_alpha, inverse;
    unsigned char bit;

    for (i = 0; i < 2*GOPPA_T; ++i)
        synd[i] = 0;

    // calculate synd[j] = ∑ word[i] / (poly[a])^2 * a^j, i:1 to n
    for (i = 0; i < GOPPA_N; ++i)
    {
        bit = word[i >> 3] & (1 << (i & 7));

        poly_alpha = evaluation(poly, list[i]);
        inverse = fq_inv(fq_mul(poly_alpha,poly_alpha));

        // synd[j] = synd[j] + (list[i])^j
        for (j = 0; j < 2*GOPPA_T; ++j)
        {
            synd[j] = (bit ? fq_add(synd[j], inverse) : synd[j]);
            inverse = fq_mul(inverse, list[i]);
        }
        
    }
}