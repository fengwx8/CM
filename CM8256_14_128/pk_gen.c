#include <string.h>

#include "pk_gen.h"
#include "evaluation.h"
#include "little_endian.h"

/* input pi: permutation                    */
/*       poly: monic irreducible polynomial */
/* output pk: public key(matrix T)          */
/* return 0 for success, 1 for failure      */
int pk_generation(unsigned char *pk, fq *pi, fq* poly)
{
    fq inv[ GOPPA_N ];
    fq g[ GOPPA_T+1 ];
    fq alpha[ GOPPA_N ];
    unsigned char matrix[ PK_ROW ][ GOPPA_N/8 ];

    int i,j,k,c,row = 0;
    unsigned char byte, bit;

    // monic irreducible polynomial
    for (i = 0; i < GOPPA_T; ++i)
        g[i] = poly[i];
    g[ GOPPA_T ] = 1;

    // get the 1/g(a) matrix
    for (i = 0; i < GOPPA_N; ++i){
        alpha[i] = elementrev(pi[i]);
        inv[i] = fq_inv(evaluation(g, alpha[i]));
    }
        

    for (i = 0;i < PK_ROW; ++i)
        for (j = 0; j < GOPPA_N/8; ++j)
            matrix[i][j] = 0;

    for (i = 0; i < GOPPA_T; ++i)
    {
        for (j = 0; j < GOPPA_N; j += 8)
        {
            for (k = 0; k < M; ++k)
            {
                byte  = (inv[j+7] >> k) & 1; byte <<= 1;
                byte |= (inv[j+6] >> k) & 1; byte <<= 1;
                byte |= (inv[j+5] >> k) & 1; byte <<= 1;
                byte |= (inv[j+4] >> k) & 1; byte <<= 1;
                byte |= (inv[j+3] >> k) & 1; byte <<= 1;
                byte |= (inv[j+2] >> k) & 1; byte <<= 1;
                byte |= (inv[j+1] >> k) & 1; byte <<= 1;
                byte |= (inv[j+0] >> k) & 1;
                matrix[ i*M + k ][ j/8 ] = byte;
            }
        }

        for (j = 0; j < GOPPA_N; ++j)
            inv[j] = fq_mul(inv[j], alpha[j]);
    }

    // gaussian elimination
    for (i = 0; i < (PK_ROW+7) / 8; ++i) {
        for (j = 0; j < 8; ++j) {

            // PK_ROW % 8 != 0
            if (row >= PK_ROW) break;

            // ensure the pivot is nonzero
            k = row + 1;
            while (!(bit = (matrix[ row ][ i ] >> j) & 1) && (k < PK_ROW)) {
                for (c = 0; c < GOPPA_N/8; ++c)
                    matrix[ row ][ c ] ^= matrix[ k ][ c ];
                ++k;
            }

            if (!bit) return 1; 
            
            // elementary row transformation
            for (k = 0; k < PK_ROW; ++k) {
                if (k == row || !((matrix[ k ][ i ] >> j) & 1)) continue;
                
                for (c = 0; c < GOPPA_N/8; ++c)
                    matrix[ k ][ c ] ^= matrix[ row ][ c ];
            }

            row++;
        }

    }

    c = (PK_ROW + 7) / 8;
    // save pk
    for (i = 0; i < PK_ROW; ++i) 
        memcpy(pk + i*PK_ROW_BYTE, matrix[i] + c, PK_ROW_BYTE);

    return 0;
}
