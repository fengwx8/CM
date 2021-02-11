#include <string.h>
#include <stdio.h>

#include "pk_gen.h"
#include "evaluation.h"
#include "little_endian.h"

int search_pivot(uint32_t row)
{
    int i;
    for (i = 31; i >= 0; --i)
    {
        if ((row >> i) & 1)
            break;
    }

    return 31 - i;
}

int column_switch(unsigned char matrix[][GOPPA_N/8], uint16_t *pi, uint32_t *col_pi)
{
    int i, j, index, row, len, tmp;
    unsigned char row_block[4];
    uint32_t mat[16]; // store 16*32 matrix
    uint32_t check, temp, ex;
    int pivot_ind[16];
    uint16_t exchange;
    
    row = PK_ROW - 16;
    index = row / 8;
    len = row & 7;

    // get the 16*32 mat to verify the matrix is semi-systematic
    for (i = 0; i < 16; ++i)
    {
        mat[i] = 0;
        for (j = 0; j < 4; ++j)
        {
            row_block[ j ] = matrix[ row + i ][ index + j ];
            mat[i] <<= 8;
            mat[i] |= byterev(row_block[ j ]);
        }
    }

    *col_pi = 0;

    // get the index of row pivot
    // ensure the index is within 0 to 31 
    // and Meet the qualification ci < ci+1
    // the indices of row pivot are stored in col_pi
    for (i = 0; i < 16; ++i)
    {
        check = mat[i];
        for (j = i+1; j < 16; ++j)
            check |= mat[j];
        
        if (check == 0)
            return 1;
        
        pivot_ind[i] = search_pivot(check);

        tmp = 31 - pivot_ind[i];

        (*col_pi) |= (1 << pivot_ind[i]);

        j = i + 1;
        while (!((mat[i] >> tmp) & 1) && (j < 16))
            mat[i] ^= mat[j++];

        for (j = i+1; j < 16; ++j)
        {
            if ((mat[j] >> tmp) & 1)
                mat[j] ^= mat[i];
        } 
    }

    // change permutation
    for (i = 0; i < 16; ++i)
    {
        exchange = pi[ row + i ] ^ pi[ row + pivot_ind[i] ];
        pi[ row + i ] ^= exchange;
        pi[ row + pivot_ind[i] ] ^= exchange;
    }

    // move column of matrix according to pivot_ind
    for (i = 0; i < PK_ROW; ++i)
    {   

        for (j = 0; j < 4; ++j)
            row_block[j] = matrix[i][ index + j ];
        
        check = load4(row_block);

        for (j = 0; j < 16; ++j)
        {
            temp  = check >> j;
			temp ^= check >> pivot_ind[ j ];
			temp &= 1;
        
			check ^= temp << pivot_ind[ j ];
			check ^= temp << j;
        }

        store4(&matrix[i][index], check);
    }
	
	printf("column swap(index begin with 0):");
	for (i = 0; i < 16; ++i)
		printf(" (%d, %d)", row + i, row + pivot_ind[i]);
	printf("\n\n");

    return 0;
}

/* input pi: permutation                    */
/*       poly: monic irreducible polynomial */
/* output pk: public key(matrix T)          */
/* return 0 for success, 1 for failure      */
int pk_generation(unsigned char *pk, fq *pi, fq* poly, uint32_t *col_pi)
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

            if (row == PK_ROW - 16)
            {
                if (column_switch(matrix, pi, col_pi))
                    return 1;
            }

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
