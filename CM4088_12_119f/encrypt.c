#include <stdio.h>

#include "encrypt.h"
#include "params.h"
#include "randnum.h"
#include "little_endian.h"

void vector_generation(unsigned char *e)
{
    unsigned char bytes[ 2*GOPPA_T * sizeof(uint16_t) ];
    int count, i, j, pass;
    fq temp, mod = (fq)(Q - 1);
    unsigned char byte;


    uint16_t elements[ 2*GOPPA_T ]; // save 2*t numbers
    fq index[ GOPPA_T ]; // save index

    while (1)
    {
        randomseq(bytes, sizeof(bytes));

        for (i = 0; i < GOPPA_T*2; ++i)
            elements[i] = load2(bytes + 2*i) & mod;

        count = i = 0; 

        while (count < GOPPA_T && i < 2*GOPPA_T)
        {
            if (elements[i] < GOPPA_N)
                index[ count++ ] = elements[i];
            i++;
        }

        if (count < GOPPA_T) continue;

        // if not all distinct, continue
        pass = 1;
        for (i = 1; i < GOPPA_T; ++i)
            for (j = 0; j < i; ++j)
                if (index[i] == index[j])
                    pass = 0;

        if (pass) break;
    }

    printf("error e in encrypt: ");

    // corresponding elements are assigned to 1
    for (j = 0; j < GOPPA_T; ++j)
    {
        temp = index[j];
        byte = 1 << (temp & 7);
        e[temp >> 3] |= byte;
        printf("%x ", temp);
    }

    printf("\n");
}


/* input pk: public key matrix T                 */
/* output c0: result after matrix multiplication */
/*        e: fixed-weight-vector                 */
void encrypt(unsigned char *c0, unsigned char *pk, unsigned char *e)
{
    unsigned char row[ GOPPA_N/8 ];
    unsigned char bit = 0;
    int i, j, len = PK_ROW & 7;

    // generate fixed-weight-vector 
    vector_generation(e);

    for (i = 0; i < CIPHERBYTE; ++i)
        c0[i] = 0;

    // matrix multiplication
    for (i = 0; i < PK_ROW; ++i)
    {
        for (j = 0; j < GOPPA_N/8; ++j)
            row[j] = 0;

        for (j = 0; j < PK_ROW_BYTE; ++j)
            row[ GOPPA_N/8 - PK_ROW_BYTE + j ] = pk[j];

        for (j = GOPPA_N/8-1; j >= GOPPA_N/8 - PK_ROW_BYTE; j--) 
			row[ j ] = (row[ j ] << len) | (row[j-1] >> (8-len));

        row[ i/8 ] |= 1 << (i & 7);
        
        bit = 0;
        for (j = 0; j < GOPPA_N/8; ++j)
            bit ^= (row[ j ] & e[ j ]);

        for (j = 4; j >=1 ; j >>= 1){
            bit ^= (bit >> j);
        }
        bit &= 1;

        c0[ i/8 ] |= (bit << (i & 7));

        pk += PK_ROW_BYTE;
    }
}