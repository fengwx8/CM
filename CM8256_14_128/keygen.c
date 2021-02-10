#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "fq.h"
#include "little_endian.h"
#include "poly_gen.h"
#include "field_ordering.h"
#include "controlbits.h"
#include "pk_gen.h"
#include "randnum.h"
#include "shake256.h"

void keygen()
//int main() 
{
    unsigned char pk[ PK_SIZE ];
    unsigned char sk[ SK_SIZE ];
    int eg_size = (GOPPA_N + sigma2*Q + sigma1*GOPPA_T + 256) / 8;
    int i;

    unsigned char *random_pointer, *sk_pointer;
    unsigned char seed[ HASHBYTE+1 ] = {64};
    unsigned char eg[ eg_size ];

    

    fq beta[ GOPPA_T ]; // an element in fq(2^mt)
    fq poly[ GOPPA_T ]; // monic irreducible polynomial
    uint32_t a[ Q ]; // q numbers to generate permutation pi
    fq pi[ Q ]; // permutation pi

    FILE *pkFP, *skFP;

    randomseq(seed+1, 32);

    while(1)
    {
        // store for the new σ
        random_pointer = eg + eg_size - 32;
        sk_pointer = sk;

        // E = G(σ)
        DRBGshake(eg, eg_size, seed, 33);

        // save σ in sk
        memcpy(sk_pointer, seed+1, HASHBYTE);

        sk_pointer += (HASHBYTE + 4);
        // update σ
        memcpy(seed+1, random_pointer, HASHBYTE);

        // generating irreducible polynomial
        for (i = GOPPA_T-1; i >= 0; --i){
            random_pointer -= 2;
            beta[i] = load2(random_pointer);
        }

        if(poly_generation(poly, beta))
            continue;

        for (i = 0; i < GOPPA_T; ++i){
            store2(sk_pointer, poly[i]);
            sk_pointer += 2;
        }
        
        // field-ordering algorithm
		// generating permutation
        for (i = Q-1; i >= 0; --i){
            random_pointer -= 4;
            a[i] = load4(random_pointer);
        }
        
        if (field_ordering(pi, a)) continue;

        // generate public key pk
        if (pk_generation(pk, pi, poly)) continue;

        // generate controlbits
        controlbits(sk_pointer, pi, M, Q);
        sk_pointer += CONTROLBYTES;

        // store the string s
        memcpy(sk_pointer, eg, GOPPA_N/8);

        store4(sk + HASHBYTE, 0xffff);

        break;
    }

    pkFP = fopen("pk.txt","w");
    skFP = fopen("sk.txt","w");

    if ((pkFP == NULL) || (skFP == NULL))
    {
        printf("Can not save keys.\n");
        exit(1);
    }

    for (i = 0; i < PK_SIZE; ++i)
        fprintf(pkFP, "%02x", pk[i]);

    for (i = 0; i < SK_SIZE; ++i)
        fprintf(skFP, "%02x", sk[i]);

    fclose(pkFP);
    fclose(skFP);
    
    //free(pk);
    //free(sk);

    //return 0;
}
