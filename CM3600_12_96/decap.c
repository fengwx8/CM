#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "decrypt.h"
#include "params.h"
#include "shake256.h"

void decap()
//int main()
{
    FILE *fp;
    unsigned char sk[ SK_SIZE ]; // secret key
    unsigned char content[ 1 + GOPPA_N/8 + CIPHERTEXTBYTE] = {0}; // b = 0
    unsigned char vector_e[ 1 + GOPPA_N/8 ] = {2};
    unsigned char *s = content + 1, *ciptext = content + 1 + GOPPA_N/8;
    unsigned char *skp = sk;
    unsigned char c1[ HASHBYTE ];
    unsigned char key[ HASHBYTE ];

    int i, different;

    if ((fp = fopen("sk.txt","r")) == NULL)
    {
        printf("nerror in open sk.txt.\n");
        exit(1);
    }
    for (i = 0; i < SK_SIZE; ++i)
        fscanf(fp, "%02hhx", &sk[i]);
    fclose(fp);

    if ((fp = fopen("ciphertext.txt","r")) == NULL)
    {
        printf("nerror in open ciphertext.txt.\n");
        exit(1);
    }
    for (i = 0; i < CIPHERTEXTBYTE; ++i)
        fscanf(fp, "%02hhx", &ciptext[i]);
    fclose(fp);

    // 3.extract s from sk
    memcpy(s, sk + 36 + POLY_BYTE + CONTROLBYTES, GOPPA_N/8);

    // 4.compute e by decrypt
    different = decrypt(vector_e + 1, sk + 36, ciptext);

    // 5.compute C1’ = hash(2,e);
    hash(c1, vector_e, sizeof(vector_e));

    // 0 for b = 1, 1 for b = 0
    // 1.ciptext + CIPHERBYTE store C1 in encrypt
    // c1 below means C1 after hash, that is C1’
    different |= memcmp(c1, ciptext + CIPHERBYTE, HASHBYTE);

    // 4、6.judge and determine whether the value of b changes
    // different = 0 for preimage = (1 | e | C)
    // otherwise preimage = (0 | s | C)
    if (!different) {
        // set b = 1
        vector_e[0] = 1;
        memcpy(content, vector_e, sizeof(vector_e));
    }

    // 7.key K = hash(b,e,C)
    hash(key, content, sizeof(content));

    // 8.output K
    printf("session key calculated after decrypt:\n");
    for (i = 0; i < HASHBYTE; ++i)
        printf("%02hhx", key[i]);

    printf("\n\n\n");

    //return 0;
}