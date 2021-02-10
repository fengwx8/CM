#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "encrypt.h"
#include "params.h"
#include "shake256.h"

void encap()
//int main()
{
    FILE *fp;
    unsigned char content[ 1 + GOPPA_N/8 + CIPHERTEXTBYTE ] = {2};
    unsigned char *c = content + 1 + GOPPA_N/8;// c = (C0,C1)
    unsigned char *c1 = c + CIPHERBYTE;
    unsigned char *e = content + 1;
    unsigned char pk[ PK_SIZE ];
    unsigned char key[HASHBYTE];

    int i;

    if ((fp = fopen("pk.txt","r")) == NULL)
    {
        printf("nerror in open pk.txt.\n");
        exit(1);
    }

    for (i = 0; i < PK_SIZE; ++i)
        fscanf(fp, "%02hhx", &pk[i]);
    fclose(fp);

    // 1、2.C0 = encrypt(e,T), T-public key
    encrypt(c, pk, e);

    // 3.C1 = hash(2,e)
    hash(c1, content, 1 + GOPPA_N/8);

    // 4.key K = hash(1,e,C), C = (C0, C1), C0 = encrypt(e,T)
    content[0] = 1;
    hash(key, content, sizeof(content));

    
    if ((fp = fopen("ciphertext.txt","w")) == NULL)
    {
        printf("nerror in open ciphertext.txt.\n");
        exit(1);
    }

    for (i = 0; i < CIPHERTEXTBYTE; ++i)
        fprintf(fp, "%02hhx", c[i]);
    fclose(fp);

    // 5.output K
    printf("session key generated by encap:\n");
    for (i = 0; i < HASHBYTE; ++i)
        printf("%02hhx", key[i]);

    printf("\n\n\n");
    

    //return 0;
}