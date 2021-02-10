#include <openssl/md5.h>
#include <openssl/aes.h>
#include <string.h>

#include "randnum.h"

// input randsize: size of buf
// output buf: pseudo-random bits
void randomseq(unsigned char *buf, int randsize)
{
    unsigned char temp[16];
    int i, j = 0;
    size_t ctr_size = sizeof(ctr);
    unsigned char out[16];
    
    AES_KEY key;

    while (randsize > 0) 
    {
        for (i = 15; i >= 0; --i)
        {
            if (ctr[i] == 0xff)
                ctr[i] = 0;
            else 
            {
                ctr[i]++;
                break;
            }
        }

        MD5(ctr, ctr_size, temp);
        
        AES_set_encrypt_key(temp, 128, &key);
        AES_encrypt(ctr, out, &key);

        memcpy(buf+j, out, (randsize > 16 ? 16 : randsize));
        j += 16;
        randsize -= 16;

        // for (i = 0; i < 16; ++i)
        //     ctr[i] ^= temp[i];
    }

}
