/*
    This file refers to the following paper
    Verified fast formulas for control bits for permutation networks
*/

#include <string.h>
#include <stdlib.h>
#include "controlbits.h"
#include "permutation.h"

#define MINMAX32(a,b) \
do{ \
    uint32_t c = b - a; \
    c >>= 31; \
    c = -c; \
    c &= a ^ b; \
    a ^= c; \
    b ^= c; \
}while(0)

/* e consists of a and b, that is, e = (a<<16)|b */
/* after sorting e = (id << 16) | (b·(a^-1)) */
/* or consists of a, b and c, e = (a << 20) | (b << 10) | c*/
/* after sorting e = (id << 20) | ((b·(a^-1)) << 10) |  (c·(a^-1))*/
void sort32(uint32_t *e,int num)
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
                MINMAX32(e[j], e[i+j]);
        j = 0;

        for (m = bound; m > i; m >>= 1)
        {
            for (; j < num - m; ++j)
            {
                if(!(j & i)){
                    uint32_t temp = e[i + j];
                    for (n = m; n > i; n >>= 1)
                        MINMAX32(temp, e[j+n]);
                    e[i + j] = temp;
                }
            }
        }
    }
}

/* output cbits: control bits                       */
/* input pos: pointer                               */
/*       gap: layer, power of 2 from 2^0 to 2^(m-1) */
/*       pi: permutation                            */
/*       m: parameter m                             */
/*       s: s = (1 << m)                            */
/*       buf: store intermediate result             */
void recursion(unsigned char *cbits, int pos, int gap, 
    uint16_t *pi, int m, int s, uint32_t *buf){
#define _for(x,a,b) for(x = a; x < b; ++x)   

    if (m == 1)
    {
        cbits[pos >> 3] ^= pi[0] << (pos & 7);
        return;
    }
    uint32_t *buf1 = buf, *buf2 = buf + s;
    uint16_t *pi_h = (uint16_t *)buf2;
    int x,i,j,k;

     _for(x,0,s)
        buf1[x] = ((pi[ x ] ^ 1) << 16)|pi[ x^1 ];
    sort32(buf1,s); /* buf1 = id << 16 | pibar, p = pibar */

    _for(x,0,s)
    {
        uint32_t pibar = buf1[x] & 0xffff;
        buf2[x] = pibar << 16;
        buf2[x] |= (x < pibar ? x : pibar);
    }
    /* buf2 = pibar << 16 | c, c = min(id,p)*/

    _for(x,0,s) 
        buf1[x] = (buf1[x] << 16) | x;
    sort32(buf1,s); /* buf1 = (id << 16) | 1/pibar */
    /* 1/pibar = q, 1/pibar is the representation of pibar^(-1) */

    _for(x,0,s)
        buf1[x] = (buf1[x] << 16) | (buf2[x] >> 16);
    /* buf1 = (q << 16) | p, or (pibar^(-1) << 16) | pibar */

    sort32(buf1,s); /* buf1 = (id << 16) | pibar^2 ,p = pibar^2*/

    if (m <= 10)
    {
        _for(x,0,s)
            buf2[x] = ((buf1[x]&0xffff) << 10) | (buf2[x]&0x3ff);
        /* buf2 = (p << 10) | c, p =pibar^2, q = p^(-1) */

        _for(i,1,m-1)
        {
            /* buf2 = (p << 10) | c */
            _for(x,0,s)
                buf1[x] = ((buf2[x]&0xfffffc00) << 6) | x;
            sort32(buf1,s); /* buf1 = (id << 16) | q */

            _for(x,0,s)
                buf1[x] = (buf1[x] << 20) | buf2[x];
            /* bubf1 = (q << 20) | (p << 10) | c */
            sort32(buf1,s);
            /* buf1 = (id << 20) | (p·p << 10) | c·p */

            _for(x,0,s){
                uint32_t c = buf2[x] & 0x3ff;
                uint32_t cp = buf1[x] & 0x3ff;
                buf2[x] = buf1[x] & 0xffc00;
                buf2[x] |= (c < cp ? c : cp);
            }
        }

        _for(x,0,s)   
            buf2[x] &= 0x3ff;     
    }
    else
    {
        _for(x,0,s)
            buf2[x] = (buf1[x] << 16) | (buf2[x] & 0xffff);

         _for(i,1,m-1){
            /* buf2 = (p << 16) | c */
             _for(x,0,s)
                buf1[x] = (buf2[x] & 0xffff0000) | x;
            sort32(buf1,s);
            /* buf1 = (id << 16) | q */

            _for(x,0,s)
                buf1[x] = (buf1[x] << 16) | (buf2[x] & 0xffff);
            /* buf1 = (q << 16) | c */

            if (i < m-2){
                _for(x,0,s)
                    buf2[x] = (buf1[x] & 0xffff0000) | (buf2[x] >> 16);
                /* buf2[x] = (q << 16) | p */
                sort32(buf2,s);
                /* buf2 = (id << 16) | p^2 */
                _for(x,0,s)
                    buf2[x] = (buf2[x] << 16) | (buf1[x] & 0xffff);
                /* buf2[x] = (p^2 << 16) | c */
            }

            sort32(buf1,s); /* buf1[x] = (id << 16) | c·p */
            _for(x,0,s){
                uint32_t cp = buf1[x] & 0xffff;
                uint32_t c = buf2[x] & 0xffff;
                buf2[x] ^= (cp < c ? (cp ^ c) : 0);
            }
        }
        _for(x,0,s)
            buf2[x] &= 0xffff;
    }
    /* buf2 = cycmin(pibar) */
    
    _for(x,0,s) buf1[x] = (((uint32_t)pi[x]) << 16) | (uint32_t)x;
    sort32(buf1,s); /* buf1 = (id << 16) | pi^(-1) */

    _for(j,0,s/2){
        int index = 2*j;
        uint32_t fj = buf2[index] & 1;// f[j] = c[2*j]%2
        uint32_t Fx = fj ^ (uint32_t)index, Fx1 = Fx ^ 1;

        cbits[pos >> 3] ^= fj << (pos & 7);
        pos += gap;

        buf2[index] = (buf1[index] << 16) | Fx;
        buf2[index+1] = (buf1[index+1] << 16) | Fx1;
    }

    sort32(buf2,s); // buf2 = (id << 16) | F·pi

    pos += (2*m-3)*gap*(s/2);

    _for(k,0,s/2){
        int y = 2*k;
        uint32_t lk = buf2[y] & 1; // l[j] = Fpi[2*j]%2
        uint32_t Ly = lk ^ (uint32_t)y, Ly1 = Ly ^ 1;

        cbits[pos >> 3] ^= lk << (pos & 7);
        pos += gap;

        buf1[y] = (Ly << 16) | (buf2[y] & 0xffff);
        buf1[y+1] = (Ly1 << 16) | (buf2[y+1] & 0xffff);
    }

    sort32(buf1,s); // buf1 = (id << 16) | M
    pos -= (2*m-2)*gap*(s/2);

    _for(j,0,s/2){
        pi_h[j] = (buf1[2*j]&0xffff) >> 1;
        pi_h[j+s/2] = (buf1[2*j+1]&0xffff) >> 1;
    }

    recursion(cbits,pos,2*gap,pi_h,m-1,s/2,buf);
    recursion(cbits,pos+gap,2*gap,pi_h+s/2,m-1,s/2,buf);
}

/* output: cbits, store controlbits generated by specific pi */
/* input: (m,s): 2^m = s, s is the size of pi, 1 <= m <= 15*/
/*         pi, a permutation from 0 to s-1*/
void controlbits(unsigned char *cbits, uint16_t *pi, int m, int s)
{
    int i,equal;
    uint32_t *buf = (uint32_t *)malloc(2 * s * sizeof(uint32_t));
    uint16_t *pi_t = (uint16_t *)malloc(s * sizeof(uint16_t));

    while (1)
    {
        memset(cbits, 0, (((2*m-1)*s/2)+7)/8);
        recursion(cbits,0,1,pi,m,s,buf);

        for (i = 0; i < s; ++i)
            pi_t[i] = i;
        permutation(pi_t,cbits,m);

        equal = 1;
        for (i = 0; i < s; ++i)
            if(pi_t[i] != pi[i]){
                equal = 0;
                break;
            }

        if (equal)
            break;            
    }

    free(buf);
    free(pi_t);
}