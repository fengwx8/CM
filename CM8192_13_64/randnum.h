#ifndef RANDNUM_H
#define RANDNUM_H

static unsigned char ctr[16] = 
    {0, 1, 2, 3, 4, 5, 6, 7,
     8, 9, 150, 11, 124, 13, 14, 15};

void randomseq(unsigned char*, int);

#endif