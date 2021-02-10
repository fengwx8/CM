#ifndef LITTLE_ENDIAN_H
#define LITTLE_ENDIAN_H

#include <stdint.h>

uint16_t load2(unsigned char*);

void store2(unsigned char*, uint16_t);

uint32_t load4(unsigned char *);

void store4(unsigned char*, uint32_t);

unsigned char byterev(unsigned char);

#endif