#include "little_endian.h"

uint16_t load2(unsigned char* s)
{
    uint16_t t = 0;
    t = s[1] << 8;
    return (t | s[0]);
}

void store2(unsigned char* s, uint16_t t)
{
    s[0] = t &0xff;
    s[1] = t >> 8;
}

uint32_t load4(unsigned char *s)
{
    uint32_t t = 0;
    t = s[3];
    t <<= 8;
    t |= s[2];
    t <<= 8;
    t |= s[1];
    t <<= 8;
    return (t | s[0]);
}
void store4(unsigned char *s, uint32_t t)
{
    int i = 0;
    for (;i < 4; ++i) {
        s[i] = t & 0xff;
        t >>= 8;
    }
}
unsigned char byterev(unsigned char byte)
{
    byte = ((byte & 0x0f) << 4) | ((byte & 0xf0) >> 4);
    byte = ((byte & 0x33) << 2) | ((byte & 0xcc) >> 2);
    byte = ((byte & 0x55) << 1) | ((byte & 0xaa) >> 1);

    return byte;
}