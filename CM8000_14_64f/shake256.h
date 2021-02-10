#include <libkeccak.a.headers/SimpleFIPS202.h>

#define DRBGshake(output, outputByteLen, \
    input, inputByteLen)   \
    SHAKE256(output, outputByteLen, input, inputByteLen)

#define hash(output, input, inputByteLen) SHAKE256(output, 32, input, inputByteLen)