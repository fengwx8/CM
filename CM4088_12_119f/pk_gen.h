#ifndef PK_GEN_H
#define PK_GEN_H

#include <stdint.h>
#include "fq.h"

int pk_generation(unsigned char*, fq*, fq*, uint32_t *);

#endif