#ifndef FQ_H
#define FQ_H

#include "params.h"


fq fq_add(fq, fq);
fq fq_mul(fq, fq);
fq fq_inv(fq);
fq fq_divi(fq, fq);
fq elementrev(fq);
void fqt_mul(fq *, fq *, fq *);

#endif