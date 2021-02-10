#include "fq.h"

// f(z) = z^15 + z + 1
// F(y) = y^119 + y^8 + 1

// input 2 elements in fq(fq1,fq2)
// return fq1 + fq2 (addition of finite field)
fq fq_add(fq fq1, fq fq2)
{
    return fq1 ^ fq2;
}

// input 2 elements in fq(fq1,fq2)
// return fq1 * fq2 (multiplication of finite field)
fq fq_mul(fq fq1, fq fq2)
{
	uint32_t mul = 0;
	uint32_t temp;
	int i;

	for (i = 0; i < M; ++i)
		mul ^= (fq1 * (fq2 & (1 << i)));

	temp = mul & 0x1FFF8000;
	mul ^= (temp >> 14) ^ (temp >> 15);

	return mul & (Q - 1);
}

// input a element fq1
// return 1/fq1 (the inversion of fq1 in finite field)
fq fq_inv(fq fq1)
{
	fq mul = 1;
	int i;
	for (i = 0; i < M-1; ++i)
	{
		mul = fq_mul(mul, fq1);
		mul = fq_mul(mul, mul);
	}

	return mul;
}

// input 2 elements in fq(fq1,fq2)
// return fq1 / fq2 (fraction of finite field)
fq fq_divi(fq fq1, fq fq2)
{   
    return fq_mul(fq1, fq_inv(fq2));
}

fq elementrev(fq element)
{
	element = ((element & 0x00FF) << 8) | ((element & 0xFF00) >> 8);
	element = ((element & 0x0F0F) << 4) | ((element & 0xF0F0) >> 4);
	element = ((element & 0x3333) << 2) | ((element & 0xCCCC) >> 2);
	element = ((element & 0x5555) << 1) | ((element & 0xAAAA) >> 1);
	
    // element = element >> (16-m)
	return element >> (16 - M);
}

void fqt_mul(fq *beta, fq *fqt1, fq *fqt2)
{
	int i,j;
	fq coeff[ 2*GOPPA_T - 1 ];

	for (i = 2*(GOPPA_T - 1); i >= 0; --i)
		coeff[i] = 0;


	// multiply
	for (i = 0; i < GOPPA_T; ++i)
		for (j = 0; j < GOPPA_T; ++j)
			coeff[ i+j ] ^= fq_mul(fqt1[i], fqt2[j]);

	// mod
	for (i = 2*(GOPPA_T - 1); i >= GOPPA_T; --i)
	{
		coeff[i - GOPPA_T + 8] ^= coeff[i];
		coeff[i - GOPPA_T ] ^= coeff[i];
	}

	for (i = 0; i < GOPPA_T; ++i)
		beta[i] = coeff[i];
}