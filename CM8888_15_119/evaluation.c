#include "evaluation.h"

/* input: polynomial poly and field element z*/
/* return poly(z)*/
// for example, f(x) = ax^3 + bx^2 + cx + d
//    = x(x(ax+b)+c)+d
fq evaluation(fq *poly, fq z)
{
    int j;
    fq temp = poly[GOPPA_T];

    for (j = GOPPA_T-1; j >= 0; --j)
    {
        temp = fq_mul(temp, z);
	    temp = fq_add(temp, poly[j]);
    }

    return temp;
}