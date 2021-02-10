#include "poly_gen.h"


/* input element: a element in fq((2^m)^t)*/
/* output poly: irreducible-polynomial */
/* return false for success and true for failure */
bool poly_generation(fq *poly, fq *element)
{
	fq matrix[ GOPPA_T+1 ][ GOPPA_T ];

	int i,j,m;
	fq inverse,coeff;

	// calculate a column vector [1 β β^2 ... β^t]
	matrix[0][0] = 1;
	for (i=1; i < GOPPA_T; ++i)
		matrix[0][i] = 0;

	for (i = 0; i < GOPPA_T; ++i)
	matrix[1][i] = element[i] & (Q - 1);

	for (i = 2; i <= GOPPA_T; ++i)
	fqt_mul(matrix[i], matrix[i-1], element);

	
	for (i = 0; i < GOPPA_T; ++i)
	{
		// ensure the pivot element nonzero
		// method: if it equals 0, add right column
		// elementary column transformation will not change the solution
		j = i + 1;
			
		// column[i] = column[i] + column[i+1]
		// if matrix[i][i] still equals 0, add column[i+2]
		// and so on, until i = t or matrix[i][i] is nonzero
		while((matrix[i][i] == 0) && (j < GOPPA_T))
		{
			for (m = i; m <= GOPPA_T; ++m)
			matrix[m][i] ^= matrix[m][j];
			++j;
		}

		// if the matrix is not GOPPAtematic, then return true
		if (matrix[i][i] == 0)
			return true;

		inverse = fq_inv(matrix[i][i]);

		// for (m = 0; m <= GOPPA_T; ++m){
		// 	for (j = 0; j < GOPPA_T; ++j)
		// 		printf("%u ",matrix[m][j]);
		// 	printf("\n");
		// }

		for (j = i; j <= GOPPA_T; ++j)
			matrix[j][i] = fq_mul(matrix[j][i],inverse);
		
		// Gauss elimination
		for (j = 0; j < GOPPA_T; ++j)
		{
			if (i == j) continue;

			coeff = matrix[i][j];

			for(m = i; m <= GOPPA_T; ++m)
				matrix[m][j] ^= fq_mul(coeff,matrix[m][i]);
		}
	}

	for(i = 0; i < GOPPA_T; ++i)
		poly[i] = matrix[GOPPA_T][i];

	return false;
}