#include <stdio.h>
#include <stdint.h>

//#define M 4
//
//#define Q (1 << M) 

typedef uint16_t fq;

int isdivide(fq a,int m, fq b, int d)
{
	if (d == 0) return 0;
	fq c = b << (m-d);
	int i;
	
	for (i = m; i >= d; --i){
		if ((a >> i) & 1)
			a ^= c;
		c >>= 1;
	} 
	
	if (a == 0)
		return 1;
	else
		return 0;
}

int recursion(fq a, int m, fq b, int d){
	if (d == m) return 0;
	fq c = b << 1;
	
//	if (recursion(a, m, c, d+1))
//		return 1;
//	if (recursion(a, m, c+1, d+1))
//		return 1;
	if (recursion(a, m, c, d+1) || recursion(a, m, c+1, d+1))
		return 1;
	if (isdivide(a, m, b, d)) 
		return 1;
	return 0;
}

int main(){
	
	// 2 <= M <= 15
	int M, Q;
	printf("input M(within 1 to 15):");
	scanf("%d",&M);
	Q = 1 << M;
	
	int i = Q;
	for (; i < 2*Q; ++i)
	{
		if(!recursion((fq)i, M, (fq)1, 0))
			printf("%x is a irreducible polynomial in hexadecimal.\n",i);
	}
	
	printf("\n\nMost significant bit means the coefficient of the highest degree of polynomial.\n" \
			"Least significant bit means the constant.\n");
	return 0;
}
