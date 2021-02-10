#include <stdio.h>

#include "keygen.h"
#include "encap.h"
#include "decap.h"

//#define SAMPLE_NUM 5

int main()
{
    int i;
    for (i = 0; i < SAMPLE_NUM; ++i)
    {
        printf("sample %d:\n",i);
        keygen();
        encap();
        decap();
    }

    return 0;
}
