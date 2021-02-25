# CM
&emsp;&emsp;每个文件夹是以参数设置的值命名的，按顺序依次是n，m，t，分别对应生成矩阵的列数量，有限域元素个数q的对数 \log_2{q} ，以及纠错能力t。参数末的f表示大小为t\*n的生成矩阵G从GF(q)映射到GF(2)的mt\*n的矩阵H在高斯消元后最后$\mu$行的主元列号不全等于行号，而且可变动范围在[0,$\nu$)。
&emsp;&emsp;编译运行在linux环境下，需要将Keccak Code Package和openssl编译生成静态库，然后执行make命令就能编译运行。
&emsp;&emsp;需要得到不同的结果，可以修改randnum.h的数组的值，每个数值不能超过unsigned char的范围[0,255]
