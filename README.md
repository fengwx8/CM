# CM
&emsp;&emsp;每个文件夹是以参数设置的值命名的，按顺序依次是n，m，t，分别对应生成矩阵的列数量，有限域元素个数q的对数$\log_2 q$，以及$\vec{a}$纠错能力t。参数末的f表示半系统形矩阵。  
&emsp;&emsp;编译运行在linux环境下，需要将Keccak Code Package和openssl编译生成静态库，然后执行make命令就能编译运行。  
&emsp;&emsp;需要得到不同的结果，可以修改randnum.h的数组的值，每个数值不能超过unsigned char的范围[0,255]
