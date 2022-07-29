# <center>Program Homework 1</center>

## 1 优化问题描述

优化如下问题

![优化问题描述](imgs/优化问题描述.png)

特别地，我们取 T = 1，则问题转化为

![简化问题](imgs/简化问题.png)

## 2 算法原理

采用 ADMM 算法，算法原理如下：

![算法原理](imgs/算法原理.png)

其中，针对三个子问题，由于解空间为n维复数空间全体，所以我们令函数导数等于0，从而得到子问题的解

## 3 程序使用指南

### 3.1 IDE

Visual Studio 2019

### 3.2 平台工具集

Visual Studio 2019 (v142)

### 3.3 编程语言

C++

### 3.4 依赖库

Eigen

### 3.5 函数调用

最简单的调用方式为

![函数调用](imgs/函数调用.png)

程序会按照所给定的维数自动初始化所有需要的变量，包括 x ；并将最终的计算结果赋值给 x 本身。也可按如下方式对所有变量进行初始化

![人工初始化调用](imgs/人工初始化调用.png)

## 4 程序测试

### 4.1 二维测试情形

A = ![A](imgs/A.png)

b = ![b](imgs/b.png)

计算结果为

![result](imgs/result.png)

### 4.2 高维测试情形

A = ![A_high](imgs/A_high.png)

b = ![b_high](imgs/b_high.png)

计算结果为

![result_high](imgs/result_high.png)

### 4.2 更高维（1000维）测试情形

from

![norm_1000x50](imgs/norm_1000x50.png)

to

![norm_result_1000_1](imgs/norm_result_1000_1.png)

![norm_result_1000_2](imgs/norm_result_1000_2.png)

## 5 结论

该稀疏优化算法针对大规模矩阵的求解有着普遍高效的收敛效率，能够在几次迭代之后迅速得到最优解，是一个高效的优化算

## 6 其他测试

用户可随意选择初始矩阵和参数，进行进一步测试

## Enjoy it~

