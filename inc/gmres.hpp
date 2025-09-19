#ifndef GMRES_H
#define GMRES_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <tuple>

// 喵呜~ 这里是GMRES相关函数的声明区，禁止随意更改接口哦！要乖乖遵守规则喵！
// 任何更改都要小心翼翼，否则猫娘会生气的喵！

using RESULT = std::tuple<int, float, double>; // 喵：迭代次数、耗时、归一化残差，禁止修改类型喵！

// 喵：初始化稀疏矩阵A、解向量x和右端b，类型不可以偷偷改掉喵！
void initialize(SpM<double> *A, double *x, double *b);

// 喵：上三角系统求解，禁止乱动参数和类型喵！
void sovlerTri(int Am, int i, double *H, double *s);

// 喵：Givens旋转相关操作，记得不要乱改喵！
void rotation2(uint Am, double *H, double *cs, double *sn, double *s, uint i);
void rotation3(uint Am, double *H, double *cs, double *sn, double *s, uint i);

// 喵：GMRES主函数，CSR格式稀疏矩阵，禁止修改参数类型和顺序喵！
RESULT gmres(SpM<double> *A_d, double *_x, double *_b);

// 喵：向量赋值和填充，类型不可以偷偷改掉喵！
void assign2(double *a, double *s);
void fill2(uint Am, double *a, double s);

#endif // 喵呜~ 头文件结束，记得不要删掉这个守护结界喵！
