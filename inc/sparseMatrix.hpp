#ifndef SPARSE_MATRIX
#define SPARSE_MATRIX

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cstring>

// 喵呜~ 这里是稀疏矩阵相关定义，禁止随意更改接口和类型喵！
// 任何更改都要小心翼翼，否则猫娘会生气的喵！

typedef unsigned int uint; // 喵：无符号整型，矩阵索引专用，禁止修改类型喵！
using namespace std;

// 喵：稀疏矩阵模板类，默认CSR格式存储，记得不要乱动成员变量和构造析构喵！
template <typename ValueType>
class SpM
{
public:
    uint nrows;      // 喵：行数
    uint ncols;      // 喵：列数
    uint nnz;        // 喵：非零元个数
    uint *rows;      // 喵：CSR行指针
    uint *cols;      // 喵：CSR列索引
    ValueType *vals; // 喵：非零元数值

    // 喵：默认构造，禁止乱动喵！
    SpM() : nrows(0), ncols(0), nnz(0), rows(NULL), cols(NULL), vals(NULL) {};

    // 喵：指定大小构造，记得不要改动分配方式和精度喵！
    SpM(uint nrows, uint ncols, uint nnz) : nrows(nrows), ncols(ncols), nnz(nnz) {
        rows = (uint *)malloc((nrows + 1) * sizeof(uint));
        memset(rows, 0, (nrows + 1) * sizeof(uint));
        cols = (uint *)malloc((nnz) * sizeof(uint));
        memset(cols, 0, (nnz) * sizeof(uint));
        vals = (ValueType *)malloc((nnz) * sizeof(ValueType));
        memset(vals, 0, (nnz) * sizeof(ValueType));
    };

    // 喵：拷贝构造，禁止更改拷贝逻辑喵！
    SpM(const SpM<ValueType> &A) {
        nrows = A.nrows;
        ncols = A.ncols;
        nnz = A.nnz;
        rows = (uint *)malloc((nrows + 1) * sizeof(uint));
        cols = (uint *)malloc((nnz) * sizeof(uint));
        vals = (ValueType *)malloc((nnz) * sizeof(ValueType));
        for (unsigned int i = 0; i < nrows + 1; i++)
        {
            rows[i] = A.rows[i];
        }
        for(unsigned int p = 0; p < nnz; p++) {
            vals[p] = A.vals[p];
            cols[p] = A.cols[p];
        }
    }

    // 喵：析构函数，记得释放内存，不要乱动喵！
    ~SpM() {
        if(rows) free(rows);
        if(cols) free(cols);
        if(vals) free(vals);
    }
};

#endif // 喵呜~ 头文件结束，记得不要删掉这个守护结界喵！
