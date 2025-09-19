#include "sparseMatrix.hpp"
#include "gmres.hpp"
#include <assert.h>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <tuple>
#include <vector>

// 喵喵～使用标准命名空间喵～
using namespace std;

// 喵喵～重启次数固定为20哦！禁止修改喵～
const int RESTART_TIMES = 20; // 喵喵～禁止修改喵～
const double REL_RESID_LIMIT =
    1e-6; // 喵喵～这是相对残差的限制值哦～禁止修改喵～
const int ITERATION_LIMIT =
    10000; // 喵喵～最多的迭代次数就是这个啦～禁止修改喵～

void applyRotation(double &dx, double &dy, double &cs, double &sn) {
  // 喵喵～应用旋转操作啦！好厉害喵～
  double temp = cs * dx + sn * dy;
  dy = (-sn) * dx + cs * dy;
  dx = temp;
}

void generateRotation(double &dx, double &dy, double &cs, double &sn) {
  // 喵喵～生成旋转参数喵～
  if (dx == double(0)) {
    cs = double(0); // 喵喵～cos是0喵～
    sn = double(1); // 喵喵～sin是1喵～
  } else {
    double scale = fabs(dx) + fabs(dy); // 喵喵～计算尺度喵～
    double norm = scale * std::sqrt(fabs(dx / scale) * fabs(dx / scale) +
                                    fabs(dy / scale) *
                                        fabs(dy / scale)); // 喵喵～计算范数喵～
    double alpha = dx / fabs(dx); // 喵喵～这是符号喵～
    cs = fabs(dx) / norm;         // 喵喵～cos值喵～
    sn = alpha * dy / norm;       // 喵喵～sin值喵～
  }
}

void rotation2(uint Am, double *H, double *cs, double *sn, double *s, uint i) {
  // 喵喵～这是旋转操作的第二种方法喵～
  for (uint k = 0; k < i; k++) {
    applyRotation(H[k * Am + i], H[(k + 1) * Am + i], cs[k],
                  sn[k]); // 喵喵～应用旋转操作喵～
  }
  generateRotation(H[i * Am + i], H[(i + 1) * Am + i], cs[i],
                   sn[i]); // 喵喵～生成新的旋转参数喵～
  applyRotation(H[i * Am + i], H[(i + 1) * Am + i], cs[i],
                sn[i]); // 喵喵～应用旋转操作喵～
  applyRotation(s[i], s[i + 1], cs[i], sn[i]); // 喵喵～更新残差向量喵～
}

double calculateNorm(const double *vec, uint N) {
  // 喵喵～计算向量的范数喵～
  double sum = 0.0;
  for (uint i = 0; i < N; ++i) {
    sum += vec[i] * vec[i]; // 喵喵～累加平方值喵～
  }
  return std::sqrt(sum); // 喵喵～最后开平方喵～
}

void spmv(const uint *rowPtr, const uint *colInd, const double *values,
          const double *x, double *y, uint numRows) {
  // 喵喵～稀疏矩阵乘以向量喵～
  for (uint i = 0; i < numRows; ++i) {
    double sum = 0.0; // 喵喵～初始化每行的结果喵～
    for (uint j = rowPtr[i]; j < rowPtr[i + 1]; ++j) {
      sum += values[j] * x[colInd[j]]; // 喵喵～累加乘积喵～
    }
    y[i] = sum; // 喵喵～保存结果喵～
  }
}

double dotProduct(const double *x, const double *y, uint N) {
  // 喵喵～计算两个向量的点积喵～
  double sum = 0.0;
  for (uint i = 0; i < N; ++i) {
    sum += x[i] * y[i]; // 喵喵～累加乘积喵～
  }
  return sum; // 喵喵～返回结果喵～
}

void daxpy(double alpha, const double *x, double *y, uint N) {
  // 喵喵～计算 y += alpha * x 喵～
  for (uint i = 0; i < N; ++i) {
    y[i] += alpha * x[i]; // 喵喵～更新y喵～
  }
}

void dscal(double alpha, double *x, uint N) {
  // 喵喵～计算 x *= alpha 喵～
  for (uint i = 0; i < N; ++i) {
    x[i] *= alpha; // 喵喵～更新x喵～
  }
}

void dcopy(const double *src, double *dst, uint N) {
  // 喵喵～复制向量喵～
  for (uint i = 0; i < N; ++i) {
    dst[i] = src[i]; // 喵喵～把src拷贝到dst喵～
  }
}

void sovlerTri(int Am, int i, double *H, double *s) {
  // 喵喵～求解上三角系统喵～
  for (int j = i; j >= 0; j--) {
    s[j] /= H[Am * j + j]; // 喵喵～更新s喵～
    for (int k = j - 1; k >= 0; k--) {
      s[k] -= H[k * Am + j] * s[j]; // 喵喵～更新上三角系统喵～
    }
  }
}

RESULT gmres(SpM<double> *A_d, double *x_d, double *_b) {
  // 喵喵～GMRES算法的核心喵～
  const uint N = A_d->nrows;

  // 喵喵～定义一些必要的向量和变量喵～
  std::vector<double> r0(N);
  std::vector<double> V((RESTART_TIMES + 1) * N);
  std::vector<double> s(RESTART_TIMES + 1, 0.0);
  std::vector<double> V0(N);
  std::vector<double> H((RESTART_TIMES + 1) * RESTART_TIMES);
  std::vector<double> cs(RESTART_TIMES);
  std::vector<double> sn(RESTART_TIMES);

  double beta = calculateNorm(_b, N); // 喵喵～计算初始残差喵～
  double RESID_LIMIT = REL_RESID_LIMIT * beta; // 喵喵～计算残差限制喵～
  double init_res = beta; // 喵喵～保存初始残差喵～

  int iteration = 0;

  auto start = std::chrono::high_resolution_clock::now(); // 喵喵～开始计时喵～

  do {
    // 喵喵～外迭代开始喵～
    spmv(A_d->rows, A_d->cols, A_d->vals, x_d, r0.data(),
         N);                            // 喵喵～计算 r0 = A*x 喵～
    daxpy(-1.0, _b, r0.data(), N);      // 喵喵～计算残差 r0 喵～
    beta = calculateNorm(r0.data(), N); // 喵喵～更新残差范数喵～

    dscal(-1.0 / beta, r0.data(), N); // 喵喵～归一化残差向量喵～
    dcopy(r0.data(), V.data(), N);    // 喵喵～复制到V喵～

    fill(s.begin(), s.end(), 0.0); // 喵喵～初始化残差向量喵～
    s[0] = beta;                   // 喵喵～设置初始值喵～

    int i = -1;
    do {
      // 喵喵～内迭代开始喵～
      i++;
      iteration++;
      std::vector<double> V_i(N);
      dcopy(V.data() + i * N, V_i.data(), N); // 喵喵～复制当前V_i喵～

      spmv(A_d->rows, A_d->cols, A_d->vals, V_i.data(), r0.data(),
           N); // 喵喵～计算 A*V_i 喵～

      for (int k = 0; k <= i; k++) {
        H[k * RESTART_TIMES + i] =
            dotProduct(r0.data(), V.data() + k * N, N); // 喵喵～计算H喵～
        daxpy(-H[k * RESTART_TIMES + i], V.data() + N * k, r0.data(),
              N); // 喵喵～更新r0喵～
      }
      H[(i + 1) * RESTART_TIMES + i] =
          calculateNorm(r0.data(), N); // 喵喵～计算H喵～
      dscal(1.0 / H[(i + 1) * RESTART_TIMES + i], r0.data(),
            N);                                    // 喵喵～归一化喵～
      dcopy(r0.data(), V.data() + N * (i + 1), N); // 喵喵～更新V喵～

      rotation2(RESTART_TIMES, H.data(), cs.data(), sn.data(), s.data(),
                i); // 喵喵～应用旋转喵～

      if (std::abs(s[i + 1]) <= RESID_LIMIT || iteration >= ITERATION_LIMIT) {
        break; // 喵喵～结束内迭代喵～
      }
    } while (i + 1 < RESTART_TIMES && iteration <= ITERATION_LIMIT);

    sovlerTri(RESTART_TIMES, i, H.data(), s.data()); // 喵喵～求解上三角系统喵～

    for (int j = 0; j <= i; j++) {
      daxpy(s[j], V.data() + j * N, x_d, N); // 喵喵～更新解喵～
    }
  } while (std::abs(s[RESTART_TIMES]) > RESID_LIMIT &&
           iteration <= ITERATION_LIMIT);

  auto stop = std::chrono::high_resolution_clock::now(); // 喵喵～结束计时喵～
  std::chrono::duration<float, std::milli> duration =
      stop - start;                   // 喵喵～计算耗时喵～
  float test_time = duration.count(); // 喵喵～保存耗时喵～

  return make_tuple(iteration, test_time,
                    std::abs(s[RESTART_TIMES]) /
                        init_res); // 喵喵～返回结果喵～
}

void initialize(SpM<double> *A, double *x, double *b) {
  // 喵喵～初始化向量喵～
  int N = A->nrows;

  for (int i = 0; i < N; i++) {
    x[i] = sin(i); // 喵喵～用sin函数初始化x喵～
  }

  double beta = calculateNorm(x, N); // 喵喵～计算范数喵～
  for (uint i = 0; i < N; i++) {
    x[i] /= beta; // 喵喵～归一化x喵～
  }

  spmv(A->rows, A->cols, A->vals, x, b, N); // 喵喵～计算b喵～

  for (uint i = 0; i < N; i++) {
    x[i] = 0.0; // 喵喵～初始化x为0喵～
  }
}
