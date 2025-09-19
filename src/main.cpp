#include "sparseMatrix.hpp"
#include "gmres.hpp"
#include <fstream>
#include <iostream>
#include <vector>

// 喵喵～定义了一个结果类型哦！里面有整数、浮点数和双精度数喵！
using RESULT = std::tuple<int, float, double>;

// 喵喵～这个函数是用来加载二进制文件到CSR格式的稀疏矩阵喵～
void loadBinToCSR(SpM<double> &A, std::string mtx_path) {
  std::ifstream binFile(mtx_path, std::ios::in | std::ios::binary);

  // 喵喵～如果文件打不开就要告诉主人呢！
  if (!binFile.is_open()) {
    std::cerr << "喵喵～文件打不开呢！请检查路径喵！" << std::endl;
    return;
  }

  // 喵喵～从文件中读取矩阵的行数、列数和非零元素个数喵！
  binFile.read(reinterpret_cast<char *>(&A.nrows), sizeof(A.nrows));
  binFile.read(reinterpret_cast<char *>(&A.ncols), sizeof(A.ncols));
  binFile.read(reinterpret_cast<char *>(&A.nnz), sizeof(A.nnz));

  // 喵喵～给矩阵的行、列和值分配内存喵～
  A.rows = (uint *)malloc((A.nrows + 1) * sizeof(uint));
  A.cols = (uint *)malloc(A.nnz * sizeof(uint));
  A.vals = (double *)malloc(A.nnz * sizeof(double));

  // 喵喵～从文件中读取矩阵的行指针、列索引和非零值喵！
  binFile.read(reinterpret_cast<char *>(A.rows), (A.nrows + 1) * sizeof(uint));
  binFile.read(reinterpret_cast<char *>(A.cols), A.nnz * sizeof(uint));
  binFile.read(reinterpret_cast<char *>(A.vals), A.nnz * sizeof(double));
  binFile.close(); // 喵喵～记得关掉文件喵！
}

int main(int argc, char *argv[]) {
  // 喵喵～从命令行参数中获取矩阵名字喵～
  string mtx_name = argv[1];
  mtx_name.erase(mtx_name.begin(), mtx_name.begin() + mtx_name.rfind('/') +
                                       1); // 喵喵～去掉路径喵！
  mtx_name.erase(mtx_name.begin() + mtx_name.find('.'),
                 mtx_name.end()); // 喵喵～去掉后缀喵！

  // 喵喵～声明一个稀疏矩阵对象喵！禁止修改哦！
  SpM<double> A_double;            // 禁止修改
  loadBinToCSR(A_double, argv[1]); // 禁止修改

  // 喵喵～输出矩阵的基本信息喵～
  std::cout << mtx_name << ": M = " << A_double.nrows
            << ", N = " << A_double.ncols << std::endl;

  // 喵喵～初始化向量x和b喵！禁止修改哦！
  std::vector<double> x_double(A_double.nrows, 0); // 禁止修改
  std::vector<double> b_double(A_double.nrows, 1); // 禁止修改
  initialize(&A_double, &x_double[0], &b_double[0]);

  // 喵喵～开始运行GMRES算法喵！好期待结果呢喵！
  std::cout << "喵喵～开始运行GMRES算法啦！" << std::endl;
  auto res = gmres(&A_double, &x_double[0],
                   &b_double[0]); // GMRES算法核心，需重点优化喵！

  // 喵喵～输出算法的迭代次数、运行时间和残差喵～
  std::cout << "喵喵～迭代次数 = " << std::get<0>(res)
            << ", 耗时 = " << std::get<1>(res)
            << "ms, 残差 = " << std::get<2>(res) << std::endl;

  // ==============喵喵～禁止修改以下结果输出部分喵==============
  std::ofstream outfile("gmres_time.txt", std::ios::app);
  if (!outfile.is_open()) {
    std::cerr << "喵喵～结果文件打不开呢！请检查喵！" << std::endl;
    return 1;
  }
  outfile << mtx_name << " " << std::get<0>(res) << " " << std::get<1>(res)
          << " " << std::get<2>(res) << "\n";
  outfile.close(); // 喵喵～记得关掉文件喵！
  // ==============喵喵～禁止修改以上结果输出部分喵==============

  return 0; // 喵喵～程序结束啦！主人辛苦了喵！
}
