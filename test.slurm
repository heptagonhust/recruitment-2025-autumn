#!/bin/bash
#SBATCH -o %j.out			    #指定作业的标准输出
#SBATCH -e %j.err			    #指定作业的错误输出
#SBATCH -J gmres			    #指定作业名称
#SBATCH -t 02:00:00			    #指定作业最大运行30分钟
#SBATCH -N 1				    #指定节点数

module purge
module load compiler/gcc/12.2.0
module load compiler/dtk/25.04

make clean && make

if [ ! -f ./gmres ]; then
    echo "Compilation failed, executable ./gmres not found."
    exit 1
fi

# 矩阵文件所在的基础目录
MATRIX_DIR="/data/Recruitment/GMRES/sparseMatrixSet"

for i in $(seq 1 80)
do
    MTX_FILE="${MATRIX_DIR}/mtx${i}.bin"
    
    if [ -f "$MTX_FILE" ]; then
        echo "================================================="
        echo "Running GMRES for mtx${i}.bin"
        echo "================================================="
        ./gmres "$MTX_FILE"
    else
        echo "Warning: Matrix file $MTX_FILE not found, skipping."
    fi
done

echo "All tests finished."
