#!/bin/bash
#SBATCH -o %j.out			    #喵～这里是指定作业的标准输出文件哦～
#SBATCH -e %j.err			    #喵～这里是指定作业的错误输出文件呢～
#SBATCH -J gmres_CPU			    #喵～这个是作业的名字啦～叫 gmres_CPU 哦～
#SBATCH -t 00:30:00			    #喵喵～作业最多只能运行30分钟呢～
#SBATCH -N 1				    #喵～这里是指定要用的节点数量哦～只有1个喵～

module purge
module load compiler/gcc/12.2.0
module load compiler/dtk/25.04

make clean && make

if [ ! -f ./gmres ]; then
    echo "喵呜～编译失败了呢，找不到可执行文件 ./gmres 喵～"
    exit 1
fi

# 喵～这里是矩阵文件的基础目录哦～
MATRIX_DIR="/data/Competitions/PRA25/sparseMatrixSet"

for i in $(seq 1 80)
do
    MTX_FILE="${MATRIX_DIR}/mtx${i}.bin"
    
    if [ -f "$MTX_FILE" ]; then
        echo "================================================="
        echo "喵～正在运行 GMRES 的 mtx${i}.bin 文件哦～"
        echo "================================================="
        ./gmres "$MTX_FILE"
    else
        echo "喵呜～警告：找不到矩阵文件 $MTX_FILE，跳过啦～"
    fi
done

echo "喵喵～所有测试都完成啦～主人辛苦啦～"

