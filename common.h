#pragma once
#include <stdio.h>
#include <iostream>
#include <cmath> 
#include <vector>
#include <bitset>
#include <chrono>
#include <string>
#include <algorithm>
#include <queue>
#include <random>
#include <fstream>
#include <iomanip>
#include <numeric>

using namespace std;
typedef long long bitmask_t;
typedef long halfmask_t;
typedef std::vector<bitmask_t> Row;
typedef vector<vector<bitmask_t>> TwoDArray;

typedef std::pair<bitmask_t, bitmask_t> Segment;
// 定义一个序列
typedef std::vector<Segment> Sequence;


#define getOneBitByPos(X,bit) ((X>>bit) & 0x01)  
#define last_k_Zero(n,k) (n & ((~0U) << k))
#define getLowKBits(X,k) (X & ((1 << k) - 1))
#define lsb32_idx(v) (N_bit[((uint32_t)((v & -v) * 0x077CB531U)) >> 27])


static const int N_bit[32] = { 0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8, 31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9 };
static const int M_bit[32] = { 0,  9,  1, 10, 13, 21,  2, 29, 11, 14, 16, 18, 22, 25,  3, 30, 8, 12, 20, 28, 15, 17, 24,  7, 19, 27, 23,  6, 26,  5,  4, 31 };


//编码
//the mappings from L1-coord to its L1-code  从L1坐标到其L1代码的映射
extern char CHM[4][2][2];

//the mappings from L1-coord to its Ln-state 从L1坐标到其Ln状态的映射
extern char CSM[4][2][2];


struct Block { int x, y, k; };
struct HilbertRange { bitmask_t min_code, max_code; };
struct NR { int Bp; bitmask_t left, right, down, up; };

