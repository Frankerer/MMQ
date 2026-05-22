#pragma once
#include "common.h"

void mirror_copy_code(std::vector<bitmask_t>& vec, int iterStartPos, int lenth);
void all_border(std::vector<bitmask_t>& result_arr, int border, bitmask_t prevCode, int k, int nType, bitmask_t x_y_grid, int last_bit, int& oldPoint);
static void computeZsub(int Z, bitmask_t X, bitmask_t Y, bitmask_t half, int Zsub[2][2]);
static NR normalize(int T, int B, bitmask_t left, bitmask_t right, bitmask_t down, bitmask_t up, bitmask_t h);
static void split(int k, bitmask_t X, bitmask_t Y, int state, bitmask_t precode, int Z);
Sequence Range_Query(int k, bitmask_t qx1, bitmask_t qy1, bitmask_t qx2, bitmask_t qy2);