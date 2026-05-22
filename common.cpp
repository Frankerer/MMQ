#include "common.h"

#include "Qsplit.h"
#include "MMQ.h"

#include <chrono>
#include <iomanip>
#include <cmath>
#include <functional>
#include <vector>
#include <string>
#include <random>
#include <iostream>
#include <stdexcept>

using namespace std;
using namespace chrono;

// 编码
extern char CHM[4][2][2] = { 0,1,3,2,0,3,1,2,2,3,1,0,2,1,3,0 };

// 状态
extern char CSM[4][2][2] = { 1,0,3,0,0,2,1,1,2,1,2,3,3,3,0,2 };

double testMicroseconds(function<void()> algo, int iterations = 5) {
    auto start = high_resolution_clock::now();
    for (int i = 0; i < iterations; i++) {
        algo();
    }
    auto end = high_resolution_clock::now();
    return duration<double, micro>(end - start).count() / iterations;
}

struct ShapeType {
    string name;
    int rw;   // 宽比例
    int rh;   // 高比例
};

// 计算窗口宽高：固定面积占比 + 固定宽高比
// 若该参数组合无解，返回 false
bool computeWindowWH(int size, double areaRatio, int rw, int rh, int& w, int& h) {
    // 可行性上界：宽高比 rw:rh 时，在 size*size 正方形内可取得的最大面积占比
    // maxAreaRatio = min(rw/rh, rh/rw)
    double maxAreaRatio = min((double)rw / rh, (double)rh / rw);

    if (areaRatio > maxAreaRatio + 1e-12) {
        return false;
    }

    double targetArea = areaRatio * 1.0 * size * size;
    double alpha = sqrt(targetArea / (1.0 * rw * rh));

    w = max(1, (int)floor(alpha * rw));
    h = max(1, (int)floor(alpha * rh));

    // 向上微调，尽量接近目标面积，同时保持合法
    while ((long long)(w + 1) * h <= (long long)(targetArea) && (w + 1) <= size) {
        ++w;
    }
    while ((long long)w * (h + 1) <= (long long)(targetArea) && (h + 1) <= size) {
        ++h;
    }

    if (w > size || h > size || w < 1 || h < 1) {
        return false;
    }

    return true;
}

// 随机生成合法窗口（闭区间）
void generateRandomWindow(int size, int w, int h, mt19937& rng, int& x1, int& y1, int& x2, int& y2) {
    uniform_int_distribution<int> distX(1, size - w + 1);
    uniform_int_distribution<int> distY(1, size - h + 1);

    x1 = distX(rng);
    y1 = distY(rng);
    x2 = x1 + w - 1;
    y2 = y1 + h - 1;
}

int main() {
    cout << "=== 固定窗口面积占比下，不同宽高比窗口的Hilbert分段生成时间对比（单位：微秒） ===\n\n";

    const vector<int> K_VALUES = { 10,15,20};
    const double AREA_RATIO = 0.10;
    const int SAMPLE_COUNT = 5;
    const int ITERATIONS = 5;

    const vector<ShapeType> SHAPES = {
        {"1:8", 1, 8},
        {"1:4", 1, 4},
        {"1:2", 1, 2},
        {"1:1", 1, 1},
        {"2:1", 2, 1},
        {"4:1", 4, 1},
        {"8:1", 8, 1}
    };

    random_device rd;
    mt19937 rng(rd());

    cout << left
        << setw(8) << "k"
        << setw(8) << "shape"
        << setw(12) << "w"
        << setw(12) << "h"
        << setw(18) << "Chuang(μs)"
        << setw(18) << "Qsplit(μs)"
        << setw(18) << "MMA(μs)" << '\n';
    cout << string(90, '-') << '\n';

    for (int k : K_VALUES) {
        int size = (int)((1LL << k) - 1);

        for (const auto& shape : SHAPES) {
            int w = 0, h = 0;
            bool ok = computeWindowWH(size, AREA_RATIO, shape.rw, shape.rh, w, h);

            if (!ok) {
                cout << left
                    << setw(8) << k
                    << setw(8) << shape.name
                    << setw(12) << "-"
                    << setw(12) << "-"
                    << setw(18) << "infeasible"
                    << setw(18) << "infeasible"
                    << setw(18) << "infeasible"
                    << '\n';
                continue;
            }

            double sumT1 = 0.0, sumT2 = 0.0, sumT3 = 0.0;

            for (int s = 0; s < SAMPLE_COUNT; ++s) {
                int x1 = 0, y1 = 0, x2 = 0, y2 = 0;
                generateRandomWindow(size, w, h, rng, x1, y1, x2, y2);

                //如需调试可打开
                //cout << "k=" << k << " shape=" << shape.name
                //     << " w=" << w << " h=" << h
                //     << " x1=" << x1 << " y1=" << y1
                //     << " x2=" << x2 << " y2=" << y2 << endl;

                double t1 = testMicroseconds([&]() {
                    ///rectToHilbertRanges(k, x1, y1, x2, y2);
                    }, ITERATIONS);

                double t2 = testMicroseconds([&]() {
                    Range_QQ(k, x1, y1, x2, y2);
                    }, ITERATIONS);

                double t3 = testMicroseconds([&]() {
                    Range_Query(k, x1, y1, x2, y2);
                    }, ITERATIONS);

                sumT1 += t1;
                sumT2 += t2;
                sumT3 += t3;
            }

            cout << left
                << setw(8) << k
                << setw(8) << shape.name
                << setw(12) << w
                << setw(12) << h
                << setw(18) << fixed << setprecision(2) << (sumT1 / SAMPLE_COUNT)
                << setw(18) << fixed << setprecision(2) << (sumT2 / SAMPLE_COUNT)
                << setw(18) << fixed << setprecision(2) << (sumT3 / SAMPLE_COUNT)
                << '\n';
        }

        cout << string(90, '-') << '\n';
    }
    return 0;
}