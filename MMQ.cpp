#include "common.h"

static bitmask_t G_pos[4];  // 边框坐标
static int       G_k[4];    // 触发阶数
static bitmask_t QX1, QX2, QY1, QY2;
static vector<bitmask_t> EntryExit;


// ORDER[state][code] = (x<<1)|y
// 按 Hilbert 序 code=0..3 给出对应象限的几何坐标 (x,y)
//0,1,3,2,0,2,3,1，3,2,0,1,3,1,0,2
static const int ORDER[4][4] = {
    {0b00, 0b01, 0b11, 0b10},  // state=0
    {0b00, 0b10, 0b11, 0b01},  // state=1
    {0b11, 0b10, 0b00, 0b01},  // state=2
    {0b11, 0b01, 0b00, 0b10},  // state=3
};


//镜像复制函数——当前数组，该层阶数，该阶状态
void mirror_copy_code(std::vector<bitmask_t>& vec, int iterStartPos, int lenth) {
    int last_bit_arr = ((1 << ((iterStartPos + 2) << 1)) - 1);
    int originalSize = vec.size();
    vec.resize(vec.size() + lenth);
    for (int i = 0; i < lenth; ++i) {
        vec[originalSize + i] = vec[originalSize - 1 - i] ^ last_bit_arr;
    }
}

//数组、边框、前缀、位标识（阶数-1）、状态、||、所在边（阶内边位）、所在边阶数（(10)2偶数位为2）、记录数（归一化到同状态）
void all_border(std::vector<bitmask_t>& result_arr, int border, bitmask_t prevCode, int k, int nType, bitmask_t x_y_grid, int last_bit, int& oldPoint) {
    //归一化操作到状态为0的状态，以减少分类
    border ^= nType;
    if (nType >= 2)  x_y_grid = ((1 << (k + 1)) - x_y_grid);
    // 底层处理 
    if (k < last_bit) {
        prevCode <<= ((k + 1) << 1);
        int mask = 1 << (k << 1);
        switch (border) {
        case 0:
            result_arr.push_back(prevCode + mask - 1);
            result_arr.push_back(prevCode + 3 * mask);
            oldPoint += 2;
            break;
        case 1:
            result_arr.push_back(prevCode + 2 * mask - 1);
            oldPoint += 1;
            break;
        case 2:
            result_arr.push_back(prevCode + mask);
            result_arr.push_back(prevCode + 3 * mask - 1);
            oldPoint += 2;
            break;
        case 3:
            result_arr.push_back(prevCode + 2 * mask);
            oldPoint += 1;
            break;
        }
        return;
    }
    int newPoint = 0, temp = oldPoint;
    bool bitX_Y = getOneBitByPos(x_y_grid, k);
    //按Hilbert序分化两个子区域
    switch (border) {
    case 0: case 2:
        //启用镜像复制
        all_border(result_arr, border, (prevCode << 2) | CHM[0][0][bitX_Y], k - 1, CSM[0][0][bitX_Y], x_y_grid, last_bit, newPoint);
        mirror_copy_code(result_arr, k - 1, newPoint);
        oldPoint = temp + (newPoint << 1);
        //不启用镜像复制
        /*all_border(result_arr, border, (prevCode << 2) | CHM[0][0][bitX_Y], k - 1, CSM[0][0][bitX_Y], x_y_grid, last_bit, oldPoint);
        all_border(result_arr, border, (prevCode << 2) | CHM[0][1][bitX_Y], k - 1, CSM[0][1][bitX_Y], x_y_grid, last_bit, oldPoint);*/
        break;
    case 1: case 3:
        all_border(result_arr, border, (prevCode << 2) | (2 * bitX_Y), k - 1, !bitX_Y, x_y_grid, last_bit, oldPoint);
        all_border(result_arr, border, (prevCode << 2) | (2 * bitX_Y + 1), k - 1, bitX_Y ? 3 : 0, x_y_grid, last_bit, oldPoint);
        break;
    }
}

// ============================================================
//  精确的 Zsub 分化
//
//  对每条活跃边框 i：
//    取 G_pos[i] 的第 (k-1) 位 bit：边框在哪一行/列
//    判断查询矩形在另一维度覆盖哪些子区域
//      水平边框(i=0,2)：判断 X 方向 → QX1<X+half(进入x=0) 且/或 QX2>X+half(进入x=1)
//      垂直边框(i=1,3)：判断 Y 方向 → QY1<Y+half(进入y=0) 且/或 QY2>Y+half(进入y=1)
//
//  结果写入 Zsub[x][y]
// ============================================================
static void computeZsub(int Z, bitmask_t X, bitmask_t Y, bitmask_t half, int Zsub[2][2]) {

    Zsub[0][0] = Zsub[0][1] = Zsub[1][0] = Zsub[1][1] = 0;

    for (int i = 0; i < 4; i++) {
        if (!(Z >> i & 1)) continue;
        if (!(i % 2)) {            // 水平边框：bit 由 Y 坐标决定
            int  ybit = (G_pos[i] >= Y + half) ? 1 : 0;
            if (QX1 < X + half)  Zsub[0][ybit] |= (1 << i);
            if (QX2 > X + half) Zsub[1][ybit] |= (1 << i);
        }
        else {            // 垂直边框：bit 由 X 坐标决定
            int  xbit = (G_pos[i] >= X + half) ? 1 : 0;
            if (QY1 < Y + half) Zsub[xbit][0] |= (1 << i);
            if (QY2 > Y + half)   Zsub[xbit][1] |= (1 << i);
        }
    }
}

// ============================================================
//  分化过程中的归一化
// ============================================================

static NR normalize(int T, int B,
    bitmask_t left, bitmask_t right, bitmask_t down, bitmask_t up, bitmask_t h)
{
    NR r; r.Bp = B ^ T; r.left = r.right = r.down = r.up = 0;
    bool B02 = (B == 0 || B == 2);
    switch (T) {
    case 0:r.left = left; r.right = right; r.down = down; r.up = up; break;
    case 1:if (B02) { r.down = left; r.up = right; }
          else { r.left = down; r.right = up; }break;
    case 2:if (B02) { r.left = h - right; r.right = h - left; }
          else { r.down = h - up; r.up = h - down; }break;
    case 3:if (B02) { r.down = h - right; r.up = h - left; }
          else { r.left = h - up; r.right = h - down; }break;
    }return r;
}

static void split(int k, bitmask_t X, bitmask_t Y, int state, bitmask_t precode, int Z) {
    //边框存在性检测
    if (!Z) return;

    //分阶中线存在与出入点获取
    //1.分阶中线识别，2.中线边框端点相对位置与归一化，3.出入点获取与暂存,4.Z更新
    bitmask_t pA = 0, pB = 0, pC = 0;
    bitmask_t h = (bitmask_t)1 << k, half = h / 2, font = precode << (2 * k), unit = (bitmask_t)1 << (2 * k - 2);
    for (int i = 0; i < 4; i++) {
        bitmask_t left = 0, right = 0, down = 0, up = 0;
        //本边框不存在或者阶数未对齐中线
        if (!(Z >> i & 1) || G_k[i] != k) continue;
        //存在一条边框位于中线
        if (!(i % 2)) {//水平边框
            left = max(QX1, X) - X;
            right = min(QX2, X + h) - X;
        }
        else {//垂直边框
            down = max(QY1, Y) - Y;
            up = min(QY2, Y + h) - Y;
        }
        //归一化
        NR nr = normalize(state, i, left, right, down, up, h);
        //出入点获取与暂存,无论分阶中线有几条，在处理某一条分阶中线时，最多只有一条边框能够恰好落在这条中线上并触发交点。过象限位置只有三个位置点，
        switch (nr.Bp) {
        case 0: //上边框
            if (nr.left == 0) pA = font + unit - 1;
            if (nr.right == h) pC = font + 3 * unit;
            break;
        case 2: //下边框
            if (nr.left == 0) pA = font + unit;
            if (nr.right == h) pC = font + 3 * unit - 1;
            break;
        case 1: //右边框
            if (nr.down <= half && nr.up > half) pB = font + 2 * unit - 1;
            break;
        case 3: //左边框
            if (nr.down <= half && nr.up > half) pB = font + 2 * unit;
            break;
        }
        //Z更新
        Z &= ~(1 << i);
    }

    //按Hilbert序分化四个子区域
    // 1.精确 Zsub 分化，2.按 ORDER 表升序递归，3.插入出入口
    // 精确 Zsub 分化
    int Zsub[2][2];
    computeZsub(Z, X, Y, half, Zsub);
    //按 ORDER 表升序递归
    for (int c = 0; c < 4; c++) {
        int xy = ORDER[state][c];
        int x = (xy >> 1) & 1;
        int y = xy & 1;
        int Zc = Zsub[x][y];
        int oldPoint = 0;

        if (Zc) {//该子区域内存在边框
            //只存在单条边框
            if (!(Zc & (Zc - 1))) {
                all_border(EntryExit, lsb32_idx(Zc), precode << 2 | CHM[state][x][y], k - 2, CSM[state][x][y], G_pos[lsb32_idx(Zc)], lsb32_idx(G_pos[lsb32_idx(Zc)]) + 1, oldPoint);
            }
            else {
                split(k - 1, X + x * half, Y + y * half, CSM[state][x][y], (precode << 2) | CHM[state][x][y], Zc);
            }
        }
        //若存在分阶中线上存在出入口，按位置插入
        if (c == 0 && pA) EntryExit.push_back(pA);
        if (c == 1 && pB) EntryExit.push_back(pB);
        if (c == 2 && pC) EntryExit.push_back(pC);
    }
}


Sequence Range_Query(int k, bitmask_t qx1, bitmask_t qy1, bitmask_t qx2, bitmask_t qy2) {

    QX1 = qx1; QX2 = qx2; QY1 = qy1; QY2 = qy2;
    EntryExit.clear();

    bitmask_t H = (bitmask_t)1 << k;
    int oldPoint = 0;

    G_pos[0] = qy2; G_k[0] = min(lsb32_idx(qy2) + 1, k);
    G_pos[1] = qx2; G_k[1] = min(lsb32_idx(qx2) + 1, k);
    G_pos[2] = qy1; G_k[2] = min(lsb32_idx(qy1) + 1, k);
    G_pos[3] = qx1; G_k[3] = min(lsb32_idx(qx1) + 1, k);

    int Z = 0b1111;
    if (qx1 == 0 && qy1 == 0) EntryExit.push_back(0);
    if (qx1 == 0) Z &= ~(1 << 3);
    if (qy1 == 0) Z &= ~(1 << 2);
    if (qx2 == H) Z &= ~(1 << 1);
    if (qy2 == H) Z &= ~(1 << 0);

    if (Z) {
        // 初始调用：先做一次非法检查（只在根节点做一次）
        for (int i = 0; i < 4; i++) {
            if (!(Z >> i & 1)) continue;
            if (!(i % 2)) { if (G_pos[i]<0 || G_pos[i]>H) Z &= ~(1 << i); }
            else { if (G_pos[i]<0 || G_pos[i]>H) Z &= ~(1 << i); }
        }
        if (Z) {
            if (!(Z & (Z - 1))) {//只有一条边,边的类型 lsb32_idx(Z)，边的位置 G_pos[lsb32_idx(Z)]，边对齐的分阶中线阶数lsb32_idx(G_pos[lsb32_idx(Z)])
                all_border(EntryExit, lsb32_idx(Z), 0, k - 1, 0, G_pos[lsb32_idx(Z)], lsb32_idx(G_pos[lsb32_idx(Z)]) + 1, oldPoint);
            }
            else {
                split(k, 0, 0, 0, 0, Z);
            }
        }
    }

    if (qx2 == H && qy1 == 0) EntryExit.push_back(H * H - 1);

    Sequence result;
    //for (int i = 0; i + 1 < (int)EntryExit.size(); i += 2)
    //    result.push_back({ EntryExit[i],EntryExit[i + 1] });
    return result;
}


