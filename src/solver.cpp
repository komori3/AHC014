#define _CRT_NONSTDC_NO_WARNINGS
#define _SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING
#include <bits/stdc++.h>
#include <random>
#include <unordered_set>
#include <array>
#include <optional>
#ifdef _MSC_VER
#include <opencv2/core.hpp>
#include <opencv2/core/utils/logger.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include <conio.h>
#include <ppl.h>
#include <filesystem>
#include <intrin.h>
//#include <boost/multiprecision/cpp_int.hpp>
int __builtin_clz(unsigned int n)
{
    unsigned long index;
    _BitScanReverse(&index, n);
    return 31 - index;
}
int __builtin_ctz(unsigned int n)
{
    unsigned long index;
    _BitScanForward(&index, n);
    return index;
}
namespace std {
    inline int __lg(int __n) { return sizeof(int) * 8 - 1 - __builtin_clz(__n); }
}
//using __uint128_t = boost::multiprecision::uint128_t;
#else
#pragma GCC target("avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
#endif

/** compro_io **/

/* tuple */
// out
namespace aux {
    template<typename T, unsigned N, unsigned L>
    struct tp {
        static void output(std::ostream& os, const T& v) {
            os << std::get<N>(v) << ", ";
            tp<T, N + 1, L>::output(os, v);
        }
    };
    template<typename T, unsigned N>
    struct tp<T, N, N> {
        static void output(std::ostream& os, const T& v) { os << std::get<N>(v); }
    };
}
template<typename... Ts>
std::ostream& operator<<(std::ostream& os, const std::tuple<Ts...>& t) {
    os << '[';
    aux::tp<std::tuple<Ts...>, 0, sizeof...(Ts) - 1>::output(os, t);
    return os << ']';
}

template<class Ch, class Tr, class Container>
std::basic_ostream<Ch, Tr>& operator<<(std::basic_ostream<Ch, Tr>& os, const Container& x);

/* pair */
// out
template<class S, class T>
std::ostream& operator<<(std::ostream& os, const std::pair<S, T>& p) {
    return os << "[" << p.first << ", " << p.second << "]";
}
// in
template<class S, class T>
std::istream& operator>>(std::istream& is, std::pair<S, T>& p) {
    return is >> p.first >> p.second;
}

/* container */
// out
template<class Ch, class Tr, class Container>
std::basic_ostream<Ch, Tr>& operator<<(std::basic_ostream<Ch, Tr>& os, const Container& x) {
    bool f = true;
    os << "[";
    for (auto& y : x) {
        os << (f ? "" : ", ") << y;
        f = false;
    }
    return os << "]";
}
// in
template <
    class T,
    class = decltype(std::begin(std::declval<T&>())),
    class = typename std::enable_if<!std::is_same<T, std::string>::value>::type
>
std::istream& operator>>(std::istream& is, T& a) {
    for (auto& x : a) is >> x;
    return is;
}

std::ostream& operator<<(std::ostream& os, const std::vector<bool>& v) {
    std::string s(v.size(), ' ');
    for (int i = 0; i < (int)v.size(); i++) s[i] = v[i] + '0';
    os << s;
    return os;
}

/* struct */
template<typename T>
auto operator<<(std::ostream& out, const T& t) -> decltype(out << t.stringify()) {
    out << t.stringify();
    return out;
}

/* setup */
struct IOSetup {
    IOSetup(bool f) {
        if (f) { std::cin.tie(nullptr); std::ios::sync_with_stdio(false); }
        std::cout << std::fixed << std::setprecision(15);
    }
} iosetup(true);

/** string formatter **/
template<typename... Ts>
std::string format(const std::string& f, Ts... t) {
    size_t l = std::snprintf(nullptr, 0, f.c_str(), t...);
    std::vector<char> b(l + 1);
    std::snprintf(&b[0], l + 1, f.c_str(), t...);
    return std::string(&b[0], &b[0] + l);
}

template<typename T>
std::string stringify(const T& x) {
    std::ostringstream oss;
    oss << x;
    return oss.str();
}

/* dump */
#define DUMPOUT std::cerr
std::ostringstream DUMPBUF;
#define dump(...) do{DUMPBUF<<"  ";DUMPBUF<<#__VA_ARGS__<<" :[DUMP - "<<__LINE__<<":"<<__FUNCTION__<<"]"<<std::endl;DUMPBUF<<"    ";dump_func(__VA_ARGS__);DUMPOUT<<DUMPBUF.str();DUMPBUF.str("");DUMPBUF.clear();}while(0);
void dump_func() { DUMPBUF << std::endl; }
template <class Head, class... Tail> void dump_func(Head&& head, Tail&&... tail) { DUMPBUF << head; if (sizeof...(Tail) == 0) { DUMPBUF << " "; } else { DUMPBUF << ", "; } dump_func(std::move(tail)...); }

/* timer */
class Timer {
    double t = 0, paused = 0, tmp;
public:
    Timer() { reset(); }
    static double time() {
#ifdef _MSC_VER
        return __rdtsc() / 3.0e9;
#else
        unsigned long long a, d;
        __asm__ volatile("rdtsc"
            : "=a"(a), "=d"(d));
        return (d << 32 | a) / 3.0e9;
#endif
    }
    void reset() { t = time(); }
    void pause() { tmp = time(); }
    void restart() { paused += time() - tmp; }
    double elapsed_ms() const { return (time() - t - paused) * 1000.0; }
};

/* rand */
struct Xorshift {
    static constexpr uint64_t M = INT_MAX;
    static constexpr double e = 1.0 / M;
    uint64_t x = 88172645463325252LL;
    Xorshift() {}
    Xorshift(uint64_t seed) { reseed(seed); }
    inline void reseed(uint64_t seed) { x = 0x498b3bc5 ^ seed; for (int i = 0; i < 20; i++) next(); }
    inline uint64_t next() { x = x ^ (x << 7); return x = x ^ (x >> 9); }
    inline int next_int() { return next() & M; }
    inline int next_int(int mod) { return next() % mod; }
    inline int next_int(int l, int r) { return l + next_int(r - l + 1); }
    inline double next_double() { return next_int() * e; }
};

/* shuffle */
template<typename T>
void shuffle_vector(std::vector<T>& v, Xorshift& rnd) {
    int n = v.size();
    for (int i = n - 1; i >= 1; i--) {
        int r = rnd.next_int(i);
        std::swap(v[i], v[r]);
    }
}

/* split */
std::vector<std::string> split(std::string str, const std::string& delim) {
    for (char& c : str) if (delim.find(c) != std::string::npos) c = ' ';
    std::istringstream iss(str);
    std::vector<std::string> parsed;
    std::string buf;
    while (iss >> buf) parsed.push_back(buf);
    return parsed;
}

template<typename A, size_t N, typename T> inline void Fill(A(&array)[N], const T& val) {
    std::fill((T*)array, (T*)(array + N), val);
}

template<typename T, typename ...Args> auto make_vector(T x, int arg, Args ...args) { if constexpr (sizeof...(args) == 0)return std::vector<T>(arg, x); else return std::vector(arg, make_vector<T>(x, args...)); }
template<typename T> bool chmax(T& a, const T& b) { if (a < b) { a = b; return true; } return false; }
template<typename T> bool chmin(T& a, const T& b) { if (a > b) { a = b; return true; } return false; }

using ll = long long;
using ld = double;
//using ld = boost::multiprecision::cpp_bin_float_quad;
using pii = std::pair<int, int>;
using pll = std::pair<ll, ll>;

using std::cin, std::cout, std::cerr, std::endl, std::string, std::vector, std::array;



// 5 6 7
// 4 x 0
// 3 2 1
constexpr int dx8[] = { 1,1,0,-1,-1,-1,0,1 };
constexpr int dy8[] = { 0,1,1,1,0,-1,-1,-1 };
constexpr int sgn2dir[3][3] = {
    {5,6,7},
    {4,8,0},
    {3,2,1}
};

struct Point { // TODO: to int8_t
    int x, y;
    Point(int x = 0, int y = 0) : x(x), y(y) {}
    Point next(int dir) const { return { x + dx8[dir], y + dy8[dir] }; }
    Point next(int dir, int distance) const { return { x + dx8[dir] * distance, y + dy8[dir] * distance }; }
    inline bool operator==(const Point& rhs) const { return x == rhs.x && y == rhs.y; }
    inline bool operator!=(const Point& rhs) const { return !(*this == rhs); }
    inline bool operator<(const Point& rhs) const { return y == rhs.y ? x < rhs.x : y < rhs.y; }
    string stringify() const { return "[" + std::to_string(x) + ", " + std::to_string(y) + "]"; }
};
std::istream& operator>>(std::istream& in, Point& p) {
    in >> p.x >> p.y;
    return in;
}

using Rect = array<Point, 4>;

// 面積の小さい順に長方形を置いていくと頑張って 39M くらい
inline int area(const Rect& rect) {
    if (rect[0].x == rect[1].x || rect[0].y == rect[1].y) {
        int d1 = abs(rect[0].x - rect[1].x) + abs(rect[0].y - rect[1].y);
        int d2 = abs(rect[0].x - rect[3].x) + abs(rect[0].y - rect[3].y);
        return d1 * d2;
    }
    int d1 = abs(rect[0].x - rect[1].x) + abs(rect[0].y - rect[1].y);
    int d2 = abs(rect[0].x - rect[3].x) + abs(rect[0].y - rect[3].y);
    return d1 * d2 / 2;
}

// "長方形の描画によって塗られる線の長さが少ない方がよい"と考えると、
// 面積より周長を評価に用いるほうが妥当 -> これで 42M くらい
inline int perimeter(const Rect & rect) {
    // return value is halved
    if (rect[0].x == rect[1].x || rect[0].y == rect[1].y) {
        int d1 = abs(rect[0].x - rect[1].x) + abs(rect[0].y - rect[1].y);
        int d2 = abs(rect[0].x - rect[3].x) + abs(rect[0].y - rect[3].y);
        return d1 + d2;
    }
    int d1 = abs(rect[0].x - rect[1].x);
    int d2 = abs(rect[0].x - rect[3].x);
    return d1 + d2;
}

std::ostream& operator<<(std::ostream& os, const Rect& r) {
    os << format("Rect [p0=(%2d,%2d), p1=(%2d,%2d), p2=(%2d,%2d), p3=(%2d,%2d), half_perimeter=%d]",
        r[0].x, r[0].y, r[1].x, r[1].y, r[2].x, r[2].y, r[3].x, r[3].y, perimeter(r)
        );
    return os;
}

// 1-indexed
struct Input;
using InputPtr = std::shared_ptr<Input>;
struct Input {

    int N;                  // 盤面サイズ (31~61 の奇数)
    vector<Point> ps;       // 印の初期配置
    int S;                  // 盤面全体の重みの総和
    int Q;                  // 印の初期配置の重みの総和
    vector<vector<int>> ws; // 重みをメモした二次元配列

    Input(std::istream& in) {
        int M;
        in >> N >> M;
        ps.resize(M);
        in >> ps;
        S = Q = 0;
        ws.resize(N + 2, vector<int>(N + 2, -1));
        for (int y = 0; y < N; y++) for (int x = 0; x < N; x++) {
            ws[y + 1][x + 1] = weight(x, y);
            S += ws[y + 1][x + 1];
        }
        for (const auto& [x, y] : ps) Q += ws[y + 1][x + 1];
        for (auto& [x, y] : ps) x++, y++;
    }

    inline int weight(int x, int y) const {
        int dx = x - N / 2, dy = y - N / 2;
        return dx * dx + dy * dy + 1;
    }

    inline int weight(const Point& p) const {
        return weight(p.x, p.y);
    }

};

int debug_count = 0;

struct State;
using StatePtr = std::shared_ptr<State>;
struct State {

    InputPtr input;
    int N;

    // 外周は印が付いているとする TODO: uint64_t?
    vector<vector<bool>> has_point; 

    // (x,y) から d 方向に進んで初めて印に衝突する印
    array<array<array<Point, 64>, 64>, 8> next_point;

    // 0: left to right
    // 1: top-left to bottom-right
    // 2: top to down
    // 3: top-right to bottom-left
    array<array<uint64_t, 128>, 4> used_bit;

    // 作れる長方形の候補
    // 便宜上 pair の first に p0->p1 の方向情報も入れている
    // p0->p1->p2->p3 は時計回りで統一する
    vector<std::pair<int, Rect>> cands;

    // 現時点での重みの総和
    int weight_sum;

    // ビームサーチ等することを考慮して、追加した長方形の情報はメンバで持たないようにする

    State(InputPtr input) : input(input), N(input->N) {
        // 外周は true
        has_point.resize(N + 2, vector<bool>(N + 2, true));
        for (int y = 1; y <= N; y++) {
            for (int x = 1; x <= N; x++) {
                has_point[y][x] = false;
            }
        }
        for (const auto& [x, y] : input->ps) {
            has_point[y][x] = true;
        }

        std::memset(used_bit.data(), 0, sizeof(uint64_t) * 4 * 128);

        // 初期化だしナイーブに計算している
        for (int sy = 1; sy <= N; sy++) {
            for (int sx = 1; sx <= N; sx++) {
                for (int d = 0; d < 8; d++) {
                    int y = sy, x = sx;
                    while (!has_point[y][x]) y += dy8[d], x += dx8[d];
                    next_point[d][sy][sx] = { x, y };
                }
            }
        }


        for (int y = 1; y <= input->N; y++) {
            for (int x = 1; x <= input->N; x++) {
                if (has_point[y][x]) continue;
                Point p0(x, y);
                for (int dir = 0; dir < 8; dir++) {
                    auto [ok, rect] = calc_rect(p0, dir);
                    if (!ok) continue;
                    cands.emplace_back(dir, rect);
                }
            }
        }

        weight_sum = input->Q;
    }

    // ナイーブに置ける場所を列挙する　デバッグ用
    vector<std::pair<int, Rect>> enum_cands_naive() const {
        vector<std::pair<int, Rect>> cands;
        for (int y = 1; y <= input->N; y++) {
            for (int x = 1; x <= input->N; x++) {
                if (has_point[y][x]) continue;
                Point p0(x, y);
                for (int dir = 0; dir < 8; dir++) {
                    auto [ok, rect] = calc_rect(p0, dir);
                    if (!ok) continue;
                    cands.emplace_back(dir, rect);
                }
            }
        }
        return cands;
    }

    // 正規化？した得点の計算
    int eval() const {
        return (int)round(1e6 * (input->N * input->N) / input->ps.size() * weight_sum / input->S);
    }

    inline bool is_inside(int x, int y) const {
        return 0 < x && x <= N && 0 < y && y <= N;
    }

    inline bool is_inside(const Point& p) const {
        return is_inside(p.x, p.y);
    }

    // p0 から dir 方向に長さ dist の線を描画する
    void draw_line(Point p0, int dir, int dist) {
        Point p1 = p0.next(dir, dist);
        if (dir >= 4) {
            std::swap(p0, p1);
            dir ^= 4;
        }
        if (dir == 0) {
            // left to right
            uint64_t mask = ((1ULL << p1.x) - 1) ^ ((1ULL << p0.x) - 1);
            used_bit[0][p0.y] |= mask;
        }
        else if (dir == 1) {
            // top-left to bottom-right
            // 6543210
            // 7
            // 8 0000000
            // 9 0111111
            // a 0122222
            // b 0123333
            // c 0123444
            //   0123455
            //   0123456
            int r = N + 1 - p0.x + p0.y, c = std::min(p0.x, p0.y);
            uint64_t mask = ((1ULL << (c + dist)) - 1) ^ ((1ULL << c) - 1);
            used_bit[1][r] |= mask;
        }
        else if (dir == 2) {
            // top to bottom
            uint64_t mask = ((1ULL << p1.y) - 1) ^ ((1ULL << p0.y) - 1);
            used_bit[2][p0.x] |= mask;
        }
        else {
            // top-right to bottom-left
            // 
            //   0123456
            //         7
            // 0000000 8
            // 1111110 9
            // 2222210 a
            // 3333210 b
            // 4443210 c
            // 5543210 
            // 6543210
            int r = p0.x + p0.y, c = std::min(N + 1 - p0.x, p0.y);
            uint64_t mask = ((1ULL << (c + dist)) - 1) ^ ((1ULL << c) - 1);
            used_bit[3][r] |= mask;
        }
    }

    // 長方形の描画
    void draw_rect(const Rect& rect) {
        for (int i = 0; i < 4; i++) {
            auto [x, y] = rect[i];
            auto [tx, ty] = rect[(i + 1) % 4];
            int dx = x < tx ? 1 : (x > tx ? -1 : 0);
            int dy = y < ty ? 1 : (y > ty ? -1 : 0);
            int dir = sgn2dir[dy + 1][dx + 1];
            draw_line(rect[i], dir, std::max(abs(x - tx), abs(y - ty)));
        }
    }

    // p0 から dir 方向に長さ dist の線を考えたときに、既に塗られている区間が存在するか？
    bool is_overlapped(Point p0, int dir, int dist) const {
        Point p1 = p0.next(dir, dist);
        if (dir >= 4) {
            std::swap(p0, p1);
            dir ^= 4;
        }
        if (dir == 0) {
            // left to right
            uint64_t mask = ((1ULL << p1.x) - 1) ^ ((1ULL << p0.x) - 1);
            return used_bit[0][p0.y] & mask;
        }
        else if (dir == 1) {
            // top-left to bottom-right
            int r = N + 1 - p0.x + p0.y, c = std::min(p0.x, p0.y);
            uint64_t mask = ((1ULL << (c + dist)) - 1) ^ ((1ULL << c) - 1);
            return used_bit[1][r] & mask;
        }
        else if (dir == 2) {
            // top to bottom
            uint64_t mask = ((1ULL << p1.y) - 1) ^ ((1ULL << p0.y) - 1);
            return used_bit[2][p0.x] & mask;
        }
        // top-right to bottom-left
        int r = p0.x + p0.y, c = std::min(N + 1 - p0.x, p0.y);
        uint64_t mask = ((1ULL << (c + dist)) - 1) ^ ((1ULL << c) - 1);
        return used_bit[3][r] & mask;
    }

    // 他の長方形との共通部分が存在するなら true
    bool is_overlapped(const Rect& rect) const {
        for (int i = 0; i < 4; i++) {
            auto [x, y] = rect[i];
            auto [tx, ty] = rect[(i + 1) & 3];
            int dx = x < tx ? 1 : (x > tx ? -1 : 0);
            int dy = y < ty ? 1 : (y > ty ? -1 : 0);
            int dir = sgn2dir[dy + 1][dx + 1];
            if (is_overlapped(rect[i], dir, std::max(abs(x - tx), abs(y - ty)))) return true;
        }
        return false;
    }

    std::pair<bool, Rect> calc_rect(const Point& p0, int dir0) const {
        assert(!has_point[p0.y][p0.x]);
        int dir1 = (dir0 + 2) & 7;
        // p0 から dir0, dir1 方向に進んで初めて衝突する印を p1, p3 とする
        // p1 から dir1 方向に伸ばした半直線と p3 から dir0 方向に伸ばした半直線の交点を p2 とする

        // 1. p1, p3 は印が付いていて、外周ではない
        Point p1 = next_point[dir0][p0.y][p0.x];
        if (!is_inside(p1)) return { 0, {} };
        assert(has_point[p1.y][p1.x]);
        Point p3 = next_point[dir1][p0.y][p0.x];
        if (!is_inside(p3)) return { 0, {} };
        assert(has_point[p3.y][p3.x]);

        // 2. p2 に印が付いている
        Point p2(p1.x + p3.x - p0.x, p1.y + p3.y - p0.y);
        if (!is_inside(p2) || !has_point[p2.y][p2.x]) return { 0, {} };

        // 3. 線分 p1-p2, p3-p2 (境界含まない) 上に印は存在してはいけない
        {
            // p1 から p2 方向に進んで初めてぶつかる印は p2 でなければならない
            auto [x, y] = p1.next(dir1);
            if (next_point[dir1][y][x] != p2) return { 0, {} };
        }
        {
            // p3 から p2 方向に進んで初めてぶつかる印は p2 でなければならない
            auto [x, y] = p3.next(dir0);
            if (next_point[dir0][y][x] != p2) return { 0, {} };
        }

        // 4. 他の長方形との共通部分は存在してはいけない
        Rect rect{ p0, p1, p2, p3 };
        if (is_overlapped(rect)) return { 0, {} };

        return { input->ws[p0.y][p0.x], rect };
    }

    std::pair<bool, Rect> check_p1(const Point& p1, int d) const {
        int nd = (d + 2) & 7;
        Point p2;
        {
            auto [x, y] = p1.next(nd);
            p2 = next_point[nd][y][x];
            if (!is_inside(p2)) return { false, {} };
        }
        nd = (nd + 2) & 7;
        Point p3;
        {
            auto [x, y] = p2.next(nd);
            p3 = next_point[nd][y][x];
            if (!is_inside(p3)) return { false, {} };
        }
        Point p0(p1.x + p3.x - p2.x, p1.y + p3.y - p2.y);
        if (!is_inside(p0) || has_point[p0.y][p0.x]) return { false, {} };
        // p0->p1, p0->p3 間に印はない
        {
            auto [x, y] = p0.next(d);
            if (next_point[d][y][x] != p1) return { false, {} };
        }
        {
            auto [x, y] = p0.next((d + 2) & 7);
            if (next_point[(d + 2) & 7][y][x] != p3) return { false, {} };
        }
        if (is_overlapped({ p0, p1, p2, p3 })) return { false, {} };
        return { true, {p0, p1, p2, p3} };
    }

    std::pair<bool, Rect> check_p2(const Point& p2, int d) const {
        int nd = (d + 4) & 7;
        Point p3;
        {
            auto [x, y] = p2.next(nd);
            p3 = next_point[nd][y][x];
            if (!is_inside(p3)) return { false, {} };
        }
        nd = (nd + 2) & 7;
        Point p1;
        {
            auto [x, y] = p2.next(nd);
            p1 = next_point[nd][y][x];
            if (!is_inside(p1)) return { false, {} };
        }
        Point p0(p1.x + p3.x - p2.x, p1.y + p3.y - p2.y);
        if (!is_inside(p0) || has_point[p0.y][p0.x]) return { false, {} };
        // p0->p1, p0->p3 間に印はない
        {
            auto [x, y] = p0.next(d);
            if (next_point[d][y][x] != p1) return { false, {} };
        }
        {
            auto [x, y] = p0.next((d + 2) & 7);
            if (next_point[(d + 2) & 7][y][x] != p3) return { false, {} };
        }
        if (is_overlapped({ p0, p1, p2, p3 })) return { false, {} };
        return { true, {p0, p1, p2, p3} };
    }

    std::pair<bool, Rect> check_p3(const Point& p3, int d) {
        int nd = d;
        Point p2;
        {
            auto [x, y] = p3.next(nd);
            p2 = next_point[nd][y][x];
            if (!is_inside(p2)) return { false, {} };
        }
        nd = (nd + 6) & 7;
        Point p1;
        {
            auto [x, y] = p2.next(nd);
            p1 = next_point[nd][y][x];
            if (!is_inside(p1)) return { false, {} };
        }
        Point p0(p3.x + p1.x - p2.x, p3.y + p1.y - p2.y);
        if (!is_inside(p0) || has_point[p0.y][p0.x]) return { false, {} };
        // p0->p1, p0->p3 間に印はない
        {
            auto [x, y] = p0.next(d);
            if (next_point[d][y][x] != p1) return { false, {} };
        }
        {
            auto [x, y] = p0.next((d + 2) & 7);
            if (next_point[(d + 2) & 7][y][x] != p3) return { false, {} };
        }
        if (is_overlapped({ p0, p1, p2, p3 })) return { false, {} };
        return { true, {p0, p1, p2, p3} };
    }

    // p に印を追加したことで構成できるようになった長方形を調べる
    void add_cands(const Point& p) {
        for (int d = 0; d < 8; d++) {
            bool ok;
            Rect rect;
            std::tie(ok, rect) = check_p1(p, d);
            if (ok) cands.emplace_back(d, rect);
            std::tie(ok, rect) = check_p2(p, d);
            if (ok) cands.emplace_back(d, rect);
            std::tie(ok, rect) = check_p3(p, d);
            if (ok) cands.emplace_back(d, rect);
        }
    }

    // invalid になった長方形を候補から除外する
    void remove_cands() {
        int new_size = 0;
        for (int i = 0; i < (int)cands.size(); i++) {
            const auto& [d, rect] = cands[i];
            if (has_point[rect[0].y][rect[0].x] || !calc_rect(rect[0], d).first) continue;
            cands[new_size++] = cands[i];
        }
        cands.erase(cands.begin() + new_size, cands.end());
    }

    // 印を追加する
    // weight_sum, next_point の更新も行う
    void add_point(const Point& p) {
        auto [x, y] = p;
        assert(!has_point[y][x]);
        has_point[y][x] = true;
        weight_sum += input->ws[p.y][p.x];
        for (int d = 0; d < 8; d++) {
            next_point[d][y][x] = { x, y };
            int rd = d ^ 4;
            int nx = x + dx8[rd], ny = y + dy8[rd];
            while (!has_point[ny][nx]) {
                next_point[d][ny][nx] = { x, y };
                nx += dx8[rd];
                ny += dy8[rd];
            }
        }
    }

    void apply_move(const Rect& rect) {
        add_point(rect[0]);
        draw_rect(rect);
        remove_cands();
        add_cands(rect[0]);
    }

    // pred に従い次に選択すべき長方形を候補から貪欲に選択する
    template<typename F>
    std::pair<bool, Rect> choose_greedy(const F& pred) const {
        if (cands.empty()) return { false, Rect() };
        Rect best = cands.front().second;
        for (int i = 1; i < (int)cands.size(); i++) {
            if (!pred(best, cands[i].second)) best = cands[i].second;
        }
        return { true, best };
    }

};

#ifdef HAVE_OPENCV_HIGHGUI
namespace NVis {

    cv::Scalar scalar(const string& hex) {
        std::stringstream ss;
        ss << std::hex << hex;
        int rgb;
        ss >> rgb;
        return cv::Scalar(rgb % 256, rgb / 256 % 256, rgb / (256 * 256));
    }

    double dist2_seg(double x, double y, double sx, double sy, double tx, double ty) {
        double dx = tx - sx, dy = ty - sy, dx2 = dx * dx, dy2 = dy * dy, r2 = dx2 + dy2;
        double tt = -(dx * (sx - x) + dy * (sy - y));
        if (tt < 0) return (sx - x) * (sx - x) + (sy - y) * (sy - y);
        if (tt > r2) return (tx - x) * (tx - x) + (ty - y) * (ty - y);
        double f1 = dx * (sy - y) - dy * (sx - x);
        return double(f1 * f1) / r2;
    }

    double dist2_rect(double x, double y, const Rect& rect) {
        double res = DBL_MAX;
        for (int i = 0; i < 4; i++) {
            auto [sx, sy] = rect[i];
            auto [tx, ty] = rect[(i + 1) & 3];
            chmin(res, dist2_seg(x, y, sx, sy, tx, ty));
        }
        return res;
    }

    struct MouseParams;
    using MouseParamsPtr = std::shared_ptr<MouseParams>;
    struct MouseParams {
        int pe, px, py, pf;
        int e, x, y, f;
        MouseParams() { e = x = y = f = pe = px = py = pf = INT_MAX; };
        inline void load(int e_, int x_, int y_, int f_) {
            pe = e; px = x; py = y; pf = f;
            e = e_; x = x_; y = y_; f = f_;
        }
        inline bool clicked_left() const { return e == 1 && pe == 0; }
        inline bool clicked_right() const { return e == 2 && pe == 0; }
        inline bool released_left() const { return e == 4; }
        inline bool released_right() const { return e == 5; }
        inline bool drugging_left() const { return e == 0 && f == 1; }
        inline bool drugging_right() const { return e == 0 && f == 2; }
        inline std::pair<int, int> coord() const { return { x, y }; }
        inline std::pair<int, int> displacement() const {
            return { abs(x - px) > 10000 ? 0 : (x - px), abs(y - py) > 10000 ? 0 : (y - py) };
        }
        string stringify() const {
            return format(
                "MouseParams [(e,x,y,f)=(%d,%d,%d,%d), (pe,px,py,pf)=(%d,%d,%d,%d)]"
                , e, x, y, f, pe, px, py, pf
            );
        }
    };

    struct ManualState;
    using ManualStatePtr = std::shared_ptr<ManualState>;
    struct ManualState {

        int frame_id;
        vector<StatePtr> states;
        vector<Rect> rects;
        std::optional<Rect> nearest_rect;

        void remove_successor_states() {
            states.erase(states.begin() + frame_id + 1, states.end());
            rects.erase(rects.begin() + frame_id + 1, rects.end());

        }

    };

    struct ManualSolver {

        const string winname = "img";
        const int board_size = 1000;
        const int grid_size;

        InputPtr input;
        //StatePtr state;

        cv::Mat_<cv::Vec3b> img;

        MouseParamsPtr mp;
        ManualStatePtr msp;

        ManualSolver(InputPtr input)
            : input(input), grid_size(board_size / (input->N + 2)) {
            mp = std::make_shared<MouseParams>();
            msp = std::make_shared<ManualState>();
            msp->frame_id = 0;
            msp->states.push_back(std::make_shared<State>(input));
            Rect sentinel;
            for (int i = 0; i < 4; i++) sentinel[i] = { -2, -2 };
            msp->rects.push_back(sentinel);
        }

        ManualSolver(InputPtr input, const vector<Rect>& rects) 
            : input(input), grid_size(board_size / (input->N + 2)) {
            mp = std::make_shared<MouseParams>();
            msp = std::make_shared<ManualState>();
            msp->frame_id = 0;
            auto state = State(input);
            msp->states.push_back(std::make_shared<State>(state));
            Rect sentinel;
            for (int i = 0; i < 4; i++) sentinel[i] = { -2, -2 };
            msp->rects.push_back(sentinel);
            for (const auto& rect : rects) {
                state.apply_move(rect);
                msp->frame_id++;
                msp->states.push_back(std::make_shared<State>(state));
                msp->rects.push_back(rect);
            }
        }

        inline cv::Point cvt(int x, int y) const {
            return { x * grid_size, y * grid_size };
        }

        inline cv::Point cvt(const Point& p) const {
            return cvt(p.x, p.y);
        }

        inline Point icvt(int x, int y) const {
            return { x / grid_size, y / grid_size };
        }

        std::optional<Rect> calc_nearest_candidate_rect(double thresh = 100.0) const {
            auto state = msp->states[msp->frame_id];
            double min_dist = DBL_MAX;
            Rect selected;
            for (const auto& [_, rect] : state->cands) {
                auto magnified = rect;
                for (auto& p : magnified) p.x *= grid_size, p.y *= grid_size;
                if (chmin(min_dist, dist2_rect(mp->x, mp->y, magnified))) {
                    selected = rect;
                }
            }
            if (min_dist == DBL_MAX || min_dist > thresh) return std::nullopt;
            return selected;
        }

        cv::Mat_<cv::Vec3b> create_board_image() const {

            auto state = msp->states[msp->frame_id];

            cv::Mat_<cv::Vec3b> base_img(board_size, board_size, cv::Vec3b(255, 255, 255));
            
            for (int y = 1; y <= input->N; y++) {
                for (int x = 1; x <= input->N; x++) {
                    cv::circle(base_img, cvt(Point(x, y)), 3, cv::Scalar(150, 150, 150), cv::FILLED);
                }
            }

            for (int i = 1; i <= msp->frame_id; i++) {
                const auto& rect = msp->rects[i];
                for (int i = 0; i < 4; i++) {
                    auto [x1, y1] = rect[i];
                    auto [x2, y2] = rect[(i + 1) & 3];
                    cv::line(base_img, cvt(rect[i]), cvt(rect[(i + 1) & 3]), cv::Scalar(0, 0, 0), 2);
                }
            }

            if (msp->nearest_rect) {
                auto rect = msp->nearest_rect.value();
                for (int i = 0; i < 4; i++) {
                    auto [x1, y1] = rect[i];
                    auto [x2, y2] = rect[(i + 1) & 3];
                    cv::line(base_img, cvt(rect[i]), cvt(rect[(i + 1) & 3]), scalar("0000ff"), 2);
                }
                cv::circle(base_img, cvt(rect[0]), 6, scalar("0000ff"), cv::FILLED);
            }

            for (int y = 1; y <= input->N; y++) {
                for (int x = 1; x <= input->N; x++) {
                    if (state->has_point[y][x]) {
                        cv::circle(base_img, cvt(Point(x, y)), 6, cv::Scalar(0, 0, 0), cv::FILLED);
                    }
                    if (msp->frame_id >= 1) {
                        cv::circle(base_img, cvt(msp->rects[msp->frame_id][0]), 6, cv::Scalar(0, 0, 255), cv::FILLED);
                    }
                }
            }

            auto rect_img = base_img.clone();

            for (const auto& [_, rect] : state->cands) {
                for (int i = 0; i < 4; i++) {
                    auto [x1, y1] = rect[i];
                    auto [x2, y2] = rect[(i + 1) & 3];
                    cv::line(rect_img, cvt(rect[i]), cvt(rect[(i + 1) & 3]), scalar("0000ff"), 2);
                }
            }

            cv::Mat_<cv::Vec3b> img;
            cv::addWeighted(base_img, 0.8, rect_img, 0.2, 0, img);

            return img;
        }

        void vis() {

            cv::namedWindow(winname, cv::WINDOW_AUTOSIZE);
            img = create_board_image();
            cv::imshow(winname, img);

            cv::createTrackbar("frame id", winname, &msp->frame_id, msp->states.size() - 1, frame_callback, this);
            cv::setTrackbarMin("frame id", winname, 0);

            cv::setMouseCallback(winname, mouse_callback, this);

            while (true) {
                int c = cv::waitKey(15);
                if (c == 27 /*ESC*/) break;
                if (c == 8 /*backspace*/) {
                    msp->frame_id = std::max(msp->frame_id - 1, 0);
                    cv::setTrackbarPos("frame id", winname, msp->frame_id);
                }
                
                img = create_board_image();
                cv::imshow(winname, img);
            }

        }

        void fit_frame_trackbar(bool move_newest = false) {
            cv::setTrackbarMax("frame id", winname, msp->states.size() - 1);
            if (move_newest) {
                msp->frame_id = msp->states.size() - 1;
                cv::setTrackbarPos("frame id", winname, msp->frame_id);
            }
        }

        void exec_mouse_action() {
            msp->nearest_rect = calc_nearest_candidate_rect();
            if (auto clicked_left = mp->clicked_left()) {
                if (msp->nearest_rect) {
                    msp->remove_successor_states();
                    auto state = *msp->states.back();
                    auto rect = msp->nearest_rect.value();
                    cerr << rect << '\n';
                    state.apply_move(rect);
                    msp->states.push_back(std::make_shared<State>(state));
                    msp->rects.push_back(rect);
                    fit_frame_trackbar(true);
                    msp->nearest_rect = std::nullopt;
                }
            }
        }

        static void mouse_callback(int e, int x, int y, int f, void* param) {
            auto vis = static_cast<ManualSolver*>(param);
            vis->mp->load(e, x, y, f);
            vis->exec_mouse_action();
        }

        static void frame_callback(int id, void* param) {
            auto vis = static_cast<ManualSolver*>(param);
            vis->msp->frame_id = id;
            cerr << format(
                "frame_id: %d, score: %d -> %d\n",
                id,
                id - 1 >= 0 ? vis->msp->states[id - 1]->eval() : 0,
                vis->msp->states[id]->eval()
            );
        }

    };


    void solve_manual(InputPtr input) {

        ManualSolver msol(input);
        msol.vis();

    }

}
#endif

// 出力用構造体
struct Output {

    vector<Rect> rects;

    // 以下に統計情報を追加してマルチテストケース実行時にサマリとして出力できるようにしている

    int score;
    double elapsed_ms;

    Output(const vector<Rect>& rects, int score, double elapsed_ms = -1.0) : rects(rects), score(score), elapsed_ms(elapsed_ms) {}

    string stringify() const {
        string ans;
        ans += format("%lld\n", rects.size());
        for (const auto& [p0, p1, p2, p3] : rects) {
            ans += format(
                "%d %d %d %d %d %d %d %d\n",
                p0.x - 1, p0.y - 1, p1.x - 1, p1.y - 1, p2.x - 1, p2.y - 1, p3.x - 1, p3.y - 1
            );
        }
        return ans;
    }
};

Output solve(InputPtr input) {

    Timer timer;

    vector<Rect> best_rects;
    int best_score = -1;

    Xorshift rnd;
    State init_state(input);

    // 時間いっぱい回して一番よかったものを採用
    int outer_loop = 0;
    constexpr int timelimit = 4900;
    while (timer.elapsed_ms() < timelimit) {
        if (timer.elapsed_ms() < timelimit) {
            auto state(init_state);
            // 貪欲のタイブレークとして randomness を加えている
            auto f = [&input, &rnd](const Rect& lhs, const Rect& rhs) {
                return 
                    std::make_pair(perimeter(lhs), rnd.next_int())
                    < std::make_pair(perimeter(rhs), rnd.next_int());
            };
            vector<Rect> rects;
            while (true) {
                auto [ok, rect] = state.choose_greedy(f);
                if (!ok) break;
                rects.push_back(rect);
                state.apply_move(rect);
                if (timer.elapsed_ms() > timelimit) break;
            }
            if (chmax(best_score, state.eval())) {
                best_rects = rects;
            }
        }
        outer_loop++;
    }

    //dump(outer_loop);
    return { best_rects, best_score, timer.elapsed_ms() };
}

#ifdef _MSC_VER
// マルチテストケース実行
void batch_test(int seed_begin = 0, int num_seed = 100, int step = 1) {

    constexpr int batch_size = 8;
    const int block_size = batch_size * step;
    int seed_end = seed_begin + num_seed;

    vector<int> scores(num_seed);

    concurrency::critical_section mtx;
    for (int batch_begin = seed_begin; batch_begin < seed_end; batch_begin += block_size) {
        int batch_end = std::min(batch_begin + block_size, seed_end);
        vector<int> seeds;
        for (int seed = batch_begin; seed < batch_end; seed += step) seeds.push_back(seed);
        concurrency::parallel_for_each(seeds.begin(), seeds.end(), [&mtx, &scores](int seed) {
            std::ifstream ifs(format("tools_win/in/%04d.txt", seed));
            std::istream& in = ifs;
            std::ofstream ofs(format("tools_win/out/%04d.txt", seed));
            std::ostream& out = ofs;

            auto input = std::make_shared<Input>(in);
            auto res = solve(input);
            {
                mtx.lock();
                out << res;
                scores[seed] = res.score;
                cerr << format("seed=%d, score=%d, elapsed_ms=%f\n", seed, scores[seed], res.elapsed_ms);
                mtx.unlock();
            }
            });
    }

    dump(std::accumulate(scores.begin(), scores.end(), 0));
}
#endif

void test() {
#ifdef _MSC_VER
    std::ifstream ifs(R"(tools_win\in\0000.txt)");
    std::istream& in = ifs;
#else
    std::istream& in = cin;
#endif
    auto input = std::make_shared<Input>(in);
    string output_raw = R"(20
9 15 12 12 15 15 12 18
15 20 12 17 15 14 18 17
23 22 19 22 19 12 23 12
23 14 22 15 21 14 22 13
10 14 10 13 12 13 12 14
11 11 12 11 12 12 11 12
18 20 15 20 15 19 18 19
19 16 22 19 21 20 18 17
12 19 12 18 15 18 15 19
15 22 12 19 15 16 18 19
14 22 15 22 15 24 14 24
15 8 18 11 15 14 12 11
10 15 9 15 9 14 10 14
11 18 12 19 10 21 9 20
22 23 20 21 21 20 23 22
21 15 18 15 18 14 21 14
15 26 13 24 15 22 17 24
20 20 16 24 14 22 18 18
21 17 18 20 15 17 18 14
11 14 10 13 11 12 12 13
)";
    
    vector<Rect> output;
    {
        std::istringstream iss(output_raw);
        int K;
        iss >> K;
        output.resize(K);
        iss >> output;
        for (auto& [p0, p1, p2, p3] : output) {
            p0.x++, p1.x++, p2.x++, p3.x++;
            p0.y++, p1.y++, p2.y++, p3.y++;
        }
    }

    {
        State state(input);
        for (const auto& rect : output) {
            state.apply_move(rect);
        }
        dump(state.eval());
    }

    {
        State state(input);
        dump(state.calc_rect({ 15, 18 }, 1));
    }

}

int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv) {

#ifdef HAVE_OPENCV_HIGHGUI
    cv::utils::logging::setLogLevel(cv::utils::logging::LogLevel::LOG_LEVEL_SILENT);
#endif

#ifdef _MSC_VER
    std::ifstream ifs(R"(tools_win\in\0005.txt)");
    std::ofstream ofs(R"(tools_win\out\0005.txt)");
    std::istream& in = ifs;
    std::ostream& out = ofs;
#else
    std::istream& in = cin;
    std::ostream& out = cout;
#endif

#if 0
    batch_test();
#else
    auto input = std::make_shared<Input>(in);

#ifdef HAVE_OPENCV_HIGHGUI
    if(false) {
        auto msol = std::make_shared<NVis::ManualSolver>(input);
        msol->vis();
        exit(1);
    }
#endif

    auto ans = solve(input);
    dump(ans.score);
    out << ans;
    dump(ans.elapsed_ms);
#endif

    return 0;
}