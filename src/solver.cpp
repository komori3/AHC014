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

//using Rect = array<Point, 4>;

struct Rect {
    union {
        int8_t data[8];
        uint64_t data64;
    };
    Rect() : data64(0) {}
    Rect(uint64_t data64) : data64(data64) {}
    Rect(const Point& p0, const Point& p1, const Point& p2, const Point& p3) {
        data[0] = p0.x; data[1] = p0.y;
        data[2] = p1.x; data[3] = p1.y;
        data[4] = p2.x; data[5] = p2.y;
        data[6] = p3.x; data[7] = p3.y;
    }
    Point operator[](int i) const {
        return { data[i << 1], data[(i << 1) + 1] };
    }
    std::tuple<Point, Point, Point, Point> get_points() const {
        return { {data[0], data[1]}, {data[2], data[3]}, {data[4], data[5]}, {data[6], data[7]} };
    }
};

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
        r[0].x-1, r[0].y-1, r[1].x-1, r[1].y-1, r[2].x-1, r[2].y-1, r[3].x-1, r[3].y-1, perimeter(r)
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

// 出力用構造体
struct Output {

    vector<Rect> rects;

    // 以下に統計情報を追加してマルチテストケース実行時にサマリとして出力できるようにしている

    int score;
    double elapsed_ms;
    int max_cands;

    Output(const vector<Rect>& rects, int score, int max_cands, double elapsed_ms = -1.0)
        : rects(rects), score(score), elapsed_ms(elapsed_ms), max_cands(max_cands) {}

    string stringify() const {
        string ans;
        ans += format("%lld\n", rects.size());
        for (const auto& rect : rects) {
            auto [p0, p1, p2, p3] = rect.get_points();
            ans += format(
                "%d %d %d %d %d %d %d %d\n",
                p0.x - 1, p0.y - 1, p1.x - 1, p1.y - 1, p2.x - 1, p2.y - 1, p3.x - 1, p3.y - 1
            );
        }
        return ans;
    }
};

struct State;
using StatePtr = std::shared_ptr<State>;
struct State {

    InputPtr input;
    int N;

    // 外周は印が付いているとする TODO: uint64_t?
    array<array<bool, 64>, 64> has_initial_point; // TODO: input に回していい
    array<array<bool, 64>, 64> has_point;

    // ある点 p に印を付けて長方形 r を描画したとき
    // 1. r[1], r[2], r[3] は p よりも前に印を付けられていなければならない
    // 2. p が他の長方形 r の辺上にあるとき、r[0] は p よりも前に印を付けられていなければならない（高々 2 つ)
    array<array<int, 64>, 64> num_children;
    array<array<array<Rect, 32>, 64>, 64> children;

    // (x,y) から d 方向に進んで初めて印に衝突する印
    array<array<array<Point, 64>, 64>, 8> next_point;

    array<array<array<Rect, 8>, 64>, 64> used;

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

    State(InputPtr input) : input(input), N(input->N) {
        // 外周は true
        for (int y = 0; y < 64; y++) {
            for (int x = 0; x < 64; x++) {
                has_initial_point[y][x] = true;
                has_point[y][x] = true;
            }
        }
        for (int y = 1; y <= N; y++) {
            for (int x = 1; x <= N; x++) {
                has_initial_point[y][x] = false;
                has_point[y][x] = false;
            }
        }
        for (const auto& [x, y] : input->ps) {
            has_initial_point[y][x] = true;
            has_point[y][x] = true;
        }

        std::memset(num_children.data(), 0, sizeof(int) * 64 * 64);
        std::memset(children.data(), 0, sizeof(Rect) * 64 * 64 * 4);

        std::memset(used.data(), 0, sizeof(Rect) * 64 * 64 * 8);
        std::memset(used_bit.data(), 0, sizeof(uint64_t) * 4 * 128);

        // 初期化だしナイーブに計算している
        update_next_point();
        update_cands();

        weight_sum = input->Q;

    }

    void update_next_point() {
        for (int sy = 1; sy <= N; sy++) {
            for (int sx = 1; sx <= N; sx++) {
                for (int d = 0; d < 8; d++) {
                    int y = sy, x = sx;
                    while (!has_point[y][x]) y += dy8[d], x += dx8[d];
                    next_point[d][sy][sx] = { x, y };
                }
            }
        }
    }

    void update_cands() {
        cands.clear();
        for (int y = 1; y <= input->N; y++) {
            for (int x = 1; x <= input->N; x++) {
                if (has_point[y][x]) continue;
                Point p0(x, y);
                for (int dir = 0; dir < 8; dir++) {
                    auto [ok, rect] = check_p0(p0, dir);
                    if (!ok) continue;
                    cands.emplace_back(dir, rect);
                }
            }
        }
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
    void toggle_line(Point p0, int dir, int dist) {
        Point p1 = p0.next(dir, dist);
        if (dir >= 4) {
            std::swap(p0, p1);
            dir ^= 4;
        }
        if (dir == 0) {
            // left to right
            uint64_t mask = ((1ULL << p1.x) - 1) ^ ((1ULL << p0.x) - 1);
            used_bit[0][p0.y] ^= mask;
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
            used_bit[1][r] ^= mask;
        }
        else if (dir == 2) {
            // top to bottom
            uint64_t mask = ((1ULL << p1.y) - 1) ^ ((1ULL << p0.y) - 1);
            used_bit[2][p0.x] ^= mask;
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
            used_bit[3][r] ^= mask;
        }
    }

    // 長方形の描画
    void toggle_rect(const Rect& rect) {
        for (int i = 0; i < 4; i++) {
            auto [x, y] = rect[i];
            auto [tx, ty] = rect[(i + 1) % 4];
            int dx = x < tx ? 1 : (x > tx ? -1 : 0);
            int dy = y < ty ? 1 : (y > ty ? -1 : 0);
            int dir = sgn2dir[dy + 1][dx + 1];
            toggle_line(rect[i], dir, std::max(abs(x - tx), abs(y - ty)));
            while (x != tx || y != ty) {
                used[y][x][dir].data64 ^= rect.data64;
                x += dx;
                y += dy;
                used[y][x][dir ^ 4].data64 ^= rect.data64;
            }
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

    std::pair<bool, Rect> check_p0(const Point& p0, int d) const {
        assert(!has_point[p0.y][p0.x]);
        int nd = (d + 2) & 7;
        // p0 から dir0, dir1 方向に進んで初めて衝突する印を p1, p3 とする
        // p1 から dir1 方向に伸ばした半直線と p3 から dir0 方向に伸ばした半直線の交点を p2 とする

        // 1. p1, p3 は印が付いていて、外周ではない
        Point p1 = next_point[d][p0.y][p0.x];
        if (!is_inside(p1)) return { 0, {} };
        assert(has_point[p1.y][p1.x]);
        Point p3 = next_point[nd][p0.y][p0.x];
        if (!is_inside(p3)) return { 0, {} };
        assert(has_point[p3.y][p3.x]);

        // 2. p2 に印が付いている
        Point p2(p1.x + p3.x - p0.x, p1.y + p3.y - p0.y);
        if (!is_inside(p2) || !has_point[p2.y][p2.x]) return { 0, {} };

        // 3. 線分 p1-p2, p3-p2 (境界含まない) 上に印は存在してはいけない
        {
            // p1 から p2 方向に進んで初めてぶつかる印は p2 でなければならない
            auto [x, y] = p1.next(nd);
            if (next_point[nd][y][x] != p2) return { 0, {} };
        }
        {
            // p3 から p2 方向に進んで初めてぶつかる印は p2 でなければならない
            auto [x, y] = p3.next(d);
            if (next_point[d][y][x] != p2) return { 0, {} };
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
            if (ok) {
                cands.emplace_back(d, rect);
            }
            std::tie(ok, rect) = check_p2(p, d);
            if (ok) {
                cands.emplace_back(d, rect);
            }
            std::tie(ok, rect) = check_p3(p, d);
            if (ok) {
                cands.emplace_back(d, rect);
            }
        }
    }

    // invalid になった長方形を候補から除外する
    void remove_cands() {
        int new_size = 0;
        for (int i = 0; i < (int)cands.size(); i++) {
            const auto& [d, rect] = cands[i];
            if (has_point[rect[0].y][rect[0].x] || !check_p0(rect[0], d).first) {
                continue;
            }
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
        auto [nx, ny] = rect[0];
        for (int i = 0; i < 4; i++) {
            auto [x, y] = rect[i];
            if (has_initial_point[y][x]) continue;
            auto& nc = num_children[y][x];
            children[y][x][nc++] = rect;
            assert(nc <= 32);
        }
        for (int dir = 0; dir < 4; dir++) {
            auto r1 = used[ny][nx][dir];
            auto r2 = used[ny][nx][dir ^ 4];
            if (r1.data64 && r1.data64 == r2.data64) {
                auto [x, y] = r1[0];
                auto& nc = num_children[y][x];
                children[y][x][nc++] = rect;
                assert(nc <= 32);
            }
        }
        toggle_rect(rect);
        remove_cands();
        add_cands(rect[0]);
    }

    // 全ての点から削除しないといけない
    void remove_point_dfs(const Point& p, vector<Rect>& del) {
        auto [x, y] = p;
        assert(has_point[y][x] && !has_initial_point[y][x]);
        // 0 番目は自身なので最後に消す
        for (int i = 1; i < num_children[y][x]; i++) {
            auto& rect = children[y][x][i];
            auto np = rect[0];
            auto [nx, ny] = np;
            rect.data64 = 0;
            if (nx == 0 || !has_point[ny][nx] || has_initial_point[ny][nx]) continue;
            remove_point_dfs(np, del);
        }
        auto& rect = children[y][x][0];
        del.push_back(rect);
        has_point[y][x] = false;
        weight_sum -= input->ws[y][x];
        toggle_rect(rect);
        rect.data64 = 0;
        num_children[y][x] = 0;
    }

    void remove_rect_naive(const Rect& rect) {
        for (int y = 1; y <= input->N; y++) {
            for (int x = 1; x <= input->N; x++) {
                if (!has_point[y][x] || has_initial_point[y][x]) continue;
                auto& numc = num_children[y][x];
                for (int i = 0; i < numc; i++) {
                    if (rect.data64 == children[y][x][i].data64) {
                        for (int j = i + 1; j < numc; j++) {
                            children[y][x][j - 1] = children[y][x][j];
                        }
                        children[y][x][numc - 1].data64 = 0;
                        numc--;
                        break;
                    }
                }
            }
        }
    }

    void remove_point(const Point& p) {
        vector<Rect> del;
        remove_point_dfs(p, del);
        for (const auto& r : del) {
            remove_rect_naive(r);
        }
        // update next_point, cands
        update_next_point();
        update_cands();
    }

    void random_remove(Xorshift& rnd) {
        vector<Point> ps;
        for (int y = 1; y <= input->N; y++) {
            for (int x = 1; x <= input->N; x++) {
                if (has_point[y][x] && !has_initial_point[y][x]) {
                    ps.emplace_back(x, y);
                }
            }
        }
        auto p = ps[rnd.next_int(ps.size())];
        remove_point(p);
    }

    Output create_output() const {
        // topological sort
        std::map<uint64_t, int> r2i;
        vector<uint64_t> i2r;
        for (int y = 1; y <= input->N; y++) {
            for (int x = 1; x <= input->N; x++) {
                if (!has_point[y][x] || has_initial_point[y][x]) continue;
                auto u = children[y][x][0].data64;
                assert(!r2i.count(u));
                r2i[u] = i2r.size();
                i2r.push_back(u);
            }
        }
        int V = r2i.size();
        vector<vector<int>> graph(V);
        vector<int> indeg(V);
        for (int y = 1; y <= input->N; y++) {
            for (int x = 1; x <= input->N; x++) {
                if (!has_point[y][x] || has_initial_point[y][x]) continue;
                auto u64 = children[y][x][0].data64;
                auto u = r2i[u64];
                auto numc = num_children[y][x];
                for (int i = 1; i < numc; i++) {
                    auto v64 = children[y][x][i].data64;
                    auto v = r2i[v64];
                    graph[u].push_back(v);
                    indeg[v]++;
                }
            }
        }
        std::queue<uint64_t> qu;
        for (int u = 0; u < V; u++) {
            if (!indeg[u]) qu.push(u);
        }
        vector<Rect> ans;
        while (!qu.empty()) {
            auto u = qu.front(); qu.pop();
            Rect r; r.data64 = i2r[u];
            ans.push_back(r);
            for (auto v : graph[u]) {
                assert(indeg[v] > 0);
                indeg[v]--;
                if (indeg[v] == 0) qu.push(v);
            }
        }
        return Output(ans, eval(), -1, -1.0);
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



Output solve(InputPtr input) {

    Timer timer;
    Xorshift rnd;

    // 貪欲のタイブレークとして randomness を加えている
    auto f = [&input, &rnd](const Rect& lhs, const Rect& rhs) {
        return
            std::make_pair(perimeter(lhs), rnd.next_int())
            < std::make_pair(perimeter(rhs), rnd.next_int());
    };

    State state(input);

    // まずは貪欲
    while (!state.cands.empty()) {
        auto [ok, rect] = state.choose_greedy(f);
        if (!ok) break;
        state.apply_move(rect);
    }
    int best_score = state.eval();
    dump(best_score);

    // 山を登る
    int outer_loop = 0;
    constexpr int timelimit = 4900;
    while (timer.elapsed_ms() < timelimit) {
        auto nstate(state);
        nstate.random_remove(rnd);
        while (!nstate.cands.empty()) {
            auto [ok, rect] = nstate.choose_greedy(f);
            if (!ok) break;
            nstate.apply_move(rect);
        }
        if (chmax(best_score, nstate.eval())) {
            state = nstate;
            dump(best_score);
        }
        outer_loop++;
    }
    dump(outer_loop);

    auto output = state.create_output();
    output.elapsed_ms = timer.elapsed_ms();

    return output;
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
                cerr << format("seed=%d, score=%d, elapsed_ms=%f, max_cands=%d\n", seed, scores[seed], res.elapsed_ms, res.max_cands);
                mtx.unlock();
            }
            });
    }

    dump(std::accumulate(scores.begin(), scores.end(), 0));
}
#endif



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
    auto ans = solve(input);
    dump(ans.score);
    out << ans;
    dump(ans.elapsed_ms);
#endif

    return 0;
}