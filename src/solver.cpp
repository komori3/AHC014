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

struct Point {
    int8_t x, y;
    Point(int x = 0, int y = 0) : x(x), y(y) {}
    Point next(int dir) const { return { x + dx8[dir], y + dy8[dir] }; }
    Point next(int dir, int distance) const { return { x + dx8[dir] * distance, y + dy8[dir] * distance }; }
    inline bool operator==(const Point& rhs) const { return x == rhs.x && y == rhs.y; }
    inline bool operator!=(const Point& rhs) const { return !(*this == rhs); }
    inline bool operator<(const Point& rhs) const { return y == rhs.y ? x < rhs.x : y < rhs.y; }
    string stringify() const { return "[" + std::to_string((int)x - 1) + ", " + std::to_string((int)y - 1) + "]"; }
};
std::istream& operator>>(std::istream& in, Point& p) {
    int x, y;
    in >> x >> y;
    p.x = x; p.y = y;
    return in;
}

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
    bool operator==(const Rect& rhs) const {
        return data64 == rhs.data64;
    }
    bool operator<(const Rect& rhs) const {
        return data64 < rhs.data64;
    }
    Point operator[](int i) const {
        return { data[i << 1], data[(i << 1) + 1] };
    }
    std::tuple<Point, Point, Point, Point> get_points() const {
        return { {data[0], data[1]}, {data[2], data[3]}, {data[4], data[5]}, {data[6], data[7]} };
    }
};

// ?????????????????????????????????????????????????????????????????? 39M ?????????
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

// "??????????????????????????????????????????????????????????????????????????????"??????????????????
// ?????????????????????????????????????????????????????? -> ????????? 42M ?????????
inline int perimeter(const Rect& rect) {
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
        r[0].x - 1, r[0].y - 1, r[1].x - 1, r[1].y - 1, r[2].x - 1, r[2].y - 1, r[3].x - 1, r[3].y - 1, perimeter(r)
    );
    return os;
}

// 1-indexed
struct Input;
using InputPtr = std::shared_ptr<Input>;
struct Input {

    int N;                          // ??????????????? (31~61 ?????????)
    vector<Point> ps;               // ??????????????????
    int S;                          // ??????????????????????????????
    int Q;                          // ????????????????????????????????????
    array<array<int, 64>, 64> ws;   // ????????????????????????????????????
    array<array<bool, 64>, 64> has_initial_point;

    Input(std::istream& in) {
        int M;
        in >> N >> M;
        ps.resize(M);
        in >> ps;
        S = Q = 0;
        memset(ws.data(), -1, sizeof(int) * 64 * 64);
        memset(has_initial_point.data(), 0, sizeof(bool) * 64 * 64);
        for (int y = 0; y < N; y++) for (int x = 0; x < N; x++) {
            ws[y + 1][x + 1] = weight(x, y);
            S += ws[y + 1][x + 1];
        }
        for (const auto& [x, y] : ps) Q += ws[y + 1][x + 1];
        for (auto& [x, y] : ps) {
            x++, y++;
            has_initial_point[y][x] = true;
        }
    }

    inline int weight(int x, int y) const {
        int dx = x - N / 2, dy = y - N / 2;
        return dx * dx + dy * dy + 1;
    }

    inline int weight(const Point& p) const {
        return weight(p.x, p.y);
    }

};

// ??????????????????
struct Output {

    vector<Rect> rects;

    // ?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????

    int score;
    int outer_loop;
    double elapsed_ms;

    Output(const vector<Rect>& rects, int score, int outer_loop, double elapsed_ms)
        : rects(rects), score(score), outer_loop(outer_loop), elapsed_ms(elapsed_ms) {}

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

int debug_count = 0;

struct State;
using StatePtr = std::shared_ptr<State>;
struct State {

    InputPtr input;
    int N;

    // ??????????????????????????????????????? TODO: uint64_t?
    array<array<bool, 64>, 64> has_point;

    // ????????? p ??????????????????????????? r ?????????????????????
    // 1. r[1], r[2], r[3] ??? p ???????????????????????????????????????????????????????????????
    // 2. p ?????????????????? r ???????????????????????????r[0] ??? p ???????????????????????????????????????????????????????????????
    array<array<int, 64>, 64> num_children;
    array<array<array<Rect, 16>, 64>, 64> children; // 16 ??????????????????????????????????????????
    array<array<int, 64>, 64> num_parents;
    array<array<array<Point, 5>, 64>, 64> parents; // ?????? 5 ???

    // (x,y) ?????? d ????????????????????????????????????????????????
    array<array<array<Point, 64>, 64>, 8> next_point;

    array<array<array<Rect, 8>, 64>, 64> used;

    // 0: left to right
    // 1: top-left to bottom-right
    // 2: top to down
    // 3: top-right to bottom-left
    array<array<uint64_t, 128>, 4> used_bit;

    // ???????????????????????????
    // ????????? pair ??? first ??? p0->p1 ?????????????????????????????????
    // p0->p1->p2->p3 ??????????????????????????????
    int num_cands;
    array<int, 1024> cand_dirs;
    array<Rect, 1024> cand_rects;

    // ??????????????????????????????
    int weight_sum;

    State() {}
    State(InputPtr input) : input(input), N(input->N) {
        // ????????? true
        for (int y = 0; y < 64; y++) {
            for (int x = 0; x < 64; x++) {
                has_point[y][x] = true;
            }
        }
        for (int y = 1; y <= N; y++) {
            for (int x = 1; x <= N; x++) {
                has_point[y][x] = false;
            }
        }
        for (const auto& [x, y] : input->ps) {
            has_point[y][x] = true;
        }

        std::memset(num_children.data(), 0, sizeof(int) * 64 * 64);
        std::memset(children.data(), 0, sizeof(Rect) * 64 * 64 * 16);

        std::memset(num_parents.data(), 0, sizeof(int) * 64 * 64);
        std::memset(parents.data(), 0, sizeof(Point) * 64 * 64 * 5);

        std::memset(used.data(), 0, sizeof(Rect) * 64 * 64 * 8);
        std::memset(used_bit.data(), 0, sizeof(uint64_t) * 4 * 128);

        num_cands = 0;

        // ????????????????????????????????????????????????
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
        num_cands = 0;
        for (int y = 1; y <= input->N; y++) {
            for (int x = 1; x <= input->N; x++) {
                if (!has_point[y][x]) continue;
                Point p(x, y);
                for (int d = 0; d < 8; d++) {
                    Rect rect;
                    rect = check_p1(p, d);
                    if (!rect.data64) continue;
                    cand_dirs[num_cands] = d;
                    cand_rects[num_cands] = rect;
                    num_cands++;
                }
            }
        }
    }

    // ?????????????????????????????????
    int eval() const {
        return (int)round(1e6 * (input->N * input->N) / input->ps.size() * weight_sum / input->S);
    }

    inline bool is_inside(int x, int y) const {
        return 0 < x && x <= N && 0 < y && y <= N;
    }

    inline bool is_inside(const Point& p) const {
        return is_inside(p.x, p.y);
    }

    // p0 ?????? dir ??????????????? dist ?????????????????????
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
            int r = p0.x + p0.y, c = std::min(N + 1 - p0.x, (int)p0.y);
            uint64_t mask = ((1ULL << (c + dist)) - 1) ^ ((1ULL << c) - 1);
            used_bit[3][r] ^= mask;
        }
    }

    // ??????????????????
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

    // p0 ?????? dir ??????????????? dist ?????????????????????????????????????????????????????????????????????????????????
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
        int r = p0.x + p0.y, c = std::min(N + 1 - p0.x, (int)p0.y);
        uint64_t mask = ((1ULL << (c + dist)) - 1) ^ ((1ULL << c) - 1);
        return used_bit[3][r] & mask;
    }

    // ?????????????????????????????????????????????????????? true
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

    Rect check_p0(const Point& p0, int d) const {
        // p0 ?????? dir0, dir1 ????????????????????????????????????????????? p1, p3 ?????????
        // p1 ?????? dir1 ????????????????????????????????? p3 ?????? dir0 ?????????????????????????????????????????? p2 ?????????

        // 1. p1, p3 ?????????????????????????????????????????????
        Point p1 = next_point[d][p0.y][p0.x];
        if (!is_inside(p1)) return 0;

        int nd = (d + 2) & 7;

        Point p3 = next_point[nd][p0.y][p0.x];
        if (!is_inside(p3)) return 0;

        // 2. p2 ????????????????????????
        Point p2(p1.x + p3.x - p0.x, p1.y + p3.y - p0.y);
        if (!is_inside(p2) || !has_point[p2.y][p2.x]) return 0;

        // 3. ?????? p1-p2, p3-p2 (??????????????????) ???????????????????????????????????????
        {
            // p1 ?????? p2 ????????????????????????????????????????????? p2 ???????????????????????????
            auto [x, y] = p1.next(nd);
            if (next_point[nd][y][x] != p2) return 0;
        }
        {
            // p3 ?????? p2 ????????????????????????????????????????????? p2 ???????????????????????????
            auto [x, y] = p3.next(d);
            if (next_point[d][y][x] != p2) return 0;
        }

        // 4. ???????????????????????????????????????????????????????????????
        Rect rect{ p0, p1, p2, p3 };
        if (is_overlapped(rect)) return 0;

        return rect;
    }

    Rect check_p1(const Point& p1, int d) const {

        int nd = (d + 2) & 7;
        Point p2 = p1.next(nd);
        p2 = next_point[nd][p2.y][p2.x];
        if (!is_inside(p2)) return 0;

        nd = (nd + 2) & 7;
        Point p3 = p2.next(nd);
        p3 = next_point[nd][p3.y][p3.x];
        if (!is_inside(p3)) return 0;

        Point p0(p1.x + p3.x - p2.x, p1.y + p3.y - p2.y);
        if (!is_inside(p0) || has_point[p0.y][p0.x]) return 0;
        // p0->p1, p0->p3 ??????????????????
        {
            auto [x, y] = p0.next(d);
            if (next_point[d][y][x] != p1) return 0;
        }
        {
            auto [x, y] = p0.next((d + 2) & 7);
            if (next_point[(d + 2) & 7][y][x] != p3) return 0;
        }
        if (is_overlapped({ p0, p1, p2, p3 })) return 0;
        return { p0, p1, p2, p3 };
    }

    Rect check_p2(const Point& p2, int d) const {
        int nd = (d + 4) & 7;
        Point p3;
        {
            auto [x, y] = p2.next(nd);
            p3 = next_point[nd][y][x];
            if (!is_inside(p3)) return 0;
        }
        nd = (nd + 2) & 7;
        Point p1;
        {
            auto [x, y] = p2.next(nd);
            p1 = next_point[nd][y][x];
            if (!is_inside(p1)) return 0;
        }
        Point p0(p1.x + p3.x - p2.x, p1.y + p3.y - p2.y);
        if (!is_inside(p0) || has_point[p0.y][p0.x]) return 0;
        // p0->p1, p0->p3 ??????????????????
        {
            auto [x, y] = p0.next(d);
            if (next_point[d][y][x] != p1) return 0;
        }
        {
            auto [x, y] = p0.next((d + 2) & 7);
            if (next_point[(d + 2) & 7][y][x] != p3) return 0;
        }
        if (is_overlapped({ p0, p1, p2, p3 })) return 0;
        return { p0, p1, p2, p3 };
    }

    Rect check_p3(const Point& p3, int d) {
        int nd = d;
        Point p2;
        {
            auto [x, y] = p3.next(nd);
            p2 = next_point[nd][y][x];
            if (!is_inside(p2)) return 0;
        }
        nd = (nd + 6) & 7;
        Point p1;
        {
            auto [x, y] = p2.next(nd);
            p1 = next_point[nd][y][x];
            if (!is_inside(p1)) return 0;
        }
        Point p0(p3.x + p1.x - p2.x, p3.y + p1.y - p2.y);
        if (!is_inside(p0) || has_point[p0.y][p0.x]) return 0;
        // p0->p1, p0->p3 ??????????????????
        {
            auto [x, y] = p0.next(d);
            if (next_point[d][y][x] != p1) return 0;
        }
        {
            auto [x, y] = p0.next((d + 2) & 7);
            if (next_point[(d + 2) & 7][y][x] != p3) return 0;
        }
        if (is_overlapped({ p0, p1, p2, p3 })) return 0;
        return { p0, p1, p2, p3 };
    }

    void add_cand(int dir, Rect rect) {
        cand_dirs[num_cands] = dir;
        cand_rects[num_cands] = rect;
        num_cands++;
    }

    // p ????????????????????????????????????????????????????????????????????????????????????
    void add_cands(const Point& p) {
        for (int d = 0; d < 8; d++) {
            Rect rect;
            rect = check_p1(p, d);
            if (rect.data64) add_cand(d, rect);
            rect = check_p2(p, d);
            if (rect.data64) add_cand(d, rect);
            rect = check_p3(p, d);
            if (rect.data64) add_cand(d, rect);
        }
    }

    // invalid ????????????????????????????????????????????????
    void remove_cands() {
        int new_size = 0;
        for (int i = 0; i < num_cands; i++) {
            auto d = cand_dirs[i];
            auto rect = cand_rects[i];
            if (has_point[rect[0].y][rect[0].x] || !check_p0(rect[0], d).data64) {
                continue;
            }
            cand_dirs[new_size] = cand_dirs[i];
            cand_rects[new_size] = cand_rects[i];
            new_size++;
        }
        num_cands = new_size;
    }

    // ??????????????????
    // weight_sum, next_point ??????????????????
    void add_point(const Point& p) {
        auto [x, y] = p;
        has_point[y][x] = true;
        weight_sum += input->ws[p.y][p.x];
        for (int d = 0; d < 8; d++) {
            next_point[d][y][x] = { x, y };
            int rd = d ^ 4;
            int nx = x - dx8[d], ny = y - dy8[d];
            while (!has_point[ny][nx]) {
                next_point[d][ny][nx] = { x, y };
                nx -= dx8[d];
                ny -= dy8[d];
            }
        }
    }

    bool update_dependencies(const Rect& rect) {
        auto [nx, ny] = rect[0];
        for (int i = 0; i < 4; i++) {
            auto [x, y] = rect[i];
            if (input->has_initial_point[y][x]) continue;
            int& numc = num_children[y][x];
            if (numc == 16) return false;
            children[y][x][numc++] = rect;
            if (i) {
                int& nump = num_parents[ny][nx];
                parents[ny][nx][nump++] = { x, y };
            }
        }
        // ???????????????????????????????????????
        for (int dir = 0; dir < 4; dir++) {
            auto r1 = used[ny][nx][dir];
            auto r2 = used[ny][nx][dir ^ 4];
            if (r1.data64 && r1.data64 == r2.data64) {
                auto [x, y] = r1[0];
                int& numc = num_children[y][x];
                if (numc == 16) return false;
                children[y][x][numc++] = rect;
                int& nump = num_parents[ny][nx];
                parents[ny][nx][nump++] = { x, y };
            }
        }
        return true;
    }

    bool apply_move(const Rect& rect) {
        add_point(rect[0]);
        if (!update_dependencies(rect)) return false;
        toggle_rect(rect);
        remove_cands();
        add_cands(rect[0]);
        return true;
    }

    // ????????????????????????????????????????????????
    void remove_point_dfs(const Point& p) {
        auto [x, y] = p;
        // 0 ???????????????????????????????????????
        for (int i = 1; i < num_children[y][x]; i++) {
            auto& rect = children[y][x][i];
            auto np = rect[0];
            auto [nx, ny] = np;
            rect.data64 = 0;
            if (nx == 0 || !has_point[ny][nx] || input->has_initial_point[ny][nx]) continue;
            remove_point_dfs(np);
        }
        auto& rect = children[y][x][0];
        for (int i = 0; i < num_parents[y][x]; i++) {
            auto [px, py] = parents[y][x][i];
            auto& pcs = children[py][px];
            int& numc = num_children[py][px];
            for (int j = 0; j < numc; j++) {
                if (pcs[j] == rect) {
                    for (int k = j + 1; k < numc; k++) {
                        pcs[k - 1] = pcs[k];
                    }
                    numc--;
                    break;
                }
            }
        }
        has_point[y][x] = false;
        weight_sum -= input->ws[y][x];
        for (int d = 0; d < 8; d++) {
            int tx = x, ty = y;
            while (!has_point[ty][tx]) tx += dx8[d], ty += dy8[d];
            int nx = tx - dx8[d], ny = ty - dy8[d];
            while (!has_point[ny][nx]) {
                next_point[d][ny][nx] = { tx, ty };
                nx -= dx8[d];
                ny -= dy8[d];
            }
        }
        toggle_rect(rect);
        rect.data64 = 0;
        num_children[y][x] = 0;
        num_parents[y][x] = 0;
    }

    void remove_point(const Point& p) {
        remove_point_dfs(p);
        update_cands();
    }

    void random_remove(Xorshift& rnd) {
        vector<Point> ps;
        for (int y = 1; y <= input->N; y++) {
            for (int x = 1; x <= input->N; x++) {
                if (has_point[y][x] && !input->has_initial_point[y][x]) {
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
                if (!has_point[y][x] || input->has_initial_point[y][x]) continue;
                auto u = children[y][x][0].data64;
                r2i[u] = i2r.size();
                i2r.push_back(u);
            }
        }
        int V = r2i.size();
        vector<vector<int>> graph(V);
        vector<int> indeg(V);
        for (int y = 1; y <= input->N; y++) {
            for (int x = 1; x <= input->N; x++) {
                if (!has_point[y][x] || input->has_initial_point[y][x]) continue;
                auto u64 = children[y][x][0].data64;
                auto u = r2i[u64];
                for (int i = 1; i < num_children[y][x]; i++) {
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
                indeg[v]--;
                if (indeg[v] == 0) qu.push(v);
            }
        }
        return Output(ans, eval(), -1, -1.0);
    }

    // pred ???????????????????????????????????????????????????????????????????????????
    template<typename F>
    Rect choose_greedy(const F& pred) const {
        if (!num_cands) return 0;
        Rect best = cand_rects[0];
        for (int i = 1; i < num_cands; i++) {
            if (!pred(best, cand_rects[i])) best = cand_rects[i];
        }
        return best;
    }

    template<typename F>
    bool solve_greedy(const F& pred) {
        while (num_cands) {
            auto rect = choose_greedy(pred);
            auto ok = apply_move(rect);
            if (!ok) return false;
        }
        return true;
    }

};



Output solve(InputPtr input) {

    Timer timer;
    Xorshift rnd;

    // ???????????????????????????????????? randomness ??????????????????
    auto f = [&input, &rnd](const Rect& lhs, const Rect& rhs) {
        return
            std::make_pair(perimeter(lhs), rnd.next_int())
            < std::make_pair(perimeter(rhs), rnd.next_int());
    };

    State init_state(input);

    State state;
    while (true) {
        state = init_state;
        bool ok = state.solve_greedy(f);
        if (ok) break;
    }

    State best_state(state);
    int best_score = best_state.eval();

    auto get_temp = [](double startTemp, double endTemp, double t, double T) {
        return endTemp + (startTemp - endTemp) * (T - t) / T;
    };

    // ??????
    int outer_loop = 0;
    constexpr double end_time = 19900;
    double start_time = timer.elapsed_ms(), now_time;
    while ((now_time = timer.elapsed_ms()) < end_time) {

        auto nstate(state);
        int prev_score = nstate.eval();
        nstate.random_remove(rnd);
        bool ok = nstate.solve_greedy(f);
        if (!ok) continue;
        int now_score = nstate.eval();

        int diff = now_score - prev_score;
        double temp = get_temp(20000.0, 0.0, now_time - start_time, end_time - start_time);
        double prob = exp(diff / temp);

        if (rnd.next_double() < prob) {
            state = nstate;
            if (chmax(best_score, now_score)) {
                best_state = state;
                dump(best_score);
            }
        }
        outer_loop++;
    }
    dump(outer_loop);

    auto output = best_state.create_output();
    output.elapsed_ms = timer.elapsed_ms();
    output.outer_loop = outer_loop;

    return output;
}

#ifdef _MSC_VER
// ?????????????????????????????????
void batch_test(int seed_begin = 0, int num_seed = 100, int step = 1) {

    constexpr int batch_size = 5;
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
                cerr << format("seed=%d, score=%d, outer_loop=%d, elapsed_ms=%f\n", seed, scores[seed], res.outer_loop, res.elapsed_ms);
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
    std::ifstream ifs(R"(tools_win\in\0000.txt)");
    std::ofstream ofs(R"(tools_win\out\0000.txt)");
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
    dump(ans.score, ans.outer_loop, ans.elapsed_ms);
    out << ans;
#endif

    return 0;
}