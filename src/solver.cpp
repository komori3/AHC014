#define _CRT_NONSTDC_NO_WARNINGS
#define _SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING
#include <bits/stdc++.h>
#include <random>
#include <unordered_set>
#include <array>
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
    int x, y;
    Point(int x = 0, int y = 0) : x(x), y(y) {}
    Point next(int dir) const { return { x + dx8[dir], y + dy8[dir] }; }
    Point next(int dir, int distance) const { return { x + dx8[dir] * distance, y + dy8[dir] * distance }; }
    inline bool operator==(const Point& rhs) const { return x == rhs.x && y == rhs.y; }
    inline bool operator!=(const Point& rhs) const { return !(*this == rhs); }
    inline bool operator<(const Point& rhs) const {
        return y == rhs.y ? x < rhs.x : y < rhs.y;
    }
    string stringify() const {
        return "[" + std::to_string(x) + ", " + std::to_string(y) + "]";
    }
};
std::istream& operator>>(std::istream& in, Point& p) {
    in >> p.x >> p.y;
    return in;
}

inline int weight(int x, int y, int N) {
    int dx = x - N / 2, dy = y - N / 2;
    return dx * dx + dy * dy + 1;
}

inline int weight(const Point& p, int N) {
    return weight(p.x, p.y, N);
}

using Rect = array<Point, 4>;

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

struct Input;
using InputPtr = std::shared_ptr<Input>;
struct Input {
    // 1-indexed
    int N;
    vector<Point> ps;
    int S, Q;
    vector<vector<int>> ws;
    Input(std::istream& in) {
        int M;
        in >> N >> M;
        ps.resize(M);
        in >> ps;
        S = Q = 0;
        ws.resize(N + 2, vector<int>(N + 2, -1));
        for (int y = 0; y < N; y++) for (int x = 0; x < N; x++) {
            ws[y + 1][x + 1] = weight(x, y, N);
            S += ws[y + 1][x + 1];
        }
        for (const auto& [x, y] : ps) Q += ws[y + 1][x + 1];
        for (auto& [x, y] : ps) x++, y++;
    }
};

int debug_count = 0;

struct State;
using StatePtr = std::shared_ptr<State>;
struct State {

    InputPtr input;
    int N;

    vector<vector<bool>> has_point; // 外周は印が付いているとする
    vector<vector<Point>> next_point[8]; // (x,y) から d 方向に進んで初めて印に衝突する点
    vector<vector<array<bool, 8>>> used;

    // 0: left to right
    // 1: bottom-left to top-right
    // 2: top to down
    // 3: bottom-right to top-left
    uint64_t used_bit[4][128];

    vector<std::pair<int, Rect>> cands; // (dir, rect)
    int weight_sum;

    State(InputPtr input) : input(input), N(input->N) {
        has_point.resize(N + 2, vector<bool>(N + 2, true));
        for (int y = 1; y <= N; y++) for (int x = 1; x <= N; x++) has_point[y][x] = false;
        used.resize(N + 2, vector<array<bool, 8>>(N + 2));
        std::memset(used_bit, 0, sizeof(uint64_t) * 4 * 128);
        weight_sum = input->Q;
        for (const auto& [x, y] : input->ps) {
            has_point[y][x] = true;
        }
        for (int d = 0; d < 8; d++) {
            next_point[d].resize(N + 2, vector<Point>(N + 2));
        }
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
    }

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

    int eval() const {
        return (int)round(1e6 * (input->N * input->N) / input->ps.size() * weight_sum / input->S);
    }

    inline bool is_inside(int x, int y) const {
        return 0 < x && x <= N && 0 < y && y <= N;
    }

    inline bool is_inside(const Point& p) const {
        return is_inside(p.x, p.y);
    }

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
        else if (dir == 2) {
            // top to bottom
            uint64_t mask = ((1ULL << p1.y) - 1) ^ ((1ULL << p0.y) - 1);
            used_bit[2][p0.x] |= mask;
        }
        else if (dir == 1) {
            // top-left to bottom-right
            // 0000000
            // 0111111
            // 0122222
            // 0123333
            // 0123444
            // 0123455
            // 0123456
            int r = N + 1 - p0.x + p0.y, c = std::min(p0.x, p0.y);
            uint64_t mask = ((1ULL << (c + dist)) - 1) ^ ((1ULL << c) - 1);
            used_bit[1][r] |= mask;
        }
        else {
            // top-right to bottom-left
            // 0000000
            // 1111110
            // 2222210
            // 3333210
            // 4443210
            // 5543210
            // 6543210
            int r = p0.x + p0.y, c = std::min(N + 1 - p0.x, p0.y);
            uint64_t mask = ((1ULL << (c + dist)) - 1) ^ ((1ULL << c) - 1);
            used_bit[3][r] |= mask;
        }
    }

    bool is_overlap(Point p0, int dir, int dist) const {
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
        else if (dir == 2) {
            // top to bottom
            uint64_t mask = ((1ULL << p1.y) - 1) ^ ((1ULL << p0.y) - 1);
            return used_bit[2][p0.x] & mask;
        }
        else if (dir == 1) {
            // top-left to bottom-right
            int r = N + 1 - p0.x + p0.y, c = std::min(p0.x, p0.y);
            uint64_t mask = ((1ULL << (c + dist)) - 1) ^ ((1ULL << c) - 1);
            return used_bit[1][r] & mask;
        }
        // top-right to bottom-left
        int r = p0.x + p0.y, c = std::min(N + 1 - p0.x, p0.y);
        uint64_t mask = ((1ULL << (c + dist)) - 1) ^ ((1ULL << c) - 1);
        return used_bit[3][r] & mask;
    }

    bool is_valid_rect(const Rect& rect) const {
        //assert(is_valid_rect1(rect) == is_valid_rect2(rect));
        return is_valid_rect2(rect);
    }

    bool is_valid_rect1(const Rect& rect) const {
        // 他の長方形との共通部分が存在しないなら true
        for (int i = 0; i < 4; i++) {
            auto [x, y] = rect[i];
            auto [tx, ty] = rect[(i + 1) & 3];
            int dx = x < tx ? 1 : (x > tx ? -1 : 0);
            int dy = y < ty ? 1 : (y > ty ? -1 : 0);
            int dir = sgn2dir[dy + 1][dx + 1];
            while (x != tx || y != ty) {
                if (used[y][x][dir]) return false;
                x += dx; y += dy;
            }
        }
        return true;
    }

    bool is_valid_rect2(const Rect& rect) const {
        // 他の長方形との共通部分が存在しないなら true
        for (int i = 0; i < 4; i++) {
            auto [x, y] = rect[i];
            auto [tx, ty] = rect[(i + 1) & 3];
            int dx = x < tx ? 1 : (x > tx ? -1 : 0);
            int dy = y < ty ? 1 : (y > ty ? -1 : 0);
            int dir = sgn2dir[dy + 1][dx + 1];
            if (is_overlap(rect[i], dir, std::max(abs(x - tx), abs(y - ty)))) return false;
        }
        return true;
    }

    std::pair<bool, Rect> calc_rect(const Point& p0, int dir0) const {
        assert(!has_point[p0.y][p0.x]);
        int dir1 = (dir0 + 2) & 7;
        // p0 から dir0, dir1 方向に進んで初めて衝突する点を p1, p3 とする
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
            // p1 から p2 方向に進んで初めてぶつかる点は p2 でなければならない
            auto [x, y] = p1.next(dir1);
            if (next_point[dir1][y][x] != p2) return { 0, {} };
        }
        {
            // p3 から p2 方向に進んで初めてぶつかる点は p2 でなければならない
            auto [x, y] = p3.next(dir0);
            if (next_point[dir0][y][x] != p2) return { 0, {} };
        }

        // 4. 他の長方形との共通部分は存在してはいけない
        Rect rect{ p0, p1, p2, p3 };
        if (!is_valid_rect(rect)) return { 0, {} };

        return { input->ws[p0.y][p0.x], rect };
    }

    std::pair<bool, Rect> check_p1(const Point& p1, int d) const {
        // 構成可能なら p0 の座標を返す
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
        if (!is_valid_rect({ p0, p1, p2, p3 })) return { false, {} };
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
        if (!is_valid_rect({ p0, p1, p2, p3 })) return { false, {} };
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
        if (!is_valid_rect({ p0, p1, p2, p3 })) return { false, {} };
        return { true, {p0, p1, p2, p3} };
    }

    void add_cands(const Point& p) {
        // p に印を追加したことで構成できるようになった長方形を調べる
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

    void remove_cands() {
        // invalid になった長方形を弾く
        int new_size = 0;
        for (int i = 0; i < (int)cands.size(); i++) {
            const auto& [d, rect] = cands[i];
            if (has_point[rect[0].y][rect[0].x] || !calc_rect(rect[0], d).first) continue;
            cands[new_size++] = cands[i];
        }
        cands.erase(cands.begin() + new_size, cands.end());
    }

    void add_point(const Point& p) {
        auto [x, y] = p;
        assert(!has_point[y][x]);
        has_point[y][x] = true;
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
        weight_sum += input->ws[rect[0].y][rect[0].x];
        for (int i = 0; i < 4; i++) {
            auto [x, y] = rect[i];
            auto [tx, ty] = rect[(i + 1) % 4];
            int dx = x < tx ? 1 : (x > tx ? -1 : 0);
            int dy = y < ty ? 1 : (y > ty ? -1 : 0);
            int dir = sgn2dir[dy + 1][dx + 1];
            draw_line(rect[i], dir, std::max(abs(x - tx), abs(y - ty)));
            while (x != tx || y != ty) {
                used[y][x][dir] = true;
                x += dx;
                y += dy;
                used[y][x][dir ^ 4] = true;
            }
        }
        remove_cands();
        add_cands(rect[0]);
    }

    // ある印を付けた際の候補点の増減
    int calc_move_score(const Rect& rect) const {
        State state(*this);
        int prev_cands = state.cands.size();
        state.apply_move(rect);
        int now_cands = state.cands.size();
        return now_cands - prev_cands;
    }

    vector<Rect> solve_greedy() {
        vector<Rect> ans;
        while (!cands.empty()) {
            Rect best_rect;
            int best_diff = INT_MIN;
            for (const auto& [_, rect] : cands) {
                int diff = calc_move_score(rect);
                if (chmax(best_diff, diff)) {
                    best_rect = rect;
                }
            }
            apply_move(best_rect);
            ans.push_back(best_rect);
        }
        return ans;
    }

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

std::tuple<int, string, State> compute_score(InputPtr input, const vector<Rect>& out) {
    State state(input);
    for (int t = 0; t < (int)out.size(); t++) {
        const auto& rect = out[t];
        state.apply_move(out[t]);
    }
    int num = 0;
    for (const auto& [x, y] : input->ps) {
        num += input->ws[y][x];
    }
    for (const auto& rect : out) {
        num += input->ws[rect[0].y][rect[0].x];
    }
    int den = 0;
    for (int y = 1; y <= input->N; y++) {
        for (int x = 1; x <= input->N; x++) {
            den += input->ws[y][x];
        }
    }
    int score = (int)round(1e6 * (input->N * input->N) / input->ps.size() * num / den);
    return { score, "", state };
}

struct Output {
    vector<Rect> rects;
    double elapsed_ms;
    Output(const vector<Rect>& rects, double elapsed_ms = -1.0) : rects(rects), elapsed_ms(elapsed_ms) {}
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

    int outer_loop = 0;
    while (timer.elapsed_ms() < 4900) {
        if (timer.elapsed_ms() < 4900) {
            auto state(init_state);
            auto f = [&input, &rnd](const Rect& lhs, const Rect& rhs) {
                return std::make_pair(-input->ws[lhs[0].y][lhs[0].x], rnd.next_int()) < std::make_pair(-input->ws[rhs[0].y][rhs[0].x], rnd.next_int());
            };
            vector<Rect> rects;
            while (true) {
                auto [ok, rect] = state.choose_greedy(f);
                if (!ok) break;
                rects.push_back(rect);
                state.apply_move(rect);
                if (timer.elapsed_ms() > 4900) break;
            }
            if (chmax(best_score, state.eval())) {
                best_rects = rects;
            }
        }
        if (timer.elapsed_ms() < 4900) {
            auto state(init_state);
            auto f = [&input, &rnd](const Rect& lhs, const Rect& rhs) {
                return std::make_pair(area(lhs), rnd.next_int()) < std::make_pair(area(rhs), rnd.next_int());
            };
            vector<Rect> rects;
            while (true) {
                auto [ok, rect] = state.choose_greedy(f);
                if (!ok) break;
                rects.push_back(rect);
                state.apply_move(rect);
                if (timer.elapsed_ms() > 4900) break;
            }
            if (chmax(best_score, state.eval())) {
                best_rects = rects;
            }
        }
        outer_loop++;
    }
    dump(outer_loop);
    return { best_rects, timer.elapsed_ms() };
}

#ifdef _MSC_VER
void batch_test(int seed_begin = 0, int num_seed = 100, int step = 1) {

    constexpr int batch_size = 8;
    const int block_size = batch_size * step;
    int seed_end = seed_begin + num_seed;

    vector<int> scores(num_seed);
#if 1
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
                scores[seed] = std::get<0>(compute_score(input, res.rects));
                cerr << format("seed=%d, score=%d, elapsed_ms=%f\n", seed, scores[seed], res.elapsed_ms);
                mtx.unlock();
            }
            });
    }
#else
    for (int seed = seed_begin; seed < seed_begin + num_seed; seed++) {
        std::ifstream ifs(format("tools/in/%04d.txt", seed));
        std::istream& in = ifs;
        std::ofstream ofs(format("tools/out/%04d.txt", seed));
        std::ostream& out = ofs;

        auto input = std::make_shared<Input>(in);
        Solver solver(input);
        auto ret = solver.solve();
        print_answer(out, ret);

        scores[seed] = calc_score(input, ret);
        cerr << seed << ": " << scores[seed] << ", " << solver.timer.elapsed_ms() << '\n';
    }
#endif

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

        dump(std::get<0>(compute_score(input, output)));
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
    auto ans = solve(input);
    dump(std::get<0>(compute_score(input, ans.rects)));
    out << ans;
    dump(ans.elapsed_ms);
#endif

    return 0;
}