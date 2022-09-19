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



constexpr int dx8[] = { 1,1,0,-1,-1,-1,0,1 };
constexpr int dy8[] = { 0,1,1,1,0,-1,-1,-1 };

struct Point {
    int x, y;
    Point(int x = 0, int y = 0) : x(x), y(y) {}
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

struct Input;
using InputPtr = std::shared_ptr<Input>;
struct Input {
    int N;
    vector<Point> ps;
    int S, Q;
    Input(std::istream& in) {
        int M;
        in >> N >> M;
        ps.resize(M);
        in >> ps;
        S = Q = 0;
        for (int y = 0; y < N; y++) for (int x = 0; x < N; x++) S += weight(x, y, N);
        for (const auto& [x, y] : ps) Q += weight(x, y, N);
    }
};

struct State;
using StatePtr = std::shared_ptr<State>;
struct State {

    InputPtr input;
    vector<vector<bool>> has_point;
    vector<vector<array<bool, 8>>> used;
    vector<Rect> rects;
    int weight_sum;

    State(InputPtr input) : input(input) {
        has_point.resize(input->N, std::vector<bool>(input->N, false));
        used.resize(input->N, std::vector<array<bool, 8>>(input->N));
        weight_sum = input->Q;
        for (const auto& [x, y] : input->ps) {
            has_point[y][x] = true;
        }
    }

    int eval() const {
        return (int)round(1e6 * (input->N * input->N) / input->ps.size() * weight_sum / input->S);
    }
    
    string check_move(const Rect& rect) const {
        for (int i = 1; i < 4; i++) {
            if (!has_point[rect[i].y][rect[i].x]) {
                return format("%s does not contain a dot", rect[i].stringify().c_str());
            }
        }
        if (has_point[rect[0].y][rect[0].x]) {
            return format("%s already contains a dot", rect[0].stringify().c_str());
        }
        int dx01 = rect[1].x - rect[0].x, dy01 = rect[1].y - rect[0].y;
        int dx03 = rect[3].x - rect[0].x, dy03 = rect[3].y - rect[0].y;
        if (dx01 * dx03 + dy01 * dy03 != 0) {
            return "Illegal rectangle";
        }
        if (dx01 != 0 && dy01 != 0 && abs(dx01) != abs(dy01)) {
            return "Illegal rectangle";
        }
        if (rect[1].x + dx03 != rect[2].x || rect[1].y + dy03 != rect[2].y) {
            return "Illegal rectangle";
        }
        for (int i = 0; i < 4; i++) {
            auto [x, y] = rect[i];
            auto [tx, ty] = rect[(i + 1) % 4];
            int dx = x < tx ? 1 : (x > tx ? -1 : 0);
            int dy = y < ty ? 1 : (y > ty ? -1 : 0);
            int dir = -1;
            for (dir = 0; dir < 8; dir++) if (dx8[dir] == dx && dy8[dir] == dy) break;
            assert(dir != -1);
            while (x != tx || y != ty) {
                if ((rect[i].x != x || rect[i].y != y) && has_point[y][x]) {
                    return format("There is an obstacle at [%d, %d]", x, y);
                }
                if (used[y][x][dir]) {
                    return "Overlapped rectangles";
                }
                x += dx;
                y += dy;
                if (used[y][x][dir ^ 4]) {
                    return "Overlapped rectangles";
                }
            }
        }
        return "";
    }

    int check_move_fast(const Rect& rect) const {
        for (int i = 1; i < 4; i++) if (!has_point[rect[i].y][rect[i].x]) return 0;
        if (has_point[rect[0].y][rect[0].x]) return 0;
        int dx01 = rect[1].x - rect[0].x, dy01 = rect[1].y - rect[0].y;
        int dx03 = rect[3].x - rect[0].x, dy03 = rect[3].y - rect[0].y;
        if (dx01 * dx03 + dy01 * dy03 != 0) return 0;
        if (dx01 != 0 && dy01 != 0 && abs(dx01) != abs(dy01)) return 0;
        if (rect[1].x + dx03 != rect[2].x || rect[1].y + dy03 != rect[2].y) return 0;
        for (int i = 0; i < 4; i++) {
            auto [x, y] = rect[i];
            auto [tx, ty] = rect[(i + 1) % 4];
            int dx = x < tx ? 1 : (x > tx ? -1 : 0);
            int dy = y < ty ? 1 : (y > ty ? -1 : 0);
            int dir = -1;
            for (dir = 0; dir < 8; dir++) if (dx8[dir] == dx && dy8[dir] == dy) break;
            assert(dir != -1);
            while (x != tx || y != ty) {
                if ((rect[i].x != x || rect[i].y != y) && has_point[y][x]) return 0;
                if (used[y][x][dir]) return 0;
                x += dx;
                y += dy;
                if (used[y][x][dir ^ 4]) return 0;
            }
        }
        return weight(rect[0], input->N);
    }

    void apply_move(const Rect& rect) {
        rects.push_back(rect);
        has_point[rect[0].y][rect[0].x] = true;
        weight_sum += weight(rect[0], input->N);
        for (int i = 0; i < 4; i++) {
            auto [x, y] = rect[i];
            auto [tx, ty] = rect[(i + 1) % 4];
            int dx = x < tx ? 1 : (x > tx ? -1 : 0);
            int dy = y < ty ? 1 : (y > ty ? -1 : 0);
            int dir = -1;
            for (dir = 0; dir < 8; dir++) if (dx8[dir] == dx && dy8[dir] == dy) break;
            assert(dir != -1);
            while (x != tx || y != ty) {
                used[y][x][dir] = true;
                x += dx;
                y += dy;
                used[y][x][dir ^ 4] = true;
            }
        }
    }

    vector<std::pair<int, Rect>> enum_moves_with_weight(int count_limit = INT_MAX) const {
        vector<std::pair<int, Rect>> res;
        vector<Point> ps;
        for (int y = 0; y < input->N; y++) {
            for (int x = 0; x < input->N; x++) {
                if (has_point[y][x]) {
                    ps.emplace_back(x, y);
                }
            }
        }
        for (int idx1 = 0; idx1 + 1 < (int)ps.size(); idx1++) {
            auto p1 = ps[idx1];
            for (int idx3 = idx1 + 1; idx3 < (int)ps.size(); idx3++) {
                auto p3 = ps[idx3];
                if (p1.x != p3.x && p1.y != p3.y) {
                    // not axis-aligned
                    Point p0(p1.x, p3.y), p2(p3.x, p1.y);
                    // 一方のみ true なら置ける
                    if (has_point[p0.y][p0.x] ^ has_point[p2.y][p2.x]) {
                        if (has_point[p0.y][p0.x]) std::swap(p0, p2);
                        // p0: false
                        if (int w = check_move_fast({ p0,p1,p2,p3 })) {
                            res.emplace_back(w, Rect{p0, p1, p2, p3});
                            if (count_limit == (int)res.size()) return res;
                        }
                    }
                }
                int dx = abs(p1.x - p3.x), dy = abs(p1.y - p3.y);
                if (dx != dy && (dx + dy) % 2 == 0) {
                    // not diagonal-aligned
                    // p0: y-y1=x-x1, y-y3=-(x-x3)
                    // y=y1-x1+x=y3+x3-x, 2x=(y3+x3)-(y1-x1)
                    // p2: y-y1=-(x-x1), y-y3=x-x3
                    // y=y1+x1-x=y3-x3+x, 2x=(y1+x1)-(y3-x3)
                    int x0 = (p3.y + p3.x - p1.y + p1.x) / 2, y0 = p1.y - p1.x + x0;
                    if (x0 < 0 || x0 >= input->N || y0 < 0 || y0 >= input->N) continue;
                    int x2 = (p1.y + p1.x - p3.y + p3.x) / 2, y2 = p1.y + p1.x - x2;
                    if (x2 < 0 || x2 >= input->N || y2 < 0 || y2 >= input->N) continue;
                    Point p0(x0, y0), p2(x2, y2);
                    if (has_point[p0.y][p0.x] ^ has_point[p2.y][p2.x]) {
                        if (has_point[p0.y][p0.x]) std::swap(p0, p2);
                        // p0: false
                        if (int w = check_move_fast({ p0,p1,p2,p3 })) {
                            res.emplace_back(w, Rect{ p0, p1, p2, p3 });
                            if (count_limit == (int)res.size()) return res;
                        }
                    }
                }
            }
        }
        return res;
    }

};

std::tuple<int, string, State> compute_score(InputPtr input, const vector<Rect>& out) {
    State state(input);
    for (int t = 0; t < (int)out.size(); t++) {
        const auto& rect = out[t];
        auto err = state.check_move(rect);
        if (!err.empty()) {
            return { 0, format("%s (turn: %d)", err, t), state };
        }
        state.apply_move(out[t]);
    }
    int num = 0;
    for (const auto& [x, y] : input->ps) {
        num += weight(x, y, input->N);
    }
    for (const auto& rect : out) {
        num += weight(rect[0].x, rect[0].y, input->N);
    }
    int den = 0;
    for (int i = 0; i < input->N; i++) {
        for (int j = 0; j < input->N; j++) {
            den += weight(i, j, input->N);
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
            ans += format("%d %d %d %d %d %d %d %d\n", p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, p3.x, p3.y);
        }
        return ans;
    }
};

Output solve(InputPtr input) {
    Timer timer;
    auto state = std::make_shared<State>(input);
    while (true) {
        auto cands = state->enum_moves_with_weight();
        if (cands.empty()) break;
        std::stable_sort(cands.begin(), cands.end(), [](const auto& lhs, const auto& rhs) { return lhs.first > rhs.first; });
        state->apply_move(cands.front().second);
    }
    return { state->rects, timer.elapsed_ms() };
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
    }

    {
        State state(input);
        dump(state.eval());
        auto res = state.enum_moves_with_weight();
        dump(res.size());
    }

    auto [score, err, state] = compute_score(input, output);
    dump(score);
}



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
    dump(std::get<0>(compute_score(input, ans.rects)));
    out << ans;
    dump(ans.elapsed_ms);
#endif

    return 0;
}