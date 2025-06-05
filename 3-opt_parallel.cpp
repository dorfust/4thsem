#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <chrono>
#ifdef _OPENMP
  #include <omp.h>
#endif

namespace fs = std::filesystem;
using Clock = std::chrono::high_resolution_clock;

// Представление точки (x, y)
struct Coord {
    double x, y;
};

// Загрузка точек из файла TSP
static std::vector<Coord> readTSP(const fs::path& path) {
    std::ifstream in(path);
    if (!in) {
        std::cerr << "Failed to open: " << path << "\n";
        std::exit(1);
    }
    int n;
    in >> n;
    std::vector<Coord> pts(n);
    for (int i = 0; i < n; ++i) {
        in >> pts[i].x >> pts[i].y;
    }
    return pts;
}

// Простая реализация 2-opt для локального улучшения
static void apply2Opt(const std::vector<std::vector<double>>& dist, std::vector<int>& tour) {
    int n = (int)tour.size();
    bool improved = true;
    while (improved) {
        improved = false;
        for (int i = 0; i < n - 1 && !improved; ++i) {
            for (int j = i + 2; j < n && !improved; ++j) {
                int a = tour[i], b = tour[i + 1];
                int c = tour[j], d = tour[(j + 1) % n];
                double before = dist[a][b] + dist[c][d];
                double after  = dist[a][c] + dist[b][d];
                if (after + 1e-9 < before) {
                    std::reverse(tour.begin() + i + 1, tour.begin() + j + 1);
                    improved = true;
                }
            }
        }
    }
}

// Полный 3-opt с параллельным поиском лучшей перемены
static double solve3Opt(const std::vector<Coord>& pts) {
    int n = (int)pts.size();
    // Построение матрицы расстояний
    std::vector<std::vector<double>> D(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double dx = pts[i].x - pts[j].x;
            double dy = pts[i].y - pts[j].y;
            D[i][j] = D[j][i] = std::sqrt(dx * dx + dy * dy);
        }
    }
    // Начальный тур – ближайший сосед + 2-opt
    std::vector<int> tour;
    tour.reserve(n);
    std::vector<bool> used(n, false);
    int current = 0;
    tour.push_back(current);
    used[current] = true;
    for (int step = 1; step < n; ++step) {
        int next = -1;
        double bestDist = 1e18;
        for (int cand = 0; cand < n; ++cand) {
            if (!used[cand] && D[current][cand] < bestDist) {
                bestDist = D[current][cand];
                next = cand;
            }
        }
        tour.push_back(next);
        used[next] = true;
        current = next;
    }
    apply2Opt(D, tour);

    // Построение списков K ближайших соседей
    const int K = 10;
    std::vector<std::vector<int>> nearList(n);
    for (int i = 0; i < n; ++i) {
        std::vector<std::pair<double,int>> tmp;
        tmp.reserve(n - 1);
        for (int j = 0; j < n; ++j) {
            if (j != i) {
                tmp.emplace_back(D[i][j], j);
            }
        }
        std::sort(tmp.begin(), tmp.end());
        int limit = std::min(K, (int)tmp.size());
        nearList[i].reserve(limit);
        for (int t = 0; t < limit; ++t) {
            nearList[i].push_back(tmp[t].second);
        }
    }

    bool anyChange = true;
    while (anyChange) {
        anyChange = false;
        double bestDelta = 0.0;
        int bestI = -1, bestJ = -1, bestK = -1;

        #pragma omp parallel
        {
            double localDelta = 0.0;
            int locI = -1, locJ = -1, locK = -1;

            #pragma omp for schedule(dynamic)
            for (int i = 0; i < n; ++i) {
                for (int jj = 0; jj < (int)nearList[i].size(); ++jj) {
                    int j = nearList[i][jj];
                    if (j <= i) continue;
                    for (int kk = 0; kk < (int)nearList[i].size(); ++kk) {
                        int k = nearList[i][kk];
                        if (k <= j) continue;
                        int a = tour[i], b = tour[i + 1];
                        int c = tour[j], d = tour[j + 1];
                        int e = tour[k], f = tour[(k + 1) % n];
                        double currLen = D[a][b] + D[c][d] + D[e][f];
                        double altLen  = D[a][c] + D[b][e] + D[d][f];
                        double delta = altLen - currLen;
                        if (delta < localDelta) {
                            localDelta = delta;
                            locI = i; locJ = j; locK = k;
                        }
                    }
                }
            }

            #pragma omp critical
            {
                if (localDelta < bestDelta) {
                    bestDelta = localDelta;
                    bestI = locI; bestJ = locJ; bestK = locK;
                }
            }
        }

        if (bestDelta < 0.0) {
            anyChange = true;
            std::reverse(tour.begin() + bestI + 1, tour.begin() + bestJ + 1);
            std::reverse(tour.begin() + bestJ + 1, tour.begin() + bestK + 1);
        }
    }

    double totalLen = 0.0;
    for (int i = 0; i < n; ++i) {
        totalLen += D[tour[i]][tour[(i + 1) % n]];
    }
    return totalLen;
}

// Измерение времени выполнения
template<typename Func, typename... Args>
static auto timeIt(Func&& fn, Args&&... args) {
    auto t0 = Clock::now();
    auto result = fn(std::forward<Args>(args)...);
    auto t1 = Clock::now();
    long long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    return std::make_pair(result, ms);
}

int main() {
    std::ofstream out("results_omp.csv");
    out << "instance,algorithm,length_ms,time_ms\n";

    // Обход всех TSP-файлов
    for (auto& entry : fs::directory_iterator("data")) {
        if (!entry.is_regular_file()) continue;
        auto pts = readTSP(entry.path());
        auto [tourLength, duration] = timeIt(solve3Opt, pts);
        out << entry.path().filename().string()
            << ",3opt_openmp," << tourLength << "," << duration << "\n";
    }

    out.close();
    std::cout << "Done: results.csv generated\n";
    return 0;
}
