#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <cmath>
#include <algorithm>

namespace fs = std::filesystem;

// Structure to hold a point's coordinates
struct Point {
    double x, y;
};

// Compute Euclidean distance between two points
inline double computeDistance(const Point &p1, const Point &p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return std::sqrt(dx * dx + dy * dy);
}

// Calculate the total length of a tour (closed or open)
double calculateTourLength(const std::vector<int> &tour,
                           const std::vector<Point> &points,
                           bool closeCycle)
{
    double total = 0.0;
    int sz = static_cast<int>(tour.size());
    for (int i = 0; i + 1 < sz; ++i) {
        total += computeDistance(points[tour[i]], points[tour[i + 1]]);
    }
    if (closeCycle && sz > 1) {
        total += computeDistance(points[tour.back()], points[tour.front()]);
    }
    return total;
}

// One 2-opt iteration: try all (i, j) pairs with reversal
bool improveTwoOpt(std::vector<int> &tour,
                   const std::vector<Point> &pts,
                   bool closeCycle)
{
    int n = static_cast<int>(tour.size());
    double bestLen = calculateTourLength(tour, pts, closeCycle);
    bool foundBetter = false;

    for (int i = 1; i < n - 2; ++i) {
        for (int j = i + 1; j < n - (closeCycle ? 0 : 1); ++j) {
            std::reverse(tour.begin() + i, tour.begin() + j);
            double newLen = calculateTourLength(tour, pts, closeCycle);
            if (newLen < bestLen) {
                bestLen = newLen;
                foundBetter = true;
            } else {
                std::reverse(tour.begin() + i, tour.begin() + j);  // revert
            }
        }
    }
    return foundBetter;
}

// Run 2-opt until no improvement; return (length, elapsed_ms)
std::pair<double, long long> performTwoOpt(std::vector<int> &tour,
                                           const std::vector<Point> &pts,
                                           bool closeCycle)
{
    auto t0 = std::chrono::high_resolution_clock::now();
    while (improveTwoOpt(tour, pts, closeCycle)) { }
    auto t1 = std::chrono::high_resolution_clock::now();
    long long elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    double length = calculateTourLength(tour, pts, closeCycle);
    return { length, elapsed };
}

// One 3-opt iteration: try (i, j, k) triples with several reversals
bool improveThreeOpt(std::vector<int> &tour,
                     const std::vector<Point> &pts,
                     bool closeCycle)
{
    int n = static_cast<int>(tour.size());
    double bestLen = calculateTourLength(tour, pts, closeCycle);

    for (int i = 0; i < n - 2; ++i) {
        for (int j = i + 1; j < n - 1; ++j) {
            for (int k = j + 1; k < n; ++k) {
                std::vector<int> candidate = tour;
                double localBest = bestLen;
                std::vector<int> bestCandidate = tour;

                auto tryModification = [&](int mode) {
                    candidate = tour;
                    if (mode == 0) {
                        std::reverse(candidate.begin() + (i + 1), candidate.begin() + (j + 1));
                    } else if (mode == 1) {
                        std::reverse(candidate.begin() + (j + 1), candidate.begin() + (k + 1));
                    } else if (mode == 2) {
                        std::reverse(candidate.begin() + (i + 1), candidate.begin() + (j + 1));
                        std::reverse(candidate.begin() + (j + 1), candidate.begin() + (k + 1));
                    } else if (mode == 3) {
                        std::reverse(candidate.begin() + (i + 1), candidate.begin() + (k + 1));
                    } else if (mode == 4) {
                        std::reverse(candidate.begin() + (i + 1), candidate.begin() + (j + 1));
                        std::reverse(candidate.begin() + (i + 1), candidate.begin() + (k + 1));
                    } else if (mode == 5) {
                        std::reverse(candidate.begin() + (j + 1), candidate.begin() + (k + 1));
                        std::reverse(candidate.begin() + (i + 1), candidate.begin() + (k + 1));
                    } else if (mode == 6) {
                        std::reverse(candidate.begin() + (i + 1), candidate.begin() + (j + 1));
                        std::reverse(candidate.begin() + (j + 1), candidate.begin() + (k + 1));
                        std::reverse(candidate.begin() + (i + 1), candidate.begin() + (k + 1));
                    }

                    double len = calculateTourLength(candidate, pts, closeCycle);
                    if (len < localBest) {
                        localBest = len;
                        bestCandidate = candidate;
                    }
                };

                for (int mode = 0; mode < 7; ++mode) {
                    tryModification(mode);
                }

                if (localBest < bestLen) {
                    tour = bestCandidate;
                    return true;
                }
            }
        }
    }
    return false;
}

// Run 3-opt until no improvement; return (length, elapsed_ms)
std::pair<double, long long> performThreeOpt(std::vector<int> &tour,
                                             const std::vector<Point> &pts,
                                             bool closeCycle)
{
    auto t0 = std::chrono::high_resolution_clock::now();
    while (improveThreeOpt(tour, pts, closeCycle)) { }
    auto t1 = std::chrono::high_resolution_clock::now();
    long long elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    double length = calculateTourLength(tour, pts, closeCycle);
    return { length, elapsed };
}

int main() {
    const std::string dataDir = "data";
    const std::string outputCsv = "results_2opt_3opt.csv";

    // Print the current working directory and which folder we try to open
    std::cout << "Current working directory: " << fs::current_path() << std::endl;
    std::cout << "Attempting to open directory: " << dataDir << std::endl;

    if (!fs::exists(dataDir) || !fs::is_directory(dataDir)) {
        std::cerr << "Error: cannot find or open folder \"" << dataDir << "\"" << std::endl;
        return 1;
    }

    std::ofstream csvFile(outputCsv, std::ios::out | std::ios::trunc);
    if (!csvFile.is_open()) {
        std::cerr << "Error: cannot open \"" << outputCsv << "\" for writing" << std::endl;
        return 1;
    }

    // Write CSV header
    csvFile << "filename,dist_2opt,time_2opt_ms,dist_3opt,time_3opt_ms\n";

    std::cout << "Starting to iterate through files in \"" << dataDir << "\"" << std::endl;
    for (const auto &entry : fs::directory_iterator(dataDir)) {
        if (!entry.is_regular_file()) {
            std::cout << "  -> Skipping " << entry.path().filename().string()
                      << " (not a regular file)" << std::endl;
            continue;
        }

        std::string filename = entry.path().string();
        std::cout << "Found file: " << entry.path().filename().string() << std::endl;

        std::ifstream inFile(filename);
        if (!inFile.is_open()) {
            std::cerr << "  Error: cannot open file: " << filename << std::endl;
            continue;
        }

        int nPoints;
        inFile >> nPoints;
        std::cout << "  Read number of points = " << nPoints << std::endl;
        if (nPoints <= 0) {
            std::cout << "  -> Skipping file, invalid number of points" << std::endl;
            inFile.close();
            continue;
        }

        std::vector<Point> allPoints(nPoints);
        for (int i = 0; i < nPoints; ++i) {
            double xx, yy;
            inFile >> xx >> yy;
            allPoints[i] = {xx, yy};
        }
        inFile.close();

        // Initialize the tour as 0,1,2,...,nPoints-1
        std::vector<int> tour(nPoints);
        for (int i = 0; i < nPoints; ++i) {
            tour[i] = i;
        }

        bool isClosed = true; // Treat as a closed tour

        // Run 2-opt
        auto [length2, time2] = performTwoOpt(tour, allPoints, isClosed);
        std::cout << "  2-opt result: length=" << length2 << ", time=" << time2 << " ms" << std::endl;

        // Run 3-opt on the improved tour
        auto [length3, time3] = performThreeOpt(tour, allPoints, isClosed);
        std::cout << "  3-opt result: length=" << length3 << ", time=" << time3 << " ms" << std::endl;

        // Write results to CSV
        csvFile << fs::path(filename).filename().string() << ","
                << std::fixed << std::setprecision(6) << length2 << ","
                << time2 << ","
                << length3 << ","
                << time3 << "\n";
    }

    csvFile.close();
    std::cout << "Done. Results have been saved to \"" << outputCsv << "\"" << std::endl;
    return 0;
}
