#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <chrono>
#include <numeric>
#include <iomanip>
#include <omp.h>



struct Point {
    double x, y;
    Point(double x = 0, double y = 0) : x(x), y(y) {}
};


// ��������� �������� ������
std::vector<Point> generateRandomCities(int numCities, double minX = 0, double maxX = 100, double minY = 0, double maxY = 100) {
    std::vector<Point> cities;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> xDist(minX, maxX);
    std::uniform_real_distribution<double> yDist(minY, maxY);

    for (int i = 0; i < numCities; i++) {
        cities.emplace_back(xDist(gen), yDist(gen));
    }

    return cities;
}


double distance(const Point& a, const Point& b) {
    return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

class TSPSolver {
private:
    std::vector<Point> cities;
    std::vector<std::vector<double>> distMatrix;
    std::mt19937 rng;

    void precomputeDistances() {
        int n = cities.size();
        distMatrix.resize(n, std::vector<double>(n));

        #pragma omp parallel for collapse(2)
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    distMatrix[i][j] = 0;
                } else {
                    distMatrix[i][j] = distance(cities[i], cities[j]);
                }
            }
        }
    }

    double calculateRouteDistance(const std::vector<int>& route) const {
        double totalDistance = 0;
        int n = route.size();

        for (int i = 0; i < n; i++) {
            int from = route[i];
            int to = route[(i + 1) % n]; // ��������� �������
            totalDistance += distMatrix[from][to];
        }

        return totalDistance;
    }

    std::vector<int> nearestNeighborHeuristic() const {
        int n = cities.size();
        std::vector<int> route;
        std::vector<bool> visited(n, false);

        int current = 0;
        route.push_back(current);
        visited[current] = true;

        for (int i = 1; i < n; i++) {
            int nearest = -1;
            double minDist = std::numeric_limits<double>::max();

            for (int j = 0; j < n; j++) {
                if (!visited[j] && distMatrix[current][j] < minDist) {
                    minDist = distMatrix[current][j];
                    nearest = j;
                }
            }

            route.push_back(nearest);
            visited[nearest] = true;
            current = nearest;
        }

        return route;
    }

    void swapSections(std::vector<int>& route, int start1, int end1, int start2, int end2) {
        std::vector<int> section1(route.begin() + start1, route.begin() + end1 + 1);
        std::vector<int> section2(route.begin() + start2, route.begin() + end2 + 1);

        int idx = start1;
        for (int city : section2) {
            route[idx++] = city;
        }
        for (int i = end1 + 1; i < start2; i++) {
            route[idx++] = route[i];
        }
        for (int city : section1) {
            route[idx++] = city;
        }
    }

    void swapSectionsReverse(std::vector<int>& route, int start1, int end1, int start2, int end2) {
        std::vector<int> section1(route.begin() + start1, route.begin() + end1 + 1);
        std::vector<int> section2(route.begin() + start2, route.begin() + end2 + 1);
        std::reverse(section2.begin(), section2.end());

        int idx = start1;
        for (int city : section2) {
            route[idx++] = city;
        }
        for (int i = end1 + 1; i < start2; i++) {
            route[idx++] = route[i];
        }
        for (int city : section1) {
            route[idx++] = city;
        }
    }

public:
    TSPSolver(const std::vector<Point>& cities) : cities(cities), rng(std::chrono::steady_clock::now().time_since_epoch().count()) {
        precomputeDistances();
    }

    std::vector<int> solve(int maxIterations = 10000, double initialTemp = 1000.0, double coolingRate = 0.995) {
        int n = cities.size();
        if (n < 2) return {};

        std::vector<int> bestRoute = nearestNeighborHeuristic();
        double bestDistance = calculateRouteDistance(bestRoute);

        std::cout << "nearest neighbor: " << std::fixed << std::setprecision(2) << bestDistance << std::endl;

        const int numThreads = omp_get_max_threads();
        std::vector<std::vector<int>> threadBestRoutes(numThreads, bestRoute);
        std::vector<double> threadBestDistances(numThreads, bestDistance);

        #pragma omp parallel
        {
            int threadId = omp_get_thread_num();
            std::mt19937 localRng(rng() + threadId);
            std::uniform_real_distribution<double> uniform(0.0, 1.0);
            std::uniform_int_distribution<int> cityDist(0, n - 1);

            std::vector<int> currentRoute = bestRoute;
            double currentDistance = bestDistance;
            double temperature = initialTemp;

            for (int iteration = 0; iteration < maxIterations; iteration++) {
                std::vector<int> newRoute = currentRoute;

                // ��������� ����� ���� �������
                int mutationType = localRng() % 5;

                switch (mutationType) {
                    case 0: { // 2-opt
                        int i = cityDist(localRng);
                        int j = cityDist(localRng);
                        if (i > j) std::swap(i, j);
                        if (i != j) {
                            std::reverse(newRoute.begin() + i, newRoute.begin() + j + 1);
                        }
                        break;
                    }
                    case 1: { // Swap ���� �������
                        int i = cityDist(localRng);
                        int j = cityDist(localRng);
                        if (i != j) {
                            std::swap(newRoute[i], newRoute[j]);
                        }
                        break;
                    }
                    case 2: { // ������� ������ � ������ �������
                        int from = cityDist(localRng);
                        int to = cityDist(localRng);
                        if (from != to) {
                            int city = newRoute[from];
                            newRoute.erase(newRoute.begin() + from);
                            newRoute.insert(newRoute.begin() + to, city);
                        }
                        break;
                    }
                    case 3: { // Or-opt (����������� ��������)
                        int segmentSize = (localRng() % 3) + 1;
                        int start = cityDist(localRng);
                        if (start + segmentSize <= n) {
                            int insertPos = cityDist(localRng);
                            if (insertPos < start || insertPos > start + segmentSize) {
                                std::vector<int> segment(newRoute.begin() + start, newRoute.begin() + start + segmentSize);
                                newRoute.erase(newRoute.begin() + start, newRoute.begin() + start + segmentSize);
                                if (insertPos > start) insertPos -= segmentSize;
                                newRoute.insert(newRoute.begin() + insertPos, segment.begin(), segment.end());
                            }
                        }
                        break;
                    }
                    case 4: { // 3-opt ����������
                        int i = cityDist(localRng);
                        int j = cityDist(localRng);
                        int k = cityDist(localRng);

                        std::vector<int> indices = {i, j, k};
                        std::sort(indices.begin(), indices.end());

                        if (indices[0] != indices[1] && indices[1] != indices[2]) {
                            if (localRng() % 2) {
                                std::reverse(newRoute.begin() + indices[0], newRoute.begin() + indices[1] + 1);
                            } else {
                                std::reverse(newRoute.begin() + indices[1], newRoute.begin() + indices[2] + 1);
                            }
                        }
                        break;
                    }
                }

                double newDistance = calculateRouteDistance(newRoute);

                // �������� ��� ���������� ������� (Simulated Annealing)
                double deltaE = newDistance - currentDistance;
                if (deltaE < 0 || (temperature > 0 && uniform(localRng) < exp(-deltaE / temperature))) {
                    currentRoute = newRoute;
                    currentDistance = newDistance;

                    // ���������� ������� ������� ��� ������
                    if (currentDistance < threadBestDistances[threadId]) {
                        threadBestRoutes[threadId] = currentRoute;
                        threadBestDistances[threadId] = currentDistance;
                    }
                }

                temperature *= coolingRate;

                // ������������� ����� ��������� ����� ��������
                if (iteration % 1000 == 0) {
                    #pragma omp critical
                    {
                        for (int t = 0; t < numThreads; t++) {
                            if (threadBestDistances[t] < bestDistance) {
                                bestRoute = threadBestRoutes[t];
                                bestDistance = threadBestDistances[t];
                            }
                        }

                        // ��������� ��� ������ ������ ��������� ��������
                        for (int t = 0; t < numThreads; t++) {
                            if (threadBestDistances[t] > bestDistance * 1.1) { // ���� ������� ������ ������� ����
                                threadBestRoutes[t] = bestRoute;
                                threadBestDistances[t] = bestDistance;
                                if (t == threadId) {
                                    currentRoute = bestRoute;
                                    currentDistance = bestDistance;
                                }
                            }
                        }
                    }
                }
            }
        }

        // ��������� ����� ������� �������
        for (int t = 0; t < numThreads; t++) {
            if (threadBestDistances[t] < bestDistance) {
                bestRoute = threadBestRoutes[t];
                bestDistance = threadBestDistances[t];
            }
        }

        std::cout << "Result: " << std::fixed << std::setprecision(2) << bestDistance << std::endl;
        std::cout << "Using thread: " << numThreads << std::endl;

        return bestRoute;
    }
};

int main() {
    // ������������� ���������� �������
    omp_set_num_threads(omp_get_max_threads());

    std::cout << "Threads available: " << omp_get_max_threads() << std::endl << std::endl;

    // ��������� �������� �������
    int numCities = 50;
    std::vector<Point> cities = generateRandomCities(numCities);

    auto start = std::chrono::high_resolution_clock::now();

    TSPSolver solver(cities);
    std::vector<int> solution = solver.solve(20000, 2000.0, 0.9995);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Execution time: " << duration.count() << " ��" << std::endl;

    // ����� ��������
    std::cout << "\n Optimal: ";
    for (size_t i = 0; i < solution.size(); i++) {
        std::cout << solution[i];
        if (i < solution.size() - 1) std::cout << " -> ";
    }
    std::cout << " -> " << solution[0] << std::endl;  // ������� � ������

    return 0;
}
