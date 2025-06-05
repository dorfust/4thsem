#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

struct Point {
    double x, y;
};

double euclideanDistance(const Point& a, const Point& b) {
    return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

double calculateTotalDistance(const vector<Point>& points, const vector<int>& path) {
    double totalDistance = 0.0;
    for (size_t i = 0; i < path.size() - 1; ++i) {
        totalDistance += euclideanDistance(points[path[i]], points[path[i + 1]]);
    }
    totalDistance += euclideanDistance(points[path.back()], points[path.front()]); // Closing the cycle
    return totalDistance;
}

void twoOptSwap(vector<int>& path, int i, int k) {
    reverse(path.begin() + i, path.begin() + k + 1);
}

bool twoOpt(vector<Point>& points, vector<int>& path) {
    int n = path.size();
    bool improved = false;
    double bestDistance = calculateTotalDistance(points, path);

    for (int i = 1; i < n - 1; ++i) {
        for (int k = i + 1; k < n; ++k) {
            twoOptSwap(path, i, k);
            double newDistance = calculateTotalDistance(points, path);
            if (newDistance < bestDistance) {
                bestDistance = newDistance;
                improved = true;
            } else {
                twoOptSwap(path, i, k); // Revert the swap
            }
        }
    }
    return improved;
}

void solveTSP(const string& filename, std::ofstream& outFile) {
    ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    int n;
    inputFile >> n;
    vector<Point> points(n);
    for (int i = 0; i < n; ++i) {
        inputFile >> points[i].x >> points[i].y;
    }
    inputFile.close();

    vector<int> path(n);
    for (int i = 0; i < n; ++i) {
        path[i] = i;
    }

    bool improved = true;
    while (improved) {
        improved = twoOpt(points, path);
    }


    //std::ofstream outFile("output.txt");


    double totalDistance = calculateTotalDistance(points, path);
    cout << totalDistance << " 0" << endl;
    outFile << totalDistance << " 0" <<endl;
    for (int i = 0; i < n; ++i) {
        cout << path[i] << " ";
        outFile << path[i] << " ";
    }
    cout << endl;
    outFile << std::endl;
    //outFile.close();
}

int main() {
    std::ofstream outFile("output.txt");
    string dataPath = "data";
    for (const auto& entry : fs::directory_iterator(dataPath)) {
        if (entry.is_regular_file()) {
            solveTSP(entry.path().string(), outFile);
        }
    }
    outFile.close();
    return 0;
}
