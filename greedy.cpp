#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <filesystem>
#include <chrono>

namespace fs = std::filesystem;

// Represents one item in the knapsack instance
struct PackItem {
    long long value;
    long long weight;
};

// Holds the outcome of running a greedy heuristic on one file
struct PackOutcome {
    std::string filename;
    std::string method; // "ratio" or "value"
    long long totalValue = 0;
    long long totalWeight = 0;
    long long timeMillis = 0;
    std::vector<long long> chosenValues;
    std::vector<long long> chosenWeights;
};

// Greedy knapsack: either sort by (value/weight) descending, or by value descending
PackOutcome solveGreedy(
    std::vector<PackItem> items,
    long long capacity,
    const std::string &fileId,
    const std::string &methodId
) {
    auto t_start = std::chrono::high_resolution_clock::now();

    if (methodId == "ratio") {
        std::sort(items.begin(), items.end(),
            [](const PackItem &a, const PackItem &b) {
                long double ra = (long double)a.value / a.weight;
                long double rb = (long double)b.value / b.weight;
                if (ra == rb) return a.value > b.value;
                return ra > rb;
            }
        );
    }
    else { // methodId == "value"
        std::sort(items.begin(), items.end(),
            [](const PackItem &a, const PackItem &b) {
                return a.value > b.value;
            }
        );
    }

    PackOutcome result;
    result.filename = fileId;
    result.method = methodId;

    long long currWeight = 0;
    for (const auto &it : items) {
        if (currWeight + it.weight <= capacity) {
            currWeight += it.weight;
            result.totalValue += it.value;
            result.totalWeight += it.weight;
            result.chosenValues.push_back(it.value);
            result.chosenWeights.push_back(it.weight);
        }
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    result.timeMillis = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();
    return result;
}

int main(int argc, char *argv[]) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    // Determine directory containing input files
    fs::path dirPath = "data";
    if (argc >= 2) {
        dirPath = argv[1];
    }
    if (!fs::exists(dirPath) || !fs::is_directory(dirPath)) {
        std::cerr << "Cannot access directory: " << dirPath << "\n";
        return 1;
    }

    std::vector<PackOutcome> allResults;

    // Iterate over every file in that directory
    for (auto &entry : fs::directory_iterator(dirPath)) {
        if (!entry.is_regular_file()) continue;

        std::ifstream inFile(entry.path());
        if (!inFile.is_open()) {
            std::cerr << "Could not open: " << entry.path() << "\n";
            continue;
        }

        size_t N;
        long long W;
        inFile >> N >> W;
        std::vector<PackItem> items(N);
        for (size_t i = 0; i < N; ++i) {
            inFile >> items[i].value >> items[i].weight;
        }
        inFile.close();

        std::string baseName = entry.path().filename().string();
        allResults.push_back(solveGreedy(items, W, baseName, "ratio"));
        allResults.push_back(solveGreedy(items, W, baseName, "value"));
    }

    // 1) Console-friendly output
    for (const auto &r : allResults) {
        std::cout << "Filename: " << r.filename
                  << " Approach: " << r.method << "\n";
        std::cout << " Total Value: " << r.totalValue
                  << " Total Weight: " << r.totalWeight
                  << " Time Taken: " << r.timeMillis << " ms\n";

        std::cout << " Values Chosen: ";
        for (size_t i = 0; i < r.chosenValues.size(); ++i) {
            if (i) std::cout << " ";
            std::cout << r.chosenValues[i];
        }
        std::cout << "\n";

        std::cout << " Weights Chosen: ";
        for (size_t i = 0; i < r.chosenWeights.size(); ++i) {
            if (i) std::cout << " ";
            std::cout << r.chosenWeights[i];
        }
        std::cout << "\n\n";
    }

    // 2) Write same info into results.csv
    std::ofstream outCsv("results_greedy_knapsack.csv");
    if (!outCsv.is_open()) {
        std::cerr << "Failed to write to results.csv\n";
        return 1;
    }

    for (const auto &r : allResults) {
        outCsv << "Filename: " << r.filename
               << " Approach: " << r.method << "\n";
        outCsv << " Total Value: " << r.totalValue
               << " Total Weight: " << r.totalWeight
               << " Time Taken: " << r.timeMillis << " ms\n";

        outCsv << " Values Chosen: ";
        for (size_t i = 0; i < r.chosenValues.size(); ++i) {
            if (i) outCsv << " ";
            outCsv << r.chosenValues[i];
        }
        outCsv << "\n";

        outCsv << " Weights Chosen: ";
        for (size_t i = 0; i < r.chosenWeights.size(); ++i) {
            if (i) outCsv << " ";
            outCsv << r.chosenWeights[i];
        }
        outCsv << "\n\n";
    }

    std::cerr << "Output written to results.csv\n";
    return 0;
}
