#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <algorithm>
#include <queue>
#include <filesystem>

namespace fs = std::filesystem;

// -------------------- Data Structures -------------------- //
struct Obj {
    int weight;
    int value;
};

struct NodeBB {
    int index;                   // current level in sorted list
    int currentWeight;
    int currentValue;
    double upperBound;
    std::vector<int> taken;      // 0/1 flags in sorted order
};

struct ResultBB {
    int bestValue;
    std::vector<int> bestTaken;  // 0/1 flags in sorted order
};

struct Individual {
    std::vector<int> genes; // 0/1 flags in original order
    int fitness;
};

struct ResultGA {
    int bestValue;
    std::vector<int> bestGenes; // 0/1 flags in original order
};

// -------------------- Time Utility -------------------- //
long long nowMillis() {
    using namespace std::chrono;
    return duration_cast<milliseconds>(steady_clock::now().time_since_epoch()).count();
}

// -------------------- GA Repair for Overweight -------------------- //
void repairOverweight(Individual &ind, const std::vector<Obj> &objects, int capacity) {
    int totalW = 0;
    for (size_t i = 0; i < ind.genes.size(); ++i) {
        if (ind.genes[i] == 1) {
            totalW += objects[i].weight;
        }
    }
    if (totalW <= capacity) return;

    // Collect indices of chosen items
    std::vector<int> chosen;
    for (size_t i = 0; i < ind.genes.size(); ++i) {
        if (ind.genes[i] == 1) {
            chosen.push_back((int)i);
        }
    }
    // Sort by increasing value/weight ratio
    std::sort(chosen.begin(), chosen.end(),
              [&](int a, int b) {
                  double ra = (double)objects[a].value / objects[a].weight;
                  double rb = (double)objects[b].value / objects[b].weight;
                  return ra < rb;
              });
    // Remove lowest-ratio items until under capacity
    for (int idx : chosen) {
        if (totalW <= capacity) break;
        ind.genes[idx] = 0;
        totalW -= objects[idx].weight;
    }
}

// -------------------- GA: Fitness and Operators -------------------- //
int computeFitness(const Individual &ind, const std::vector<Obj> &objects, int capacity) {
    int sumW = 0, sumV = 0;
    for (size_t i = 0; i < ind.genes.size(); ++i) {
        if (ind.genes[i] == 1) {
            sumW += objects[i].weight;
            sumV += objects[i].value;
        }
    }
    if (sumW > capacity) return 0;
    return sumV;
}

Individual makeRandomIndividual(int n) {
    Individual ind;
    ind.genes.resize(n);
    for (int i = 0; i < n; ++i) {
        ind.genes[i] = rand() % 2;
    }
    ind.fitness = 0;
    return ind;
}

void crossover(Individual &a, Individual &b) {
    int L = (int)a.genes.size();
    int cut = rand() % L;
    for (int i = cut; i < L; ++i) {
        std::swap(a.genes[i], b.genes[i]);
    }
}

void mutate(Individual &ind, double mutationRate) {
    for (size_t i = 0; i < ind.genes.size(); ++i) {
        if ((double)rand() / RAND_MAX < mutationRate) {
            ind.genes[i] = 1 - ind.genes[i];
        }
    }
}

Individual tournamentSelection(const std::vector<Individual> &population) {
    int sz = (int)population.size();
    int i1 = rand() % sz;
    int i2 = rand() % sz;
    return (population[i1].fitness > population[i2].fitness) ? population[i1] : population[i2];
}

ResultGA runGeneticAlgorithm(const std::vector<Obj> &objects, int capacity) {
    const int POP_SIZE = 60;
    const int MAX_GEN = 250;
    const double MUTATION_RATE = 0.04;

    int n = (int)objects.size();
    std::vector<Individual> population(POP_SIZE);

    // Initialize population
    for (int i = 0; i < POP_SIZE; ++i) {
        population[i] = makeRandomIndividual(n);
        repairOverweight(population[i], objects, capacity);
        population[i].fitness = computeFitness(population[i], objects, capacity);
    }

    int bestFitness = 0;
    std::vector<int> bestGenes(n, 0);

    // Evolve
    for (int gen = 0; gen < MAX_GEN; ++gen) {
        std::vector<Individual> newPop;

        // Elitism: keep the best in this generation
        Individual champion = population[0];
        for (auto &ind : population) {
            if (ind.fitness > champion.fitness) {
                champion = ind;
            }
        }
        if (champion.fitness > bestFitness) {
            bestFitness = champion.fitness;
            bestGenes = champion.genes;
        }
        newPop.push_back(champion);

        // Fill new population
        while ((int)newPop.size() < POP_SIZE) {
            Individual p1 = tournamentSelection(population);
            Individual p2 = tournamentSelection(population);
            crossover(p1, p2);
            mutate(p1, MUTATION_RATE);
            mutate(p2, MUTATION_RATE);
            repairOverweight(p1, objects, capacity);
            repairOverweight(p2, objects, capacity);
            p1.fitness = computeFitness(p1, objects, capacity);
            p2.fitness = computeFitness(p2, objects, capacity);
            newPop.push_back(p1);
            if ((int)newPop.size() < POP_SIZE) {
                newPop.push_back(p2);
            }
        }
        population = std::move(newPop);
    }

    // Final check
    for (auto &ind : population) {
        if (ind.fitness > bestFitness) {
            bestFitness = ind.fitness;
            bestGenes = ind.genes;
        }
    }

    return { bestFitness, bestGenes };
}

// -------------------- BnB: Compute Upper Bound -------------------- //
double computeUpperBound(const NodeBB &node, int n, int capacity, const std::vector<Obj> &sortedObjs) {
    if (node.currentWeight >= capacity) return 0.0;
    double bound = (double)node.currentValue;
    int wSoFar = node.currentWeight;
    int idx = node.index;

    while (idx < n && wSoFar + sortedObjs[idx].weight <= capacity) {
        wSoFar += sortedObjs[idx].weight;
        bound += sortedObjs[idx].value;
        ++idx;
    }
    if (idx < n) {
        bound += (capacity - wSoFar) * (double)sortedObjs[idx].value / sortedObjs[idx].weight * 0.95;
    }
    return bound;
}

struct CompareBB {
    bool operator()(const NodeBB &a, const NodeBB &b) {
        return a.upperBound < b.upperBound;
    }
};

ResultBB runBranchAndBound(const std::vector<Obj> &origObjs, int capacity) {
    const long long TIME_LIMIT_MS = 120000; // 120 seconds
    const int MAX_LEVEL = 80;

    long long startTime = nowMillis();
    int n = (int)origObjs.size();

    // Sort objects by descending value/weight
    std::vector<Obj> sortedObjs = origObjs;
    std::sort(sortedObjs.begin(), sortedObjs.end(),
              [&](const Obj &a, const Obj &b) {
                  return (double)a.value / a.weight > (double)b.value / b.weight;
              });

    // Priority queue (max-heap by upper bound)
    std::priority_queue<NodeBB, std::vector<NodeBB>, CompareBB> pq;
    NodeBB root;
    root.index = 0;
    root.currentWeight = 0;
    root.currentValue = 0;
    root.taken.assign(n, 0);
    root.upperBound = computeUpperBound(root, n, capacity, sortedObjs);
    pq.push(root);

    int bestValue = 0;
    std::vector<int> bestTaken(n, 0);

    while (!pq.empty()) {
        if (nowMillis() - startTime > TIME_LIMIT_MS) break;
        NodeBB node = pq.top();
        pq.pop();

        if (node.index >= MAX_LEVEL || node.index >= n) continue;
        if (node.upperBound <= (double)bestValue) continue;

        // Option 1: Take sortedObjs[node.index]
        NodeBB withItem = node;
        withItem.index = node.index + 1;
        withItem.currentWeight = node.currentWeight + sortedObjs[node.index].weight;
        withItem.currentValue = node.currentValue + sortedObjs[node.index].value;
        withItem.taken = node.taken;
        withItem.taken[node.index] = 1;
        if (withItem.currentWeight <= capacity && withItem.currentValue > bestValue) {
            bestValue = withItem.currentValue;
            bestTaken = withItem.taken;
        }
        withItem.upperBound = computeUpperBound(withItem, n, capacity, sortedObjs);
        if (withItem.upperBound > (double)bestValue) {
            pq.push(withItem);
        }

        // Option 2: Skip sortedObjs[node.index]
        NodeBB skipItem = node;
        skipItem.index = node.index + 1;
        skipItem.taken = node.taken;
        skipItem.taken[node.index] = 0;
        skipItem.upperBound = computeUpperBound(skipItem, n, capacity, sortedObjs);
        if (skipItem.upperBound > (double)bestValue) {
            pq.push(skipItem);
        }
    }

    return { bestValue, bestTaken };
}

// -------------------- MAIN: Process All Files -------------------- //
int main() {
    const std::string dataFolder = "data";
    const std::string summaryCsv = "summary_knapsack_results.csv";

    std::cout << "Current working directory: " << fs::current_path() << std::endl;
    std::cout << "Scanning folder: " << dataFolder << std::endl;

    if (!fs::exists(dataFolder) || !fs::is_directory(dataFolder)) {
        std::cerr << "Error: cannot open directory \"" << dataFolder << "\"" << std::endl;
        return 1;
    }

    std::ofstream summaryOut(summaryCsv, std::ios::out | std::ios::trunc);
    if (!summaryOut.is_open()) {
        std::cerr << "Error: cannot create file \"" << summaryCsv << "\"" << std::endl;
        return 1;
    }
    summaryOut << "File,BestValue_BnB,Weight_BnB,Time_BnB_ms,BestValue_GA,Weight_GA,Time_GA_ms\n";

    // Create detail-results directory if it doesn't exist
    if (!fs::exists("knapsack_individual_results")) {
        fs::create_directory("knapsack_individual_results");
    }

    for (const auto &entry : fs::directory_iterator(dataFolder)) {
        if (!entry.is_regular_file()) continue;
        std::string inputPath = entry.path().string();
        std::cout << "Processing file: " << entry.path().filename().string() << std::endl;

        std::ifstream fin(inputPath);
        if (!fin.is_open()) {
            std::cerr << "  Cannot open file: " << inputPath << std::endl;
            continue;
        }

        int N, CAP;
        fin >> N >> CAP;
        if (N <= 0) {
            fin.close();
            continue;
        }
        std::vector<Obj> objects(N);
        for (int i = 0; i < N; ++i) {
            fin >> objects[i].weight >> objects[i].value;
        }
        fin.close();

        // ----- Branch & Bound ----- //
        long long t0 = nowMillis();
        ResultBB resBB = runBranchAndBound(objects, CAP);
        long long t1 = nowMillis();
        int valueBB = resBB.bestValue;
        int weightBB = 0;
        {
            // Recreate the same sorted order
            std::vector<Obj> sortedObjs = objects;
            std::sort(sortedObjs.begin(), sortedObjs.end(),
                      [&](const Obj &a, const Obj &b) {
                          return (double)a.value / a.weight > (double)b.value / b.weight;
                      });
            for (size_t i = 0; i < resBB.bestTaken.size(); ++i) {
                if (resBB.bestTaken[i] == 1) {
                    weightBB += sortedObjs[i].weight;
                }
            }
        }
        long long timeBB = t1 - t0;
        std::cout << "  [BnB] Value=" << valueBB << ", Weight=" << weightBB
                  << ", Time=" << timeBB << " ms" << std::endl;

        // ----- Genetic Algorithm ----- //
        long long t2 = nowMillis();
        ResultGA resGA = runGeneticAlgorithm(objects, CAP);
        long long t3 = nowMillis();
        int valueGA = resGA.bestValue;
        int weightGA = 0;
        for (int i = 0; i < N; ++i) {
            if (resGA.bestGenes[i] == 1) {
                weightGA += objects[i].weight;
            }
        }
        long long timeGA = t3 - t2;
        std::cout << "  [GA ] Value=" << valueGA << ", Weight=" << weightGA
                  << ", Time=" << timeGA << " ms" << std::endl;

        // ----- Append to Summary CSV ----- //
        summaryOut << entry.path().filename().string() << ","
                   << valueBB << "," << weightBB << "," << timeBB << ","
                   << valueGA << "," << weightGA << "," << timeGA << "\n";

        // ----- Write Detailed CSV per File ----- //
        std::string baseName = entry.path().stem().string();
        std::string detailPath = "knapsack_individual_results/" + baseName + ".csv";
        std::ofstream fout(detailPath, std::ios::out | std::ios::trunc);
        if (fout.is_open()) {
            // Write BnB selection weights
            fout << "BnB,";
            {
                std::vector<Obj> sortedObjs = objects;
                std::sort(sortedObjs.begin(), sortedObjs.end(),
                          [&](const Obj &a, const Obj &b) {
                              return (double)a.value / a.weight > (double)b.value / b.weight;
                          });
                for (size_t i = 0; i < resBB.bestTaken.size(); ++i) {
                    if (resBB.bestTaken[i] == 1) {
                        fout << sortedObjs[i].weight;
                    } else {
                        fout << 0;
                    }
                    if (i + 1 < resBB.bestTaken.size()) fout << ",";
                }
            }
            fout << "\n";

            // Write GA selection weights
            fout << "GA,";
            for (int i = 0; i < N; ++i) {
                if (resGA.bestGenes[i] == 1) {
                    fout << objects[i].weight;
                } else {
                    fout << 0;
                }
                if (i + 1 < N) fout << ",";
            }
            fout << "\n";
            fout.close();
        }
    }

    summaryOut.close();
    std::cout << "Finished. Summary saved to \"" << summaryCsv << "\"" << std::endl;
    return 0;
}
