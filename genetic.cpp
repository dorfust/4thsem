#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <ctime>
#include <random>
#include <chrono>
#include <filesystem>
#include <string>
#include <iomanip>

using namespace std;
using namespace std::chrono;


// ��������� ��� �������� �����������
struct TestResult {
    string filename;
    double gaTime;
    int gaResult;
};



// ������� ��� ��������� ������ ������ � ����������
vector<string> getFilesInDirectory(const string& directoryPath) {
    vector<string> files;

    try {
        for (const auto& entry : filesystem::directory_iterator(directoryPath)) {
            if (entry.is_regular_file()) {
                files.push_back(entry.path().string());
            }
        }
        sort(files.begin(), files.end());
    } catch (const filesystem::filesystem_error& ex) {
        cerr << "Reading error: " << ex.what() << endl;
    }

    return files;
}



// ������� ��� ���������� ����������� � CSV ����
void saveResultsToCSV(const vector<TestResult>& results, const string& filename) {
    ofstream file(filename);
    if (!file) {
        cerr << "������ ��� �������� ����� �����������!" << endl;
        return;
    }

    file << "����, ����� � ��������,GA ���������" << endl;

    for (const auto& result : results) {
        // ��������� ������ ��� ����� ��� ������� ����
        string shortFilename = filesystem::path(result.filename).filename().string();
        file << shortFilename << "," << fixed << setprecision(4) << result.gaTime
             << "," << result.gaResult << endl;
    }

    file.close();
    cout << "Results in file: " << filename << endl;
}



// ������� ��� ������ ����������� � ������� � ���� �������
void printResultsTable(const vector<TestResult>& results) {
    cout << "\n" << string(80, '=') << endl;
    cout << "Results of genetic algorithm" << endl;
    cout << string(80, '=') << endl;

    cout << left << setw(20) << "Filename"
         << setw(20) << "time (s)"
         << setw(15) << "result" << endl;
    cout << string(80, '-') << endl;

    for (const auto& result : results) {
        string shortFilename = filesystem::path(result.filename).filename().string();
        cout << left << setw(20) << shortFilename
             << setw(20) << fixed << setprecision(4) << result.gaTime
             << setw(15) << result.gaResult << endl;
    }
    cout << string(80, '=') << endl;
}


// ����� ��� ������������� ��������
class Item
{
private:
    int weight_;
    int value_;

public:
    Item() : weight_(0), value_(0) {}
    Item(int weight, int value) : weight_(weight), value_(value) {}

    int getWeight() const { return weight_; }
    int getValue() const { return value_; }

    void setWeight(int weight) { weight_ = weight; }
    void setValue(int value) { value_ = value; }
};

// ����� ��� ������������� ����� � ���������
class Individual
{
private:
    vector<int> chromosome_;
    int fitness_;

public:
    Individual() : fitness_(0) {}
    Individual(int size) : chromosome_(size, 0), fitness_(0) {}

    vector<int>& getChromosome() { return chromosome_; }
    const vector<int>& getChromosome() const { return chromosome_; }

    int getFitness() const { return fitness_; }
    void setFitness(int fitness) { fitness_ = fitness; }

    void randomizeChromosome()
    {
        for (auto& gene : chromosome_)
        {
            gene = rand() % 2;
        }
    }

    void mutate(double mutationRate)
    {
        for (auto& gene : chromosome_)
        {
            if ((double)rand() / RAND_MAX < mutationRate)
            {
                gene = 1 - gene;
            }
        }
    }

    int calculateTotalWeight(const vector<Item>& items) const
    {
        int totalWeight = 0;
        for (size_t i = 0; i < chromosome_.size(); i++)
        {
            if (chromosome_[i] == 1)
            {
                totalWeight += items[i].getWeight();
            }
        }
        return totalWeight;
    }

    int calculateTotalValue(const vector<Item>& items) const
    {
        int totalValue = 0;
        for (size_t i = 0; i < chromosome_.size(); i++)
        {
            if (chromosome_[i] == 1)
            {
                totalValue += items[i].getValue();
            }
        }
        return totalValue;
    }

    void evaluateFitness(const vector<Item>& items, int capacity)
    {
        int totalWeight = calculateTotalWeight(items);
        int totalValue = calculateTotalValue(items);

        if (totalWeight > capacity)
        {
            fitness_ = 0;
        }
        else
        {
            fitness_ = totalValue;
        }
    }

    bool operator>(const Individual& other) const
    {
        return fitness_ > other.fitness_;
    }
};

struct FitnessComparator {
    bool operator()(const Individual& a, const Individual& b) const {
        return a.getFitness() > b.getFitness();
    }
};

// ����� ��� ���������� ������������� ���������
class geneticResult
{
private:
    int bestValue_;
    vector<int> bestChromosome_;
    double executionTime_;

public:
    geneticResult() : bestValue_(0), executionTime_(0.0) {}
    geneticResult(int value, const vector<int>& chromosome, double time = 0.0)
        : bestValue_(value), bestChromosome_(chromosome), executionTime_(time) {}

    int getBestValue() const { return bestValue_; }
    const vector<int>& getBestChromosome() const { return bestChromosome_; }
    double getExecutionTime() const { return executionTime_; }

    void setBestValue(int value) { bestValue_ = value; }
    void setBestChromosome(const vector<int>& chromosome) { bestChromosome_ = chromosome; }
    void setExecutionTime(double time) { executionTime_ = time; }
};

// �������� ����� ������������� ���������
class GeneticAlgorithm
{
private:
    int populationSize_;
    int maxGenerations_;
    double mutationRate_;
    int tournamentSize_;
    vector<Individual> population_;

    Individual tournamentSelection()
    {
        vector<Individual> candidates;
        for (int i = 0; i < tournamentSize_; i++)
        {
            candidates.push_back(population_[rand() % population_.size()]);
        }

        Individual best = candidates[0];
        for (const auto& candidate : candidates)
        {
            if (candidate > best)
            {
                best = candidate;
            }
        }
        return best;
    }

    pair<Individual, Individual> crossover(Individual parent1, Individual parent2)
    {
        int size = parent1.getChromosome().size();
        int point = rand() % (size - 1);

        for (int i = point; i < size; i++)
        {
            swap(parent1.getChromosome()[i], parent2.getChromosome()[i]);
        }

        return make_pair(parent1, parent2);
    }

    void initializePopulation(int chromosomeSize)
    {
        population_.clear();
        population_.resize(populationSize_);

        for (auto& individual : population_)
        {
            individual = Individual(chromosomeSize);
            individual.randomizeChromosome();
        }
    }

    void evaluatePopulation(const vector<Item>& items, int capacity)
    {
        for (auto& individual : population_)
        {
            individual.evaluateFitness(items, capacity);
        }
    }

    void sortPopulation() {
        FitnessComparator comparator;
        std::sort(population_.begin(), population_.end(), comparator);
    }

public:
    GeneticAlgorithm(int popSize = 50, int maxGen = 200,
                     double mutRate = 0.05, int tournSize = 3)
        : populationSize_(popSize), maxGenerations_(maxGen),
          mutationRate_(mutRate), tournamentSize_(tournSize) {}

    geneticResult solve(const vector<Item>& items, int capacity)
    {
        auto start = high_resolution_clock::now();

        initializePopulation(items.size());
        evaluatePopulation(items, capacity);

        for (int generation = 0; generation < maxGenerations_; generation++)
        {
            vector<Individual> newPopulation;
            sortPopulation();

            // �������: ��������� ������ ������
            newPopulation.push_back(population_[0]);
            newPopulation.push_back(population_[1]);

            // ������� ����� ���������
            while (newPopulation.size() < populationSize_)
            {
                Individual parent1 = tournamentSelection();
                Individual parent2 = tournamentSelection();

                auto children = crossover(parent1, parent2);

                children.first.mutate(mutationRate_);
                children.second.mutate(mutationRate_);

                children.first.evaluateFitness(items, capacity);
                children.second.evaluateFitness(items, capacity);

                newPopulation.push_back(children.first);
                if (newPopulation.size() < populationSize_)
                {
                    newPopulation.push_back(children.second);
                }
            }

            population_ = newPopulation;
        }

        sortPopulation();

        auto end = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(end - start);
        double timeInSeconds = duration.count() / 1000000.0;

        geneticResult result(population_[0].getFitness(), population_[0].getChromosome(), timeInSeconds);
        return result;
    }
};

// ����� ��� ������������� ������ � �������
class KnapsackProblem
{
private:
    vector<Item> items_;
    int capacity_;

public:
    KnapsackProblem() : capacity_(0) {}
    KnapsackProblem(const vector<Item>& items, int capacity)
        : items_(items), capacity_(capacity) {}

    bool loadFromFile(const string& filename)
    {
        ifstream input(filename);
        if (!input)
        {
            cerr << "Error: Could not open file " << filename << endl;
            return false;
        }

        int n;
        input >> n >> capacity_;

        items_.clear();
        items_.resize(n);

        for (int i = 0; i < n; i++)
        {
            int weight, value;
            input >> weight >> value;
            items_[i] = Item(weight, value);
        }

        input.close();
        return true;
    }

    geneticResult solveWithGA(int popSize = 50, int maxGen = 200,
                         double mutRate = 0.05, int tournSize = 3)
    {
        GeneticAlgorithm ga(popSize, maxGen, mutRate, tournSize);
        return ga.solve(items_, capacity_);
    }

    void printSolution(const geneticResult& result) const
    {
        Individual solution(items_.size());
        solution.getChromosome() = result.getBestChromosome();

        int totalWeight = solution.calculateTotalWeight(items_);

        cout << "Best Profit: " << result.getBestValue() << endl;
        cout << "Total Weight: " << totalWeight << endl;
        cout << "Execution Time: " << fixed << setprecision(4) << result.getExecutionTime() << " seconds" << endl;
    }

    const vector<Item>& getItems() const { return items_; }
    int getCapacity() const { return capacity_; }
};

int main()
{
    srand((unsigned)time(nullptr));

    string dataDirectory = "data";
    vector<string> files = getFilesInDirectory(dataDirectory);

    if (files.empty()) {
        cerr << "There's no files in directory: " << dataDirectory << endl;
        return -1;
    }

    vector<TestResult> results;

    for (const string& filepath : files) {
        string filename = filesystem::path(filepath).filename().string();
        cout << "Analise: " << filename << "..." << endl;

        KnapsackProblem problem;

        if (!problem.loadFromFile(filepath)) {
            cerr << "File skip: " << filename << endl;
            continue;
        }

        geneticResult gaResult = problem.solveWithGA(50, 200, 0.05, 3);

        TestResult testResult;
        testResult.filename = filepath;
        testResult.gaTime = gaResult.getExecutionTime();
        testResult.gaResult = gaResult.getBestValue();

        results.push_back(testResult);

        cout << "  time: " << fixed << setprecision(4) << gaResult.getExecutionTime()
             << " sec, result: " << gaResult.getBestValue() << endl;
    }

    printResultsTable(results);

    saveResultsToCSV(results, "genetic_results.csv");

    return 0;
}
