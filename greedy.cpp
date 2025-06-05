#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <filesystem>

struct Item {
    int value;
    int weight;
    double value_per_weight;

    Item(int v, int w) : value(v), weight(w) {
        value_per_weight = static_cast<double>(v) / w;
    }
};

// Функция для чтения данных из файла
bool readDataFromFile(const std::string& filename, int& N, int& W, std::vector<Item>& items) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Не удалось открыть файл: " << filename << std::endl;
        return false;
    }

    file >> N >> W;
    for (int i = 0; i < N; ++i) {
        int v, w;
        file >> v >> w;
        items.emplace_back(v, w);
    }

    file.close();
    return true;
}

// Жадный алгоритм для задачи о рюкзаке
int greedyKnapsack(int W, std::vector<Item>& items) {
    // Сортируем предметы по убыванию ценности на единицу веса
    std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
        return a.value_per_weight > b.value_per_weight;
    });

    int totalValue = 0;
    int currentWeight = 0;

    for (const auto& item : items) {
        if (currentWeight + item.weight <= W) {
            currentWeight += item.weight;
            totalValue += item.value;
        } else {
            int remain = W - currentWeight;
            totalValue += item.value_per_weight * remain;
            break;
        }
    }

    return totalValue;
}

int main() {
    std::string dataPath = "data";
    for (const auto& entry : std::filesystem::directory_iterator(dataPath)) {
        if (entry.is_regular_file()) {
            std::string filename = entry.path().string();
            int N, W;
            std::vector<Item> items;

            if (readDataFromFile(filename, N, W, items)) {
                int maxValue = greedyKnapsack(W, items);
                std::cout << "File: " << filename << " Maximum cost: " << maxValue << std::endl;
            }
        }
    }

    return 0;
}
