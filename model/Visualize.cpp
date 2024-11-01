// Visualize.cpp
#include "Visualize.h"
#include <fstream>

Visualize::Visualize(std::string id) {
    this->id = id;
}

Visualize& Visualize::init() {
    std::cout << id << std::endl;
    return *this;
}

void Visualize::addTrace(std::vector<my_Vector> config) {
    this->traces.push_back(config);
}

void Visualize::extendsTraceByVec(my_Vector vec) {
    if (!traces.empty()) {
        traces.back().push_back(vec);
    }
}

void Visualize::printTraces() {
    for (const auto& trace : traces) {
        for (const auto& vec : trace) {
            vec.print();
        }
        std::cout << "---" << std::endl;
    }
}

void Visualize::saveTraceToFile(const std::string& filename) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Не удалось открыть файл " << filename << " для записи." << std::endl;
        return;
    }

    if (traces.empty()) {
        std::cerr << "Нет данных для сохранения." << std::endl;
        return;
    }

    // Предполагаем, что хотим сохранить первую траекторию
    const auto& trace = traces[0];
    for (const auto& vec : trace) {
        if (vec.vec.size() >= 3) {
            outfile << vec.get(0) << " " << vec.get(1) << " " << vec.get(2) << std::endl;
        }
    }

    outfile.close();
}

const std::vector<std::vector<my_Vector>>& Visualize::getTraces() const {
    return traces;
}
