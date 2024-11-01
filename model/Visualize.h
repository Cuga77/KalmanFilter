// Visualize.h
#pragma once

#include <vector>
#include <string>
#include <iostream>
#include "vector.h"

class Visualize {
    std::vector<std::vector<my_Vector>> traces;
    std::string id;
public:
    Visualize(std::string id);
    Visualize& init();
    void addTrace(std::vector<my_Vector> config);
    void extendsTraceByVec(my_Vector vec);
    void printTraces();
    void saveTraceToFile(const std::string& filename) const;
    const std::vector<std::vector<my_Vector>>& getTraces() const;
};
