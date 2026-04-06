#pragma once
#include "ComputeWorker.h"
#include <string>

// ML predictor stub — will be wired to the trained model later.
// Currently returns empty results.

class MLPredictor {
public:
    MLPredictor() = default;
    ~MLPredictor() = default;

    bool loadModels(const std::string& /*modelDir*/) { return false; }
    bool isLoaded() const { return false; }
};
