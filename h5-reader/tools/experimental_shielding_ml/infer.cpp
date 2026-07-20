#include <torch/cuda.h>
#include <torch/script.h>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace fs = std::filesystem;

namespace {

const auto kStartTime = std::chrono::steady_clock::now();

void Log(const std::string& message) {
    const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - kStartTime);
    std::cerr << "[shielding-infer +" << elapsed.count() << "ms] "
              << message << std::endl;
}

std::map<std::string, std::int64_t> LoadShapes(const fs::path& inputDir) {
    std::ifstream input(inputDir / "shapes.txt");
    if (!input)
        throw std::runtime_error("could not open shapes.txt");

    std::map<std::string, std::int64_t> shapes;
    std::string key;
    std::int64_t value = 0;
    while (input >> key >> value) {
        if (!shapes.emplace(key, value).second)
            throw std::runtime_error("duplicate shape key: " + key);
    }
    for (const std::string& required :
         {"N", "E", "C1", "C2", "C0", "A", "label_count", "radial_dim"}) {
        if (shapes.find(required) == shapes.end())
            throw std::runtime_error("missing shape key: " + required);
    }
    if (shapes.at("N") <= 0 || shapes.at("E") < 0 || shapes.at("C1") <= 0
        || shapes.at("C2") <= 0 || shapes.at("C0") <= 0
        || shapes.at("A") <= 0 || shapes.at("label_count") <= 0
        || shapes.at("radial_dim") <= 0) {
        throw std::runtime_error("invalid non-positive input shape");
    }
    return shapes;
}

std::int64_t CheckedProduct(std::initializer_list<std::int64_t> dimensions) {
    std::int64_t product = 1;
    for (const std::int64_t dimension : dimensions) {
        if (dimension < 0
            || (dimension != 0
                && product > std::numeric_limits<std::int64_t>::max() / dimension)) {
            throw std::runtime_error("input shape overflows");
        }
        product *= dimension;
    }
    return product;
}

template <typename T>
std::vector<T> ReadBinaryVector(const fs::path& path, std::int64_t count) {
    const std::int64_t bytes =
        CheckedProduct({count, static_cast<std::int64_t>(sizeof(T))});
    if (!fs::is_regular_file(path))
        throw std::runtime_error("input is absent: " + path.string());
    if (fs::file_size(path) != static_cast<std::uintmax_t>(bytes)) {
        throw std::runtime_error(
            "input byte count differs for " + path.filename().string());
    }

    std::vector<T> data(static_cast<std::size_t>(count));
    if (bytes == 0)
        return data;
    std::ifstream input(path, std::ios::binary);
    if (!input)
        throw std::runtime_error("could not open " + path.string());
    input.read(
        reinterpret_cast<char*>(data.data()),
        static_cast<std::streamsize>(bytes));
    if (!input || input.gcount() != static_cast<std::streamsize>(bytes))
        throw std::runtime_error("short read from " + path.string());
    return data;
}

torch::Tensor LoadFloatTensor(
    const fs::path& inputDir,
    const std::string& name,
    std::initializer_list<std::int64_t> dimensions,
    const torch::Device& device) {
    const std::int64_t count = CheckedProduct(dimensions);
    Log("read " + name);
    std::vector<float> data =
        ReadBinaryVector<float>(inputDir / (name + ".bin"), count);
    return torch::from_blob(
               data.data(),
               std::vector<std::int64_t>(dimensions),
               torch::kFloat32)
        .clone()
        .to(device);
}

torch::Tensor LoadLongTensor(
    const fs::path& inputDir,
    const std::string& name,
    std::initializer_list<std::int64_t> dimensions,
    const torch::Device& device) {
    const std::int64_t count = CheckedProduct(dimensions);
    Log("read " + name);
    std::vector<std::int64_t> data =
        ReadBinaryVector<std::int64_t>(inputDir / (name + ".bin"), count);
    return torch::from_blob(
               data.data(),
               std::vector<std::int64_t>(dimensions),
               torch::kInt64)
        .clone()
        .to(device);
}

void WriteOutput(const torch::Tensor& tensor, const fs::path& path) {
    const torch::Tensor output =
        tensor.to(torch::kCPU, torch::kFloat32).contiguous();
    std::ofstream file(path, std::ios::binary);
    if (!file)
        throw std::runtime_error("could not create " + path.string());
    file.write(
        reinterpret_cast<const char*>(output.data_ptr<float>()),
        static_cast<std::streamsize>(output.numel() * sizeof(float)));
    if (!file)
        throw std::runtime_error("short write to " + path.string());
}

std::string LowerAscii(std::string text) {
    std::transform(
        text.begin(),
        text.end(),
        text.begin(),
        [](unsigned char value) {
            return static_cast<char>(std::tolower(value));
        });
    return text;
}

void LoadAcceleratorBackend(const std::string& requested) {
    const std::string device = LowerAscii(requested);
    if (device != "rocm" && device != "cuda")
        return;
#ifdef _WIN32
    if (!LoadLibraryW(L"torch_hip.dll")) {
        throw std::runtime_error(
            "could not load torch_hip.dll (Windows error "
            + std::to_string(GetLastError()) + ")");
    }
#else
    (void)requested;
#endif
}

torch::Device ResolveDevice(const std::string& requested) {
    const std::string device =
        LowerAscii(requested.empty() ? "cpu" : requested);
    if (device == "cpu")
        return torch::Device(torch::kCPU);
    if (device == "rocm" || device == "cuda") {
        if (!torch::cuda::is_available()) {
            throw std::runtime_error(
                device + " requested, but the accelerator is unavailable");
        }
        return torch::Device(torch::kCUDA);
    }
    throw std::runtime_error("device must be cpu, rocm, or cuda");
}

int Run(int argc, char** argv) {
    try {
        if (argc < 4 || argc > 5) {
            std::cerr
                << "usage: infer <model.ts> <input_dir> <output.bin> "
                   "[cpu|rocm]\n";
            return 64;
        }

        const fs::path modelPath(argv[1]);
        const fs::path inputDir(argv[2]);
        const fs::path outputPath(argv[3]);
        const std::string requestedDevice = argc == 5 ? argv[4] : "cpu";
        LoadAcceleratorBackend(requestedDevice);
        const torch::Device device = ResolveDevice(requestedDevice);
        Log("device " + device.str());

        Log("load model");
        torch::jit::Module module =
            torch::jit::load(modelPath.string(), device);
        module.eval();

        const auto shapes = LoadShapes(inputDir);
        const std::int64_t n = shapes.at("N");
        const std::int64_t e = shapes.at("E");
        const std::int64_t c1 = shapes.at("C1");
        const std::int64_t c2 = shapes.at("C2");
        const std::int64_t c0 = shapes.at("C0");
        const std::int64_t applicability = shapes.at("A");
        const std::int64_t labelCount = shapes.at("label_count");
        const std::int64_t radialDim = shapes.at("radial_dim");

        std::vector<torch::jit::IValue> inputs = {
            LoadFloatTensor(inputDir, "pos", {n, 3}, device),
            LoadFloatTensor(inputDir, "l1", {n, c1, 3}, device),
            LoadFloatTensor(inputDir, "l1_valid", {n, c1}, device),
            LoadFloatTensor(inputDir, "l2", {n, c2, 5}, device),
            LoadFloatTensor(inputDir, "l2_valid", {n, c2}, device),
            LoadFloatTensor(inputDir, "scalars", {n, c0}, device),
            LoadFloatTensor(inputDir, "scalar_valid", {n, c0}, device),
            LoadFloatTensor(
                inputDir,
                "applicability",
                {n, applicability},
                device),
            LoadLongTensor(
                inputDir,
                "label_ids",
                {n, labelCount},
                device),
            LoadLongTensor(inputDir, "edge_src", {e}, device),
            LoadLongTensor(inputDir, "edge_dst", {e}, device),
            LoadFloatTensor(
                inputDir,
                "radial",
                {e, radialDim},
                device),
        };

        Log("forward");
        torch::InferenceMode guard;
        const torch::Tensor output = module.forward(inputs).toTensor();
        if (output.dim() != 2 || output.size(0) != n || output.size(1) != 6)
            throw std::runtime_error("model output is not N x 6");
        if (!torch::isfinite(output).all().item<bool>())
            throw std::runtime_error("model output contains non-finite values");

        WriteOutput(output, outputPath);
        Log("complete");
        std::cout << "wrote " << outputPath.string() << " shape " << n
                  << " 6 device=" << device.str() << '\n';
        return 0;
    } catch (const c10::Error& error) {
        std::cerr << "LibTorch error: " << error.what() << '\n';
        return 2;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}

}  // namespace

int main(int argc, char** argv) {
    const int result = Run(argc, argv);
    std::cout.flush();
    std::cerr.flush();
    std::_Exit(result);
}
