#include "BlessCompare.h"
#include "NpyWriter.h"

#include <gtest/gtest.h>

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <string>
#include <vector>

#include <unistd.h>

namespace fs = std::filesystem;

namespace nmr::test {
namespace {

class BlessCompareTest : public ::testing::Test {
protected:
    void SetUp() override {
        const auto nonce = std::chrono::steady_clock::now()
                               .time_since_epoch()
                               .count();
        dir_ = fs::temp_directory_path()
             / ("nmr_bless_compare_" + std::to_string(nonce));
        fs::create_directories(dir_);
    }

    void TearDown() override {
        // Avoid std::filesystem::remove_all here: libtorch exports a second
        // std::filesystem implementation in this test binary's runtime.
        static constexpr const char* kFiles[] = {
            "run.npy", "blessed.npy",
            "signed_run.npy", "signed_blessed.npy",
            "unsigned_run.npy", "unsigned_blessed.npy",
        };
        for (const char* filename : kFiles) {
            const std::string path = Path(filename).string();
            ::unlink(path.c_str());
        }
        const std::string directory = dir_.string();
        ::rmdir(directory.c_str());
    }

    fs::path Path(const std::string& filename) const {
        return dir_ / filename;
    }

    BlessPolicy NoSparsityCheck() const {
        BlessPolicy policy;
        policy.min_nonzero_fraction = 0.0;
        return policy;
    }

    fs::path dir_;
};

TEST_F(BlessCompareTest, NanAndFiniteAreAlwaysDrift) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const std::vector<double> run{nan, 2.0};
    const std::vector<double> blessed{1.0, 2.0};
    ASSERT_TRUE(NpyWriter::WriteFloat64(Path("run.npy").string(),
                                        run.data(), run.size()));
    ASSERT_TRUE(NpyWriter::WriteFloat64(Path("blessed.npy").string(),
                                        blessed.data(), blessed.size()));

    BlessPolicy policy = NoSparsityCheck();
    policy.atol = std::numeric_limits<double>::infinity();
    policy.rtol = std::numeric_limits<double>::infinity();
    const BlessResult result = CompareNpy(Path("run.npy").string(),
                                          Path("blessed.npy").string(),
                                          policy);

    EXPECT_EQ(result.verdict, BlessVerdict::Drifted);
    EXPECT_NE(result.diagnostic.find("1 non-finite mismatches"),
              std::string::npos);
}

TEST_F(BlessCompareTest, FiniteAndNanAreAlwaysDrift) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const std::vector<double> run{1.0, 2.0};
    const std::vector<double> blessed{nan, 2.0};
    ASSERT_TRUE(NpyWriter::WriteFloat64(Path("run.npy").string(),
                                        run.data(), run.size()));
    ASSERT_TRUE(NpyWriter::WriteFloat64(Path("blessed.npy").string(),
                                        blessed.data(), blessed.size()));

    const BlessResult result = CompareNpy(Path("run.npy").string(),
                                          Path("blessed.npy").string(),
                                          NoSparsityCheck());

    EXPECT_EQ(result.verdict, BlessVerdict::Drifted);
}

TEST_F(BlessCompareTest, NanSentinelsAgreeWithNanSentinels) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const std::vector<double> run{nan, 1.0 + 1.0e-9};
    const std::vector<double> blessed{nan, 1.0};
    ASSERT_TRUE(NpyWriter::WriteFloat64(Path("run.npy").string(),
                                        run.data(), run.size()));
    ASSERT_TRUE(NpyWriter::WriteFloat64(Path("blessed.npy").string(),
                                        blessed.data(), blessed.size()));

    const BlessResult result = CompareNpy(Path("run.npy").string(),
                                          Path("blessed.npy").string(),
                                          NoSparsityCheck());

    EXPECT_EQ(result.verdict, BlessVerdict::WithinTolerance);
}

TEST_F(BlessCompareTest, ReadsSignedAndUnsignedByteDtypes) {
    const std::vector<std::int8_t> signed_run{-2, 0, 3};
    const std::vector<std::int8_t> signed_blessed{-1, 0, 3};
    const std::vector<std::uint8_t> unsigned_run{0, 1, 255};
    const std::vector<std::uint8_t> unsigned_blessed{0, 1, 254};
    ASSERT_TRUE(NpyWriter::WriteInt8(Path("signed_run.npy").string(),
                                     signed_run.data(), signed_run.size()));
    ASSERT_TRUE(NpyWriter::WriteInt8(Path("signed_blessed.npy").string(),
                                     signed_blessed.data(),
                                     signed_blessed.size()));
    ASSERT_TRUE(NpyWriter::WriteUInt8(Path("unsigned_run.npy").string(),
                                      unsigned_run.data(),
                                      unsigned_run.size()));
    ASSERT_TRUE(NpyWriter::WriteUInt8(Path("unsigned_blessed.npy").string(),
                                      unsigned_blessed.data(),
                                      unsigned_blessed.size()));

    BlessPolicy policy = NoSparsityCheck();
    policy.atol = 1.0;
    policy.rtol = 0.0;
    const BlessResult signed_result = CompareNpy(
        Path("signed_run.npy").string(),
        Path("signed_blessed.npy").string(), policy);
    const BlessResult unsigned_result = CompareNpy(
        Path("unsigned_run.npy").string(),
        Path("unsigned_blessed.npy").string(), policy);

    EXPECT_EQ(signed_result.verdict, BlessVerdict::WithinTolerance)
        << signed_result.diagnostic;
    EXPECT_EQ(unsigned_result.verdict, BlessVerdict::WithinTolerance)
        << unsigned_result.diagnostic;
}

}  // namespace
}  // namespace nmr::test
