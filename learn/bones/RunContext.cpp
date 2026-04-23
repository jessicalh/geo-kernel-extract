#include "RunContext.h"
#include "TrajectoryResult.h"

namespace nmr {


RunContext::RunContext(RunConfiguration config,
                       std::filesystem::path output_dir)
    : config_(std::move(config)),
      output_dir_(std::move(output_dir)) {}


void RunContext::AttachExtraResult(std::unique_ptr<TrajectoryResult> extra) {
    extras_.push_back(std::move(extra));
}


void RunContext::SetStride(size_t s) {
    // Zero stride is nonsense; clamp to 1 (= every frame).
    stride_ = (s == 0) ? 1 : s;
}


RunOptions RunContext::PerFrameRunOptions() const {
    RunOptions opts = config_.PerFrameRunOptions();
    opts.aimnet2_model = aimnet2_model_;
    return opts;
}


std::vector<std::unique_ptr<TrajectoryResult>>
RunContext::MoveOutExtraResults() {
    std::vector<std::unique_ptr<TrajectoryResult>> out;
    out.swap(extras_);
    return out;
}

}  // namespace nmr
