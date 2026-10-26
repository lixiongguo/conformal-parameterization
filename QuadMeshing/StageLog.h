#ifndef STAGE_LOG_H
#define STAGE_LOG_H

#include <chrono>
#include <string>
#include <utility>
#include <vector>

struct StageTime {
    std::string name;
    double ms = 0;
};

/** Records named stages. run() keeps a false return from the stage body. */
class StageClock {
public:
    explicit StageClock(std::vector<StageTime>& out) : out_(out) {}

    template <typename Fn>
    bool run(const char* name, Fn&& fn) {
        const auto t0 = std::chrono::steady_clock::now();
        const bool ok = static_cast<bool>(fn());
        const double ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - t0).count();
        out_.push_back(StageTime{name, ms});
        return ok;
    }

    template <typename Fn>
    void lap(const char* name, Fn&& fn) {
        run(name, [&] { fn(); return true; });
    }

private:
    std::vector<StageTime>& out_;
};

#endif
