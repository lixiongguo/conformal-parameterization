/// Morphing sequence via optimal mass transport (sec 7.4)

#include <aht_ot/aht_ot.hpp>

#include <cstdio>
#include <cstdlib>
#include <string>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

static void ensure_dir(const char* path) {
#ifdef _WIN32
    _mkdir(path);
#else
    mkdir(path, 0755);
#endif
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr,
                     "Usage: %s <image0> <image1> [block=8] [frames=8]\n",
                     argv[0]);
        return 1;
    }

    const std::string path0 = argv[1];
    const std::string path1 = argv[2];
    aht_ot::AHTOptions opt;
    opt.p_type = aht_ot::PType::WithComparison;
    if (argc >= 4) opt.square_edge_len = std::atoi(argv[3]);
    const int num_frames = (argc >= 5) ? std::atoi(argv[4]) : 8;

    ensure_dir("output");

    std::printf("=== AHT Morphing Demo ===\n");

    try {
        auto img0 = aht_ot::imread_gray(path0);
        auto img1 = aht_ot::imread_gray(path1);
        const aht_ot::AHTResult result = aht_ot::solve(img0, img1, opt);

        auto series = aht_ot::create_morph_sequence(
            result.optimal_map, img0, img1, opt.square_edge_len, num_frames);

        for (int i = 0; i < static_cast<int>(series.size()); ++i) {
            char fname[128];
            std::snprintf(fname, sizeof(fname), "output/frame_%02d.bmp", i + 1);
            aht_ot::imwrite_gray(fname, series[static_cast<size_t>(i)]);
            std::printf("  wrote %s\n", fname);
        }
    } catch (const std::exception& e) {
        std::fprintf(stderr, "Error: %s\n", e.what());
        return 2;
    }
    return 0;
}
