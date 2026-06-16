/// AHT optimal mass transport — command-line demo
///
/// Usage:
///   run_aht_ot <image0> <image1> [block_size] [output.bmp]
///
/// Sec 7.2 initial mapping + sec 7.4 gradient descent (Haker et al. 2004)

#include <aht_ot/aht_ot.hpp>

#include <cstdio>
#include <cstdlib>
#include <string>

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr,
                     "Usage: %s <image0> <image1> [block_size=8] [output.bmp]\n",
                     argv[0]);
        return 1;
    }

    const std::string path0 = argv[1];
    const std::string path1 = argv[2];
    aht_ot::AHTOptions opt;
    if (argc >= 4) opt.square_edge_len = std::atoi(argv[3]);

    std::printf("=== AHT Optimal Mass Transport ===\n");
    std::printf("  image0: %s\n", path0.c_str());
    std::printf("  image1: %s\n", path1.c_str());
    std::printf("  block:  %d\n\n", opt.square_edge_len);

    try {
        const aht_ot::AHTResult result = aht_ot::solve_files(path0, path1, opt);

        auto img0 = aht_ot::imread_gray(path0);
        auto img1 = aht_ot::imread_gray(path1);

        const aht_ot::ImageMat warped_init =
            aht_ot::warp_image(result.initial_map, img1, opt.square_edge_len);
        const aht_ot::ImageMat warped_opt =
            aht_ot::warp_image(result.optimal_map, img1, opt.square_edge_len);

        const std::string out_init = "output_aht_init.bmp";
        const std::string out_opt =
            (argc >= 5) ? argv[4] : std::string("output_aht_opt.bmp");

        aht_ot::imwrite_gray(out_init, warped_init);
        aht_ot::imwrite_gray(out_opt, warped_opt);

        std::printf("\nInitial map warp -> %s\n", out_init.c_str());
        std::printf("Optimal map warp -> %s\n", out_opt.c_str());
        std::printf("Final mean |curl| = %.6f\n", result.final_mean_curl);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "Error: %s\n", e.what());
        return 2;
    }
    return 0;
}
