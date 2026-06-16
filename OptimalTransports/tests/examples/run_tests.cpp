// Regression tests from original run_tests.m (Haker et al. MATLAB port)
// Usage: run_tests <test_type 1-7> [image_dir]

#include <aht_ot/aht_ot.hpp>

#include <cstdio>
#include <cstdlib>
#include <string>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

static std::string g_image_dir;
static const std::string k_output_dir = "output/";

static void ensure_output_dir() {
#ifdef _WIN32
    _mkdir(k_output_dir.c_str());
#else
    mkdir(k_output_dir.c_str(), 0755);
#endif
}

static std::string img(const std::string& name) { return g_image_dir + name; }
static std::string out(const std::string& name) { return k_output_dir + name; }

static aht_ot::Tensor2D solve_ot(const std::string& p0, const std::string& p1,
                                 int block, aht_ot::PType ptype,
                                 aht_ot::Equalization eq) {
    aht_ot::AHTOptions opt;
    opt.square_edge_len = block;
    opt.p_type = ptype;
    opt.equalization = eq;
    return aht_ot::solve_files(p0, p1, opt).optimal_map;
}

static void test_1() {
    printf("\n===== Test 1: initial mapping u0 (x-axis) =====\n");
    const int block = 12;
    auto image_1 = aht_ot::imread_gray(img("test_u0_ax_2.bmp"));
    auto myu_0 = aht_ot::compute_density_map(
        aht_ot::imread_gray(img("test_u0_ax_1.bmp")), block);
    auto myu_1 = aht_ot::compute_density_map(image_1, block);
    auto u0 = aht_ot::compute_initial_mapping(myu_0, myu_1);
    aht_ot::imwrite_gray(out("test1_deformed.bmp"),
                         aht_ot::transform(u0, image_1, block));
}

static void test_2() {
    printf("\n===== Test 2: initial mapping u0 (y-axis) =====\n");
    const int block = 16;
    auto image_1 = aht_ot::imread_gray(img("test_u0_bxy_1.bmp"));
    auto myu_0 = aht_ot::compute_density_map(
        aht_ot::imread_gray(img("test_u0_bxy_2.bmp")), block);
    auto myu_1 = aht_ot::compute_density_map(image_1, block);
    auto u0 = aht_ot::compute_initial_mapping(myu_0, myu_1);
    aht_ot::imwrite_gray(out("test2_deformed.bmp"),
                         aht_ot::transform(u0, image_1, block));
}

static void test_3() {
    printf("\n===== Test 3: initial mapping u0 (y-axis, histeq) =====\n");
    const int block = 16;
    auto image_0 = aht_ot::histeq(aht_ot::imread_gray(img("test_u0_bxy_2.bmp")));
    auto image_1 = aht_ot::histeq(aht_ot::imread_gray(img("test_u0_bxy_1.bmp")));
    auto myu_0 = aht_ot::compute_density_map(image_0, block);
    auto myu_1 = aht_ot::compute_density_map(image_1, block);
    auto u0 = aht_ot::compute_initial_mapping(myu_0, myu_1);
    image_1 = aht_ot::imread_gray(img("test_u0_bxy_1.bmp"));
    aht_ot::imwrite_gray(out("test3_deformed.bmp"),
                         aht_ot::transform(u0, image_1, block));
}

static void test_4() {
    printf("\n===== Test 4: flame colour deformation series =====\n");
    const int block = 2;
    auto optimal_u = solve_ot(img("flame_2.bmp"), img("flame_1.bmp"), block,
                              aht_ot::PType::WithComparison,
                              aht_ot::Equalization::MassRatio);
    auto series = aht_ot::create_morph_sequence_rgb(
        optimal_u, aht_ot::imread_rgb(img("flame_2.bmp")),
        aht_ot::imread_rgb(img("flame_1.bmp")), block, 20);
    for (int i = 0; i < static_cast<int>(series.size()); ++i) {
        aht_ot::imwrite_rgb(out("test4_frame_" + std::to_string(i + 1) + ".bmp"),
                            series[static_cast<size_t>(i)]);
    }
}

static void test_5() {
    printf("\n===== Test 5: cloud image deformation series =====\n");
    const int block = 4;
    auto optimal_u = solve_ot(img("cloud_2.bmp"), img("cloud_1.bmp"), block,
                              aht_ot::PType::WithComparison,
                              aht_ot::Equalization::MassRatio);
    auto series = aht_ot::create_morph_sequence(
        optimal_u, aht_ot::imread_gray(img("cloud_2.bmp")),
        aht_ot::imread_gray(img("cloud_1.bmp")), block, 6);
    for (int i = 0; i < static_cast<int>(series.size()); ++i) {
        aht_ot::imwrite_gray(out("test5_frame_" + std::to_string(i + 1) + ".bmp"),
                             series[static_cast<size_t>(i)]);
    }
}

static void test_6() {
    printf("\n===== Test 6: water image deformation series =====\n");
    const int block = 2;
    auto optimal_u = solve_ot(img("water_1.jpg"), img("water_2.jpg"), block,
                              aht_ot::PType::WithComparison,
                              aht_ot::Equalization::MassRatio);
    auto series = aht_ot::create_morph_sequence(
        optimal_u, aht_ot::imread_gray(img("water_1.jpg")),
        aht_ot::imread_gray(img("water_2.jpg")), block, 6);
    for (int i = 0; i < static_cast<int>(series.size()); ++i) {
        aht_ot::imwrite_gray_jpg(out("test6_frame_" + std::to_string(i + 1) + ".jpg"),
                                 series[static_cast<size_t>(i)]);
    }
}

static void test_7() {
    printf("\n===== Test 7: cloud deformation (no intensity mixing) =====\n");
    const int block = 4;
    auto optimal_u_0 = solve_ot(img("cloud_1.bmp"), img("cloud_2.bmp"), block,
                                aht_ot::PType::WithComparison,
                                aht_ot::Equalization::MassRatio);
    auto image_1 = aht_ot::imread_gray(img("cloud_2.bmp"));
    aht_ot::imwrite_gray(out("test7_deformed_eq0.bmp"),
                         aht_ot::transform(optimal_u_0, image_1, block));

    auto optimal_u_1 = solve_ot(img("cloud_1.bmp"), img("cloud_2.bmp"), block,
                                aht_ot::PType::WithComparison,
                                aht_ot::Equalization::HistEq);
    image_1 = aht_ot::imread_gray(img("cloud_2.bmp"));
    aht_ot::imwrite_gray(out("test7_deformed_eq1.bmp"),
                         aht_ot::transform(optimal_u_1, image_1, block));
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: run_tests <test_type 1-7> [image_dir]\n");
        fprintf(stderr, "  image_dir defaults to testdata/ (Midas journal BMPs)\n");
        return 1;
    }
    g_image_dir = (argc >= 3) ? argv[2] : "testdata/";
    if (g_image_dir.back() != '/' && g_image_dir.back() != '\\') g_image_dir += '/';

    const int test_type = std::atoi(argv[1]);
    ensure_output_dir();

    try {
        switch (test_type) {
            case 1: test_1(); break;
            case 2: test_2(); break;
            case 3: test_3(); break;
            case 4: test_4(); break;
            case 5: test_5(); break;
            case 6: test_6(); break;
            case 7: test_7(); break;
            default:
                fprintf(stderr, "Invalid test type: %d (must be 1-7)\n", test_type);
                return 1;
        }
    } catch (const std::exception& e) {
        fprintf(stderr, "Error: %s\n", e.what());
        return 1;
    }
    printf("\nDone.\n");
    return 0;
}
