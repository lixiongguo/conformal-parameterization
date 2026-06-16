#pragma once

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include "stb_image.h"
#include "stb_image_write.h"

#include "gradient_flow.hpp"

namespace aht_ot {

inline ImageMat imread_gray(const std::string& path) {
    int w = 0, h = 0, ch = 0;
    unsigned char* data = stbi_load(path.c_str(), &w, &h, &ch, 1);
    if (!data) throw std::runtime_error("Cannot open image: " + path);
    ImageMat img(h, w);
    for (int r = 0; r < h; ++r) {
        for (int c = 0; c < w; ++c) {
            img(r, c) = static_cast<double>(data[r * w + c]);
        }
    }
    stbi_image_free(data);
    return img;
}

inline void imwrite_gray(const std::string& path, const ImageMat& img) {
    const int h = static_cast<int>(img.rows());
    const int w = static_cast<int>(img.cols());
    std::vector<unsigned char> buf(static_cast<size_t>(h * w));
    for (int r = 0; r < h; ++r) {
        for (int c = 0; c < w; ++c) {
            double v = std::round(img(r, c));
            v = std::max(0.0, std::min(255.0, v));
            buf[static_cast<size_t>(r * w + c)] = static_cast<unsigned char>(v);
        }
    }
    stbi_write_bmp(path.c_str(), w, h, 1, buf.data());
}

inline void imwrite_gray_jpg(const std::string& path, const ImageMat& img,
                             int quality = 90) {
    const int h = static_cast<int>(img.rows());
    const int w = static_cast<int>(img.cols());
    std::vector<unsigned char> buf(static_cast<size_t>(h * w));
    for (int r = 0; r < h; ++r) {
        for (int c = 0; c < w; ++c) {
            double v = std::round(img(r, c));
            v = std::max(0.0, std::min(255.0, v));
            buf[static_cast<size_t>(r * w + c)] = static_cast<unsigned char>(v);
        }
    }
    stbi_write_jpg(path.c_str(), w, h, 1, buf.data(), quality);
}

inline RGBImage imread_rgb(const std::string& path) {
    int w = 0, h = 0, ch = 0;
    unsigned char* data = stbi_load(path.c_str(), &w, &h, &ch, 3);
    if (!data) throw std::runtime_error("Cannot open image: " + path);
    RGBImage rgb;
    for (int k = 0; k < 3; ++k) rgb[k] = ImageMat(h, w);
    for (int r = 0; r < h; ++r) {
        for (int c = 0; c < w; ++c) {
            for (int k = 0; k < 3; ++k) {
                rgb[k](r, c) =
                    static_cast<double>(data[(r * w + c) * 3 + k]);
            }
        }
    }
    stbi_image_free(data);
    return rgb;
}

inline void imwrite_rgb(const std::string& path, const RGBImage& img) {
    const int h = static_cast<int>(img[0].rows());
    const int w = static_cast<int>(img[0].cols());
    std::vector<unsigned char> buf(static_cast<size_t>(h * w * 3));
    for (int r = 0; r < h; ++r) {
        for (int c = 0; c < w; ++c) {
            for (int k = 0; k < 3; ++k) {
                double v = std::round(img[k](r, c));
                v = std::max(0.0, std::min(255.0, v));
                buf[(r * w + c) * 3 + k] = static_cast<unsigned char>(v);
            }
        }
    }
    stbi_write_bmp(path.c_str(), w, h, 3, buf.data());
}

inline ImageMat histeq(const ImageMat& img) {
    const int rows = static_cast<int>(img.rows());
    const int cols = static_cast<int>(img.cols());
    const int N = rows * cols;
    std::vector<int> hist(256, 0);
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            int v = static_cast<int>(std::round(img(r, c)));
            v = std::max(0, std::min(255, v));
            hist[static_cast<size_t>(v)]++;
        }
    }
    std::vector<double> cdf(256, 0.0);
    double acc = 0.0;
    for (int i = 0; i < 256; ++i) {
        acc += static_cast<double>(hist[i]) / static_cast<double>(N);
        cdf[i] = acc;
    }
    ImageMat out(rows, cols);
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            int v = static_cast<int>(std::round(img(r, c)));
            v = std::max(0, std::min(255, v));
            out(r, c) = cdf[static_cast<size_t>(v)] * 255.0;
        }
    }
    return out;
}

inline void equalize_masses(ImageMat& img0, ImageMat& img1,
                            Equalization method) {
    if (method == Equalization::MassRatio) {
        const double sum0 = img0.sum();
        const double sum1 = img1.sum();
        if (sum0 > kEps) img0 *= (sum1 / sum0);
    } else if (method == Equalization::HistEq) {
        img0 = histeq(img0);
        img1 = histeq(img1);
    }
}

/// Upsample density-level map to pixel coordinates (bilinear per block)
inline Tensor2D upsample_mapping(const Tensor2D& u_density,
                                 const ImageMat& image, int block) {
    const int img_h = static_cast<int>(image.rows());
    const int img_w = static_cast<int>(image.cols());
    const int d_h = static_cast<int>(u_density.first.rows());
    const int d_w = static_cast<int>(u_density.first.cols());

    Tensor2D u_img;
    u_img.first = Eigen::MatrixXd::Zero(img_h, img_w);
    u_img.second = Eigen::MatrixXd::Zero(img_h, img_w);

    for (int i = 0; i < d_h; ++i) {
        for (int j = 0; j < d_w; ++j) {
            const int pi = i * block;
            const int pj = j * block;
            if (pi < img_h && pj < img_w) {
                u_img.first(pi, pj) =
                    1.0 + (u_density.first(i, j) - 1.0) * block;
                u_img.second(pi, pj) =
                    1.0 + (u_density.second(i, j) - 1.0) * block;
            }
        }
    }

    int ul_x = 0, ul_y = 0, x_cnt = 1, y_cnt = 1;
    for (int i = 0; i < img_h - block; ++i) {
        y_cnt = 1;
        ul_y = 0;
        for (int j = 0; j < img_w - block; ++j) {
            if (u_img.first(i, j) < kEps) {
                const double rx = static_cast<double>(i - ul_x) / block;
                const double ry = static_cast<double>(j - ul_y) / block;
                const int ulx_nx = std::min(ul_x + block, img_h - 1);
                const int uly_ny = std::min(ul_y + block, img_w - 1);
                double x_loc = (1.0 - rx) * u_img.first(ul_x, ul_y) +
                               rx * u_img.first(ulx_nx, ul_y);
                double y_loc = (1.0 - ry) * u_img.second(ul_x, ul_y) +
                               ry * u_img.second(ul_x, uly_ny);
                x_loc = std::max(1.0, std::min(static_cast<double>(img_h), x_loc));
                y_loc = std::max(1.0, std::min(static_cast<double>(img_w), y_loc));
                u_img.first(i, j) = x_loc;
                u_img.second(i, j) = y_loc;
            }
            ++y_cnt;
            if (y_cnt > block) {
                ul_y += block;
                y_cnt = 1;
            }
        }
        ++x_cnt;
        if (x_cnt > block) {
            ul_x += block;
            x_cnt = 1;
        }
    }

    for (int i = 0; i < img_h; ++i) {
        for (int j = 0; j < img_w; ++j) {
            if (u_img.first(i, j) < kEps) {
                u_img.first(i, j) = static_cast<double>(i + 1);
                u_img.second(i, j) = static_cast<double>(j + 1);
            }
        }
    }
    return u_img;
}

inline ImageMat warp_image(const Tensor2D& u_density, const ImageMat& source,
                           int block) {
    const Tensor2D u_img = upsample_mapping(u_density, source, block);
    const int h = static_cast<int>(source.rows());
    const int w = static_cast<int>(source.cols());
    ImageMat out = ImageMat::Zero(h, w);
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            const int xl = std::max(
                0, std::min(h - 1, static_cast<int>(std::round(u_img.first(i, j))) - 1));
            const int yl = std::max(
                0, std::min(w - 1, static_cast<int>(std::round(u_img.second(i, j))) - 1));
            out(i, j) = source(xl, yl);
        }
    }
    const int crop = block + 4;
    const int nh = std::max(1, h - crop);
    const int nw = std::max(1, w - crop);
    return out.block(0, 0, nh, nw);
}

inline ImageMat warp_blend(const Tensor2D& u_density, const ImageMat& img0,
                           const ImageMat& img1, int block, double t) {
    const Tensor2D u_img = upsample_mapping(u_density, img1, block);
    const int h = static_cast<int>(img1.rows());
    const int w = static_cast<int>(img1.cols());
    ImageMat out = ImageMat::Zero(h, w);
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            const int xl = std::max(
                0, std::min(h - 1, static_cast<int>(std::round(u_img.first(i, j))) - 1));
            const int yl = std::max(
                0, std::min(w - 1, static_cast<int>(std::round(u_img.second(i, j))) - 1));
            out(i, j) = (1.0 - t) * img0(i, j) + t * img1(xl, yl);
        }
    }
    const int crop = block + 4;
    return out.block(0, 0, std::max(1, h - crop), std::max(1, w - crop));
}

inline std::vector<ImageMat> create_morph_sequence(const Tensor2D& optimal_u,
                                                   const ImageMat& img0,
                                                   const ImageMat& img1,
                                                   int block, int num_frames) {
    const int s1 = static_cast<int>(optimal_u.first.rows());
    const int s2 = static_cast<int>(optimal_u.first.cols());
    std::vector<ImageMat> series;
    series.reserve(static_cast<size_t>(num_frames));

    for (int f = num_frames; f >= 1; --f) {
        Tensor2D u_frame;
        u_frame.first = Eigen::MatrixXd(s1, s2);
        u_frame.second = Eigen::MatrixXd(s1, s2);
        for (int i = 0; i < s1; ++i) {
            for (int j = 0; j < s2; ++j) {
                const double frac =
                    static_cast<double>(f - 1) / static_cast<double>(num_frames);
                u_frame.first(i, j) =
                    optimal_u.first(i, j) - frac * (optimal_u.first(i, j) - (i + 1.0));
                u_frame.second(i, j) =
                    optimal_u.second(i, j) - frac * (optimal_u.second(i, j) - (j + 1.0));
            }
        }
        const double t = static_cast<double>(f) / static_cast<double>(num_frames);
        series.push_back(warp_blend(u_frame, img0, img1, block, t));
    }
    return series;
}

/// Alias matching the original MATLAB/C++ port naming.
inline ImageMat transform(const Tensor2D& u, const ImageMat& source, int block) {
    return warp_image(u, source, block);
}

/// Grayscale OT map applied to all RGB channels.
inline std::vector<RGBImage> create_morph_sequence_rgb(const Tensor2D& optimal_u,
                                                     const RGBImage& img0,
                                                     const RGBImage& img1,
                                                     int block, int num_frames) {
    const int s1 = static_cast<int>(optimal_u.first.rows());
    const int s2 = static_cast<int>(optimal_u.first.cols());
    const int h = static_cast<int>(img1[0].rows());
    const int w = static_cast<int>(img1[0].cols());
    const int crop = block + 4;
    const int nh = std::max(1, h - crop);
    const int nw = std::max(1, w - crop);

    std::vector<RGBImage> series;
    series.reserve(static_cast<size_t>(num_frames));

    for (int f = num_frames; f >= 1; --f) {
        Tensor2D u_frame;
        u_frame.first = Eigen::MatrixXd(s1, s2);
        u_frame.second = Eigen::MatrixXd(s1, s2);
        for (int i = 0; i < s1; ++i) {
            for (int j = 0; j < s2; ++j) {
                u_frame.first(i, j) = optimal_u.first(i, j) -
                    (f - 1) * (optimal_u.first(i, j) - (i + 1.0)) / num_frames;
                u_frame.second(i, j) = optimal_u.second(i, j) -
                    (f - 1) * (optimal_u.second(i, j) - (j + 1.0)) / num_frames;
            }
        }
        const Tensor2D u_img = upsample_mapping(u_frame, img1[0], block);
        const double t = static_cast<double>(f) / static_cast<double>(num_frames);

        RGBImage frame;
        for (int k = 0; k < 3; ++k) {
            frame[k] = ImageMat::Zero(h, w);
            for (int i = 0; i < h; ++i) {
                for (int j = 0; j < w; ++j) {
                    const int xl = std::max(
                        0, std::min(h - 1,
                                    static_cast<int>(std::round(u_img.first(i, j))) - 1));
                    const int yl = std::max(
                        0, std::min(w - 1,
                                    static_cast<int>(std::round(u_img.second(i, j))) - 1));
                    frame[k](i, j) =
                        (1.0 - t) * img0[k](i, j) + t * img1[k](xl, yl);
                }
            }
            frame[k] = frame[k].block(0, 0, nh, nw);
        }
        series.push_back(std::move(frame));
    }
    return series;
}

}  // namespace aht_ot
