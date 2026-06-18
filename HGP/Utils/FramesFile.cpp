#include "FramesFile.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

namespace {

bool endsWith(const std::string& s, const std::string& suffix)
{
    return s.size() >= suffix.size()
           && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool loadTextFframes(const std::string& path, std::vector<std::complex<double>>& frames)
{
    std::ifstream in(path);
    if (!in) {
        return false;
    }

    frames.clear();
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        size_t start = 0;
        while (start < line.size() && std::isspace(static_cast<unsigned char>(line[start]))) {
            ++start;
        }
        if (start >= line.size() || line[start] == '#') {
            continue;
        }

        for (char& c : line) {
            if (c == ',') {
                c = ' ';
            }
        }

        std::istringstream iss(line.substr(start));
        double re = 0.0, im = 0.0;
        if (!(iss >> re)) {
            continue;
        }
        if (!(iss >> im)) {
            im = 0.0;
        }
        frames.emplace_back(re, im);
    }

    return !frames.empty();
}

uint32_t readU32(const std::vector<uint8_t>& buf, size_t& pos)
{
    if (pos + 4 > buf.size()) {
        return 0;
    }
    uint32_t v = static_cast<uint32_t>(buf[pos])
                 | (static_cast<uint32_t>(buf[pos + 1]) << 8)
                 | (static_cast<uint32_t>(buf[pos + 2]) << 16)
                 | (static_cast<uint32_t>(buf[pos + 3]) << 24);
    pos += 4;
    return v;
}

void align8(size_t& pos)
{
    pos = (pos + 7) & ~size_t(7);
}

bool readDataElement(const std::vector<uint8_t>& buf, size_t& pos, uint32_t& dataType, std::vector<uint8_t>& payload)
{
    if (pos + 8 > buf.size()) {
        return false;
    }

    uint32_t tag = readU32(buf, pos);
    uint32_t nbytes = readU32(buf, pos);

    if (tag == 0) {
        dataType = nbytes;
        nbytes = readU32(buf, pos);
        readU32(buf, pos);
        align8(pos);
    } else {
        dataType = tag & 0xFFFF;
        if (tag & 0xFFFF0000u) {
            nbytes = tag >> 16;
        }
    }

    if (pos + nbytes > buf.size()) {
        return false;
    }

    payload.assign(buf.begin() + static_cast<std::ptrdiff_t>(pos),
                   buf.begin() + static_cast<std::ptrdiff_t>(pos + nbytes));
    pos += nbytes;
    align8(pos);
    return true;
}

bool parseMiMatrix(const std::vector<uint8_t>& payload, std::vector<std::complex<double>>& frames)
{
    if (payload.size() < 20) {
        return false;
    }

    size_t p = 0;
    uint32_t flags = readU32(payload, p);
    readU32(payload, p);
    const uint32_t classId = flags & 0xFF;
    const bool isComplex = (flags & 0x0800u) != 0;

    if (p + 8 > payload.size()) {
        return false;
    }
    uint32_t ndim = readU32(payload, p);
    readU32(payload, p);
    if (ndim < 1 || ndim > 2) {
        return false;
    }

    std::vector<uint32_t> dims(ndim);
    for (uint32_t d = 0; d < ndim; ++d) {
        dims[d] = readU32(payload, p);
    }

    if (p + 8 > payload.size()) {
        return false;
    }
    uint32_t nameLen = readU32(payload, p);
    readU32(payload, p);
    p += nameLen;
    align8(p);

    size_t count = 1;
    for (uint32_t d : dims) {
        count *= d;
    }
    if (count == 0) {
        return false;
    }

    if (classId != 6) {
        return false;
    }

    const size_t bytesNeeded = count * sizeof(double);
    if (p + bytesNeeded > payload.size()) {
        return false;
    }

    std::vector<double> real(count);
    std::memcpy(real.data(), payload.data() + p, bytesNeeded);
    p += bytesNeeded;
    align8(p);

    std::vector<double> imag;
    if (isComplex) {
        if (p + bytesNeeded > payload.size()) {
            return false;
        }
        imag.resize(count);
        std::memcpy(imag.data(), payload.data() + p, bytesNeeded);
    }

    const uint32_t rows = dims[0];
    const uint32_t cols = ndim > 1 ? dims[1] : 1u;

    frames.clear();

    if (isComplex && cols == 1) {
        frames.resize(rows);
        for (uint32_t i = 0; i < rows; ++i) {
            frames[i] = std::complex<double>(real[i], imag[i]);
        }
        return true;
    }
    if (isComplex && rows == 1) {
        frames.resize(cols);
        for (uint32_t j = 0; j < cols; ++j) {
            frames[j] = std::complex<double>(real[j], imag[j]);
        }
        return true;
    }
    if (!isComplex && cols == 2) {
        frames.resize(rows);
        for (uint32_t i = 0; i < rows; ++i) {
            frames[i] = std::complex<double>(real[i], real[rows + i]);
        }
        return true;
    }
    if (!isComplex && rows == 2) {
        frames.resize(cols);
        for (uint32_t j = 0; j < cols; ++j) {
            frames[j] = std::complex<double>(real[2 * j], real[2 * j + 1]);
        }
        return true;
    }

    frames.resize(count);
    for (size_t i = 0; i < count; ++i) {
        frames[i] = std::complex<double>(real[i], isComplex ? imag[i] : 0.0);
    }
    return true;
}

bool loadMatV5Frames(const std::string& path, std::vector<std::complex<double>>& frames)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        return false;
    }

    std::vector<uint8_t> buf((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    if (buf.size() < 128) {
        return false;
    }
    if (buf[124] != 'I' || buf[125] != 'M') {
        std::cout << "FramesFile: only MATLAB v5/v7.2 .mat is supported (not v7.3 HDF5). "
                     "Export with save('-v7', ...) or use .fframes.\n";
        return false;
    }

    size_t pos = 128;
    while (pos + 8 <= buf.size()) {
        uint32_t dataType = 0;
        std::vector<uint8_t> payload;
        if (!readDataElement(buf, pos, dataType, payload)) {
            break;
        }
        if (dataType == 14) {
            std::vector<std::complex<double>> candidate;
            if (parseMiMatrix(payload, candidate)) {
                frames = std::move(candidate);
                return true;
            }
        }
    }

    return false;
}

} // namespace

namespace FramesFile {

bool load(const std::string& path, std::vector<std::complex<double>>& frames)
{
    frames.clear();

    if (endsWith(path, ".fframes")) {
        if (!loadTextFframes(path, frames)) {
            std::cout << "FramesFile: failed to read .fframes: " << path << "\n";
            return false;
        }
        return true;
    }

    if (endsWith(path, ".mat")) {
        if (loadMatV5Frames(path, frames)) {
            return true;
        }
        if (loadTextFframes(path, frames)) {
            return true;
        }
        std::cout << "FramesFile: failed to read .mat frames from " << path << "\n"
                  << "  Tip: in MATLAB: writematrix([real(frames), imag(frames)], 'frames.fframes');\n"
                  << "  Or: save('-v7', 'frames.mat', 'frames');\n";
        return false;
    }

    std::cout << "FramesFile: unsupported extension (use .fframes or .mat)\n";
    return false;
}

} // namespace FramesFile
