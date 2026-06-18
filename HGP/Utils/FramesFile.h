#pragma once

#include <complex>
#include <string>
#include <vector>

namespace FramesFile {

// Load per-triangle complex frames for FastHGP reduced mesh F_cb.
// Supports:
//   .fframes  — ASCII, one frame per line: "real imag" (or comma-separated)
//   .mat      — MATLAB v5/v7.2 binary, variable name "frames" (complex or Nx2 real)
bool load(const std::string& path, std::vector<std::complex<double>>& frames);

} // namespace FramesFile
