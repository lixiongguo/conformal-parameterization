#pragma once

#include <stdexcept>
#include <string>

#define FASTHGP_NOT_IMPLEMENTED(msg) \
    throw std::runtime_error(std::string("NotImplemented: ") + (msg))
