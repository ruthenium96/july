#include "Spectrum.h"

Spectrum::Spectrum(std::vector<Subspectrum>&& m) : blocks(std::move(m)) {}
