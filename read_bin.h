#pragma once

#include <string_view>

/**
    Bin::load_chunk( bin_name, data );
 */

struct ParticleDataContainter
{

};

namespace Bin
{
void to_txt(std::string name);

ParticleDataContainter load_chunk(std::string_view name);

};
