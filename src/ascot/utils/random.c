/**
 * Implements "random.h".
 */

#include "random.h"
#include "consts.h"
#include "defines.h"
#include "random123/philox.h"
#include "random123/u01fixedpt.h"
#include <math.h>
#include <stdlib.h>

/**
 * Seeded key for the random number generator.
 */
static philox4x32_key_t rng_key;

void random_seed(uint64_t seed)
{
    philox4x32_ukey_t user_key = {{(uint32_t)seed, (uint32_t)(seed >> 32)}};
    rng_key = philox4x32keyinit(user_key);
}

void random_uniform(
    real *rnd, size_t nrnd, const size_t *marker_id, size_t *rng_count,
    size_t index)
{
    size_t idx = 0;
    uint64_t marker = (uint64_t)marker_id[index];
    uint64_t count = (uint64_t)rng_count[index];

    while (idx < nrnd)
    {
        philox4x32_ctr_t counter = {
            {(uint32_t)marker, (uint32_t)(marker >> 32), (uint32_t)count,
             (uint32_t)(count >> 32)}};

        philox4x32_ctr_t result = philox4x32(counter, rng_key);

        size_t nloop = nrnd - idx;
        if (nloop > 4)
            nloop = 4;
        for (size_t i = 0; i < nloop; ++i)
        {
            rnd[idx + i] = (real)u01fixedpt_open_closed_32_float(result.v[i]);
        }

        idx += nloop;
        ++count;
    }
    rng_count[index] = (size_t)count;
}

void random_norm(
    real *rnd, size_t nrnd, const size_t *marker_id, size_t *rng_count,
    size_t index)
{
    size_t idx = 0;
    uint64_t marker = (uint64_t)marker_id[index];
    uint64_t count = (uint64_t)rng_count[index];

    while (idx < nrnd)
    {
        philox4x32_ctr_t counter = {
            {(uint32_t)marker, (uint32_t)(marker >> 32), (uint32_t)count,
             (uint32_t)(count >> 32)}};

        philox4x32_ctr_t result = philox4x32(counter, rng_key);
        real u1, u2, radius, angle, rand[4];

        u1 = u01fixedpt_open_closed_32_double(result.v[0]);
        u2 = u01fixedpt_open_closed_32_double(result.v[1]);
        radius = sqrt(-2.0 * log(u1));
        angle = CONST_2PI * u2;
        rand[0] = radius * cos(angle);
        rand[1] = radius * sin(angle);

        u1 = u01fixedpt_open_closed_32_double(result.v[2]);
        u2 = u01fixedpt_open_closed_32_double(result.v[3]);
        radius = sqrt(-2.0 * log(u1));
        angle = CONST_2PI * u2;
        rand[2] = radius * cos(angle);
        rand[3] = radius * sin(angle);

        size_t nloop = nrnd - idx;
        if (nloop > 4)
            nloop = 4;
        for (size_t i = 0; i < nloop; ++i)
        {
            rnd[idx + i] = rand[i];
        }

        idx += nloop;
        ++count;
    }
    rng_count[index] = (size_t)count;
}
