#ifndef SASSY_CONFIGURATION_H
#define SASSY_CONFIGURATION_H

#include <functional>

namespace sassy {
// hook for automorphisms found in dejavu
// void sassy_hook(int n, const int *perm, int nsupp, const int *support)
// int n:        domain size of the graph and automorphism group
// int* perm:    permutation in one-line notation, i.e., perm[i] = j if and only if the permutation maps i to j
// int nsupp:    number of vertices moved by the permutation, i.e., support of the permutation and length of int* support array
// int* support: vertices moved by the permutation
// IMPORTANT NOTE: Try to avoid sequential reads of perm, rather use the support array to only access those parts of the
// permutation that are non-trivial.
    //typedef void sassy_hook(int, const int *, int, const int *);
    typedef const std::function<void(int, const int *, int, const int *)> sassy_hook;

    struct configstruct {
        bool CONFIG_PRINT = false;
        bool CONFIG_IR_FULL_INVARIANT = false; // uses a complete invariant and no certification if enabled
        bool CONFIG_IR_IDLE_SKIP = true;  // blueprints
        bool CONFIG_IR_INDIVIDUALIZE_EARLY = false; // experimental feature, based on an idea by Adolfo Piperno
        bool CONFIG_PREP_DEACT_PROBE = false; // preprocessor: no probing
        bool CONFIG_PREP_DEACT_DEG01 = false; // preprocessor: no degree 0,1 processing
        bool CONFIG_PREP_DEACT_DEG2 = false;  // preprocessor: no degree 2   processing
        bool CONFIG_IR_REFINE_EARLYOUT_LATE = false;
        bool CONFIG_TRANSLATE_ONLY = false;
    };
}

#endif //SASSY_CONFIGURATION_H
