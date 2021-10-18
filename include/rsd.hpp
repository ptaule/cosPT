/*
   rsd.hpp

   Created by Petter Taule on 18.10.2021
   Copyright (c) 2021 Petter Taule. All rights reserved.
*/

#ifndef RSD_HPP
#define RSD_HPP

namespace rsd {
int integrand(
        __attribute__((unused)) const int *ndim,
        const double xx[],
        __attribute__((unused)) const int *ncomp,
        double ff[],
        void *userdata,
        __attribute__((unused)) const int *nvec,
        const int *core
        );
} /* namespace rsd */

#endif /* !RSD_HPP */
