#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

extern "C" {
    #include <gsl/gsl_sf.h>
}

#include "../include/utilities.hpp"
#include "../include/combinatorics.hpp"
#include "../include/parameters.hpp"
#include "../include/tables.hpp"
#include "../include/diagrams.hpp"

using std::size_t;

/* Helper function converting from unsigned int */
double gsl_sf_fact(int n) {
    return gsl_sf_fact(static_cast<unsigned int>(n));
}



void PowerSpectrumDiagram::kernel_arguments(size_t rearr_idx, size_t sign_idx)
{
    size_t n_coeffs = loop_params.n_coeffs();

    Vec1D<int>& rearrangement = rearrangements.at(rearr_idx);
    Vec1D<bool>& signs = sign_configs.at(sign_idx);

    Vec1D<int>& arguments_l = arg_configs_l.at(rearr_idx).at(sign_idx).args;
    Vec1D<int>& arguments_r = arg_configs_r.at(rearr_idx).at(sign_idx).args;

    /* First argument : q_m1 = q_m - (±q_m2) - (±q_m3) - ... */
    Vec1D<int> config(n_coeffs, 0);
    config.at(n_coeffs - 1) = 1;

    for (size_t i = 2; i <= static_cast<size_t>(m); ++i) {
        /* Note extra minus sign from kernel expression */
        config[static_cast<size_t>(rearrangement.at(i - 2))] =
            (signs.at(i - 2) ? -1 : 1);
    }

    int label = config2label(config);
    arguments_l.at(0) = label;
    arguments_r.at(0) = label;
    q_m1_labels.at(rearr_idx).at(sign_idx) = label;

    /* Reset config */
    std::fill(config.begin(), config.end(), 0);

    /* Argument indices */
    size_t index_l = 1;
    size_t index_r = 1;

    /* q_m2,q_m3,... arguments, the ordering of which is stored by the first
     * (m-1) entries of rearrangement. Note that the signs are opposite of
     * corresponding terms in the first argument */
    for (size_t i = 2; i <= static_cast<size_t>(m); ++i) {
        config.at(static_cast<size_t>(rearrangement.at(i - 2))) =
            (signs.at(i - 2) ? 1 : -1);
        arguments_l.at(index_l++) = config2label(config);
        arguments_r.at(index_r++) = config2label(config);
        std::fill(config.begin(), config.end(), 0);
    }

    /* l-loop arguments */
    for (size_t i = 0; i < static_cast<size_t>(l); ++i) {
        size_t loop_momentum_index =
            static_cast<size_t>(rearrangement.at(i + static_cast<size_t>(m - 1)));
        config.at(loop_momentum_index) = 1;
        arguments_l.at(index_l++) = config2label(config);
        config.at(loop_momentum_index) = -1;
        arguments_l.at(index_l++) = config2label(config);

        std::fill(config.begin(), config.end(), 0);
    }
    /* r-loop arguments */
    for (size_t i = 0; i < static_cast<size_t>(r); ++i) {
        size_t loop_momentum_index = static_cast<size_t>(
            rearrangement.at(i + static_cast<size_t>(m + l - 1)));
        config.at(loop_momentum_index) = 1;
        arguments_r.at(index_r++) = config2label(config);
        config.at(loop_momentum_index) = -1;
        arguments_r.at(index_r++) = config2label(config);

        std::fill(config.begin(), config.end(), 0);
    }

    /* Fill remaining spots with zero-label */
    int zero_label = loop_params.zero_label();
    size_t n_kernel_args = loop_params.n_kernel_args();
    while (index_l < n_kernel_args) {
        arguments_l.at(index_l++) = zero_label;
    }
    while (index_r < n_kernel_args) {
        arguments_r.at(index_r++) = zero_label;
    }

#if DEBUG >= 1
    /* arguments_2_kernel_index() assumes that the length of arguments[] is
     * n_kernel_args. Checking this explicitly: */
    if (arguments_l.size() != n_kernel_args ||
        arguments_r.size() != n_kernel_args
       ) {
        throw(std::logic_error(
            "PowerSpectrumDiagram::kernel_arguments(): Size of left argument "
            "vector does not equal n_kernel_args."));
    }
#endif

    arg_configs_l.at(rearr_idx).at(sign_idx).kernel_index =
        loop_params.args_2_kernel_index(arguments_l.data());
    arg_configs_r.at(rearr_idx).at(sign_idx).kernel_index =
        loop_params.args_2_kernel_index(arguments_r.data());
}



void PowerSpectrumDiagram::compute_rearrangements(int n_loops)
{
    // Allocate memory for rearrangements
    for (auto& rearr : rearrangements) {
        rearr.resize(static_cast<size_t>(n_loops));
    }

    Vec1D<int> group_sizes;

    /* Connecting loops? */
    if (m > 1) {
        group_sizes.push_back(m-1);
    }
    if (l > 0) {
        group_sizes.push_back(l);
    }
    if (r > 0) {
        group_sizes.push_back(r);
    }

    Orderings orderings(n_loops, group_sizes);

    size_t counter = 0;
    do {
        rearrangements.at(counter) = orderings.get_current();
        ++counter;
    } while (orderings.next());
}



void PowerSpectrumDiagram::compute_sign_flips() {
    // Allocate memory for sign configurations
    for (auto& sign_config : sign_configs) {
        // Signs of k2,...,km can be flipped, hence allocate (m-1)-length array
        sign_config.resize(static_cast<size_t>(m-1));
    }

    // Loop over possible sign configurations and store in sign table
    for (size_t i = 0; i < sign_configs.size(); ++i) {
        for (size_t j = 0; j < static_cast<size_t>(m-1); ++j) {
            // Add 1 so that the first sign configuration is +1,+1,...,+1
            sign_configs.at(i).at(j) =
                (static_cast<int>(i) /
                     static_cast<int>(pow(2, static_cast<int>(j))) + 1)
                % 2;
        }
    }
}



PowerSpectrumDiagram::PowerSpectrumDiagram(
        const LoopParameters& loop_params,
        int m, int l, int r
        ) : loop_params(loop_params), m(m), l(l), r(r)
{
    int n_loops = loop_params.n_loops();
    size_t n_kernel_args = loop_params.n_kernel_args();

    if (m + l + r != n_loops + 1) {
        throw(std::invalid_argument("m + l + r != n_loops + 1"));
    }

    diagram_factor_ = static_cast<int>(
        (gsl_sf_fact(2 * l + m) * gsl_sf_fact(2 * r + m)) /
        (pow(2, l + r) * gsl_sf_fact(l) * gsl_sf_fact(r) * gsl_sf_fact(m)));

    rearrangements.resize(static_cast<size_t>(
        gsl_sf_fact(n_loops) /
        (gsl_sf_fact(m - 1) * gsl_sf_fact(l) * gsl_sf_fact(r))
        ));
    sign_configs.resize(static_cast<size_t>(pow(2,m-1)));

    compute_rearrangements(n_loops);
    compute_sign_flips();

    // Allocate memory
    arg_configs_l.resize(rearrangements.size());
    arg_configs_r.resize(rearrangements.size());

    /* q_m1_label values are set in kernel_arguments function */
    q_m1_labels.resize(rearrangements.size());

    for (size_t a = 0; a < rearrangements.size(); ++a) {
        arg_configs_l.at(a).resize(sign_configs.size());
        arg_configs_r.at(a).resize(sign_configs.size());

        q_m1_labels.at(a).resize(sign_configs.size());

        for (size_t b = 0; b < sign_configs.size(); ++b) {
            arg_configs_l.at(a).at(b).args.resize(n_kernel_args);
            arg_configs_r.at(a).at(b).args.resize(n_kernel_args);

            // Initialize arguments and kernel indices for this configuration
            kernel_arguments(a, b);
        }
    }
}



std::string PowerSpectrumDiagram::argument_configuration(size_t a, size_t b) const {
    return (
        labels2string(arg_configs_l.at(a).at(b).args,
                      loop_params.n_coeffs(), loop_params.spectrum())
        + " "
        + labels2string(arg_configs_r.at(a).at(b).args,
                        loop_params.n_coeffs(), loop_params.spectrum())
    );
}



void BiSpectrumDiagram::compute_rearrangements(int n_loops) {
    // Allocate memory for rearrangements
    for (auto& rearr : rearrangements) {
        rearr.resize(static_cast<size_t>(n_loops));
    }

    Vec1D<int> group_sizes;

    /* Overall loops? */
    if (overall_loop_) {
        group_sizes.push_back(1);
    }
    /* Connecting loops? */
    if (n_ab > 1) {
        group_sizes.push_back(n_ab - 1);
    }
    if (n_bc > 1) {
        group_sizes.push_back(n_bc - 1);
    }
    if (n_ca > 1) {
        group_sizes.push_back(n_ca - 1);
    }
    /* Self loops? */
    if (n_a > 0) {
        group_sizes.push_back(n_a);
    }
    if (n_b > 0) {
        group_sizes.push_back(n_b);
    }
    if (n_c > 0) {
        group_sizes.push_back(n_c);
    }

    Orderings orderings(n_loops, group_sizes);

    size_t counter = 0;
    do {
        rearrangements.at(counter) = orderings.get_current();
        ++counter;
    } while (orderings.next());
}



Vec2D<bool> BiSpectrumDiagram::connecting_line_sign_flips(int connecting_loops)
    const
{
    if (connecting_loops == 0) {
        return Vec2D<bool>();
    }
    size_t n_sign_configs_xy = static_cast<size_t>(pow(2, connecting_loops));

    Vec2D<bool> sign_configs_xy;

    /* Allocate memory for sign configurations */
    sign_configs_xy.resize(n_sign_configs_xy);
    for (size_t i = 0; i < n_sign_configs_xy; ++i) {
        sign_configs_xy.at(i).resize(static_cast<size_t>(connecting_loops));
    }

    /* Loop over possible sign configurations and store */
    for (size_t i = 0; i < n_sign_configs_xy; ++i) {
        for (size_t j = 0; j < static_cast<size_t>(connecting_loops); ++j) {
            // Add 1 so that the first sign configuration is +1,+1,...,+1
            sign_configs_xy.at(i).at(j) =
                (static_cast<int>(i) /
                 static_cast<int>(pow(2, static_cast<int>(j))) + 1)
                % 2;
        }
    }
    return sign_configs_xy;
}



/* Alias functions for Triple<> */
#define ab() a()
#define bc() b()
#define ca() c()



void BiSpectrumDiagram::compute_sign_flips(
        const Triple<Vec2D<bool>>& sign_configs_xy
        )
{
    int tot_connecting_loops = n_connecting_loops.ab() +
        n_connecting_loops.bc() + n_connecting_loops.ca();

    size_t n_sign_configs = static_cast<size_t>(pow(2, tot_connecting_loops));

    /* Allocate memory for sign configurations */
    sign_configs.resize(n_sign_configs);
    for (size_t i = 0; i < n_sign_configs; ++i) {
        sign_configs.at(i).resize(static_cast<size_t>(tot_connecting_loops));
    }

    /* Combine sign_configs_ab,... into sign_configs */
    size_t idx = 0;
    size_t i = 0;
    do {
        size_t j = 0;
        do {
            size_t k = 0;
            do {
                if (!sign_configs_xy.ab().empty()) {
                    std::copy(sign_configs_xy.ab().at(i).begin(),
                            sign_configs_xy.ab().at(i).end(),
                            sign_configs.at(idx).begin());
                }
                if (!sign_configs_xy.bc().empty()) {
                    int shift = n_connecting_loops.ab();
                    std::copy(sign_configs_xy.bc().at(j).begin(),
                            sign_configs_xy.bc().at(j).end(),
                            sign_configs.at(idx).begin() + shift);
                }
                if (!sign_configs_xy.ca().empty()) {
                    int shift = n_connecting_loops.ab() +
                        n_connecting_loops.bc();
                    std::copy(sign_configs_xy.ca().at(k).begin(),
                              sign_configs_xy.ca().at(k).end(),
                              sign_configs.at(idx).begin() + shift);
                }
                idx++;
            } while (++k < sign_configs_xy.ca().size());
        } while (++j < sign_configs_xy.bc().size());
    } while (++i < sign_configs_xy.ab().size());
}



void BiSpectrumDiagram::kernel_arguments(
        size_t rearr_idx,       /* Rearrangement index            */
        size_t sign_idx,        /* Sign config index              */
        size_t overall_loop_idx /* Overall loop assosiation index */
        )
{
    size_t n_coeffs = loop_params.n_coeffs();
    Vec1D<int>& rearrangement = rearrangements.at(rearr_idx);

    size_t rearr_counter = 0; /* How many loop momenta have been assigned */
    size_t k_a_idx = n_coeffs - 1;
    size_t k_b_idx = n_coeffs - 2;

    Vec1D<bool>& sign_config = sign_configs.at(sign_idx);
    Triple<Vec1D<bool>> signs(
        Vec1D<bool>(sign_config.begin(),
                    sign_config.begin() + n_connecting_loops.ab()),
        Vec1D<bool>(sign_config.begin() + n_connecting_loops.ab(),
                    sign_config.begin() + n_connecting_loops.ab() +
                        n_connecting_loops.bc()),
        Vec1D<bool>(sign_config.begin() + n_connecting_loops.ab() +
                        n_connecting_loops.bc(),
                    sign_config.end())
    );

    Triple<Vec1D<int>&> args(
        arg_configs.at(rearr_idx).at(sign_idx).at(overall_loop_idx).a().args,
        arg_configs.at(rearr_idx).at(sign_idx).at(overall_loop_idx).b().args,
        arg_configs.at(rearr_idx).at(sign_idx).at(overall_loop_idx).c().args
    );

    /* How many arguments have been assigned to the kernels a,b,c? */
    Triple<size_t> args_idx(0,0,0);

    /* Determine labels of "main" connecting lines (with external/overall loop
     * momenta), and their complementary with oppisite signs */
    Triple<Vec1D<int>> config_xy(
        Vec1D<int>(n_coeffs, 0),
        Vec1D<int>(n_coeffs, 0),
        Vec1D<int>(n_coeffs, 0)
    );
    Triple<Vec1D<int>> config_xy_sign_flip = config_xy;

    /* Single loop config */
    Vec1D<int> config_single(n_coeffs, 0);

    if (overall_loop_) {
        /* Add (rearranged) first loop momenta to all connecting lines. If
         * overall_loop_idx = 0,1,2 add +Q, and if overall_loop_idx = 3,4,5 add
         * -Q */
        size_t loop_idx = static_cast<size_t>(rearrangement.at(rearr_counter++));
        if (overall_loop_idx < 3) {
            config_xy.ab().at(loop_idx) = 1;
            config_xy.bc().at(loop_idx) = 1;
            config_xy.ca().at(loop_idx) = 1;
        }
        else {
            config_xy.ab().at(loop_idx) = -1;
            config_xy.bc().at(loop_idx) = -1;
            config_xy.ca().at(loop_idx) = -1;
        }
    }
    if (n_ab == 0 || (overall_loop_ && (overall_loop_idx % 3 == 0))) {
        config_xy.bc().at(k_b_idx) = -1; /* q_bc += -k_b */
        config_xy.ca().at(k_a_idx) = 1;  /* q_bc += +k_a */
    }
    else if (n_bc == 0 || (overall_loop_ && (overall_loop_idx % 3 == 1))) {
        config_xy.ab().at(k_b_idx) = 1;  /* q_ab += k_b */
        config_xy.ca().at(k_a_idx) = 1;  /* q_ca += -k_c = k_a + k_b */
        config_xy.ca().at(k_b_idx) = 1;
    }
    else if (n_ca == 0 || (overall_loop_ && (overall_loop_idx % 3 == 2))) {
        config_xy.ab().at(k_a_idx) = -1;  /* q_ab += -k_a */
        config_xy.bc().at(k_a_idx) = -1;  /* q_ca += k_c = - k_a - k_b */
        config_xy.bc().at(k_b_idx) = -1;
    }
    else {
        throw(std::logic_error(
            "BiSpectrumDiagram::kernel_arguments(): Neither n_ab = 0, n_bc = 0, "
            "n_ca = 0, nor valid overall_loop_idx."));
    }

    /* Cache labels corresponding to q_ab, q_bc and q_ca (i.e. total momentum
     * of all connecting loops) for this rearrangement */
    q_xy_labels.at(rearr_idx).at(overall_loop_idx).ab() =
        config2label(config_xy.ab());
    q_xy_labels.at(rearr_idx).at(overall_loop_idx).bc() =
        config2label(config_xy.bc());
    q_xy_labels.at(rearr_idx).at(overall_loop_idx).ca() =
        config2label(config_xy.ca());

    /* Add connecting loops */
    if (n_ab > 0) {
        for (int n = 2; n <= n_ab; ++n) {
            size_t loop_idx = static_cast<size_t>(rearrangement.at(rearr_counter++));

            /* Opposite sign for main connecting loop */
            config_xy.ab().at(loop_idx) = (signs.ab().at(static_cast<size_t>(n - 2))) ? -1 : 1;
            config_single.at(loop_idx)  = (signs.ab().at(static_cast<size_t>(n - 2))) ? 1 : -1;

            args.b().at(args_idx.b()++) = config2label(config_single);
            config_single.at(loop_idx) *= -1;
            args.a().at(args_idx.a()++) = config2label(config_single);

            /* Reset config_single */
            config_single.at(loop_idx) *= 0;
        }

        int label = config2label(config_xy.ab());
        args.b().at(args_idx.b()++) = label;
        args.a().at(args_idx.a()++) = flip_signs(label, n_coeffs);
    }

    if (n_bc > 0) {
        for (int n = 2; n <= n_bc; ++n) {
            size_t loop_idx = static_cast<size_t>(rearrangement.at(rearr_counter++));

            /* Opposite sign for main connecting loop */
            config_xy.bc()[loop_idx] = (signs.bc().at(static_cast<size_t>(n - 2))) ? -1 : 1;
            config_single[loop_idx]  = (signs.bc().at(static_cast<size_t>(n - 2))) ? 1 : -1;

            args.c().at(args_idx.c()++) = config2label(config_single);
            config_single[loop_idx] *= -1;
            args.b().at(args_idx.b()++) = config2label(config_single);

            /* Reset config_single */
            config_single.at(loop_idx) *= 0;
        }

        int label = config2label(config_xy.bc());
        args.c().at(args_idx.c()++) = label;
        args.b().at(args_idx.b()++) = flip_signs(label, n_coeffs);
    }

    if (n_ca > 0) {
        for (int n = 2; n <= n_ca; ++n) {
            size_t loop_idx = static_cast<size_t>(rearrangement.at(rearr_counter++));

            /* Opposite sign for main connecting loop */
            config_xy.ca()[loop_idx] = (signs.ca().at(static_cast<size_t>(n - 2))) ? -1 : 1;
            config_single[loop_idx]  = (signs.ca().at(static_cast<size_t>(n - 2))) ? 1 : -1;

            args.a().at(args_idx.a()++) = config2label(config_single);
            config_single.at(loop_idx) *= -1;
            args.c().at(args_idx.c()++) = config2label(config_single);

            /* Reset config_single */
            config_single.at(loop_idx) *= 0;
        }

        int label = config2label(config_xy.ca());
        args.a().at(args_idx.a()++) = label;
        args.c().at(args_idx.c()++) = flip_signs(label, n_coeffs);
    }

    /* Cache labels corresponding to q_ab1, q_bc1 and q_ca1 for this
     * rearrangement, sign config and overall_loop setup */
    q_xy1_labels.at(rearr_idx).at(sign_idx).at(overall_loop_idx).ab() =
        config2label(config_xy.ab());
    q_xy1_labels.at(rearr_idx).at(sign_idx).at(overall_loop_idx).bc() =
        config2label(config_xy.bc());
    q_xy1_labels.at(rearr_idx).at(sign_idx).at(overall_loop_idx).ca() =
        config2label(config_xy.ca());

    /* Add self loops */
    for (int i = 0; i < n_a; ++i) {
        size_t loop_idx = static_cast<size_t>(rearrangement.at(rearr_counter++));
        config_single.at(loop_idx) = 1;
        args.a().at(args_idx.a()++) = config2label(config_single);
        config_single.at(loop_idx) = -1;
        args.a().at(args_idx.a()++) = config2label(config_single);
        config_single.at(loop_idx) = 0;
    }
    for (int i = 0; i < n_b; ++i) {
        size_t loop_idx = static_cast<size_t>(rearrangement.at(rearr_counter++));
        config_single.at(loop_idx) = 1;
        args.b().at(args_idx.b()++) = config2label(config_single);
        config_single.at(loop_idx) = -1;
        args.b().at(args_idx.b()++) = config2label(config_single);
        config_single.at(loop_idx) = 0;
    }
    for (int i = 0; i < n_c; ++i) {
        size_t loop_idx = static_cast<size_t>(rearrangement.at(rearr_counter++));
        config_single.at(loop_idx) = 1;
        args.c().at(args_idx.c()++) = config2label(config_single);
        config_single.at(loop_idx) = -1;
        args.c().at(args_idx.c()++) = config2label(config_single);
        config_single.at(loop_idx) = 0;
    }

    // Fill remaining spots with zero-label
    int zero_label = loop_params.zero_label();
    size_t n_kernel_args = loop_params.n_kernel_args();
    while (args_idx.a() < n_kernel_args) {
        args.a().at(args_idx.a()++) = zero_label;
    }
    while (args_idx.b() < n_kernel_args) {
        args.b().at(args_idx.b()++) = zero_label;
    }
    while (args_idx.c() < n_kernel_args) {
        args.c().at(args_idx.c()++) = zero_label;
    }

#if DEBUG >= 1
    /* arguments_2_kernel_index() assumes that the length of arguments[] is
     * n_kernel_args. Checking this explicitly: */
    if (args.a().size() != n_kernel_args ||
        args.b().size() != n_kernel_args ||
        args.c().size() != n_kernel_args
       ) {
        throw(std::logic_error(
            "BiSpectrumDiagram::kernel_arguments(): Size of left argument "
            "vector does not equal n_kernel_args."));
    }
#endif

    arg_configs.at(rearr_idx).at(sign_idx).at(overall_loop_idx).a().kernel_index
        = loop_params.args_2_kernel_index(args.a().data());
    arg_configs.at(rearr_idx).at(sign_idx).at(overall_loop_idx).b().kernel_index
        = loop_params.args_2_kernel_index(args.b().data());
    arg_configs.at(rearr_idx).at(sign_idx).at(overall_loop_idx).c().kernel_index
        = loop_params.args_2_kernel_index(args.c().data());
}



BiSpectrumDiagram::BiSpectrumDiagram(
        const LoopParameters& loop_params,
        int n_ab, int n_bc, int n_ca,
        int n_a, int n_b, int n_c
        ) :
    loop_params(loop_params),
    n_ab(n_ab), n_bc(n_bc), n_ca(n_ca),
    n_a(n_a), n_b(n_b), n_c(n_c)
{
    int n_loops = loop_params.n_loops();
    size_t n_kernel_args = loop_params.n_kernel_args();

    if (n_ab + n_bc + n_ca + n_a + n_b + n_c != n_loops + 2) {
        throw(std::invalid_argument(
            "BiSpectrumDiagram::BiSpectrumDiagram(): n_ab + n_bc + n_ca + n_a "
            "+ n_b + n_c != n_loops + 2"));
    }
    if ((n_ab > 0) + (n_bc > 0) + (n_ca > 0) < 2) {
        throw(std::invalid_argument(
            "BiSpectrumDiagram::BiSpectrumDiagram(): More than one number out "
            "of {n_ab, n_bc, n_ca} is zero."));
    }
    overall_loop_ = true;
    if (n_ab == 0 || n_bc == 0 || n_ca == 0) {
        overall_loop_ = false;
    }

    diagram_factor_ = static_cast<int>(
        gsl_sf_fact(2 * n_a + n_ab + n_ca) *
        gsl_sf_fact(2 * n_b + n_bc + n_ab) *
        gsl_sf_fact(2 * n_c + n_ca + n_bc) /
        (
         pow(2,n_a + n_b + n_c) *
         gsl_sf_fact(n_a) * gsl_sf_fact(n_b) * gsl_sf_fact(n_c) *
         gsl_sf_fact(n_ab) * gsl_sf_fact(n_bc) * gsl_sf_fact(n_ca)
         ));

    rearrangements.resize(static_cast<size_t>(
            gsl_sf_fact(n_loops) /
        (
         gsl_sf_fact(n_ab > 0 ? n_ab - 1 : 0) *
         gsl_sf_fact(n_bc > 0 ? n_bc - 1 : 0) *
         gsl_sf_fact(n_ca > 0 ? n_ca - 1 : 0) *
         gsl_sf_fact(n_a) * gsl_sf_fact(n_b) * gsl_sf_fact(n_c)
        )));

    compute_rearrangements(n_loops);

    n_connecting_loops.ab() = n_ab > 1 ? n_ab - 1 : 0;
    n_connecting_loops.bc() = n_bc > 1 ? n_bc - 1 : 0;
    n_connecting_loops.ca() = n_ca > 1 ? n_ca - 1 : 0;

    Triple<Vec2D<bool>> sign_configs_xy(
        connecting_line_sign_flips(n_connecting_loops.ab()),
        connecting_line_sign_flips(n_connecting_loops.bc()),
        connecting_line_sign_flips(n_connecting_loops.ca())
    );

    compute_sign_flips(sign_configs_xy);

    /* Allocate memory */
    arg_configs.resize(rearrangements.size());
    q_xy1_labels.resize(rearrangements.size());
    q_xy_labels.resize(rearrangements.size());

    for (size_t i = 0; i < rearrangements.size(); ++i) {
        arg_configs.at(i).resize(sign_configs.size());
        q_xy1_labels.at(i).resize(sign_configs.size());

        for (size_t j = 0; j < sign_configs.size(); ++j) {
            /* There are three ways to associate the loop momenta to each
             * connecting lines, and 2 sign flips for each, hence 6 total */
            size_t overall_loop_assosiations = overall_loop_ ? 6 : 1;
            arg_configs.at(i).at(j).resize(overall_loop_assosiations);
            q_xy1_labels.at(i).at(j).resize(overall_loop_assosiations);

            q_xy_labels.at(i).resize(overall_loop_assosiations);

            for (size_t k = 0; k < overall_loop_assosiations; ++k) {
                arg_configs.at(i).at(j).at(k).a().args.resize(n_kernel_args);
                arg_configs.at(i).at(j).at(k).b().args.resize(n_kernel_args);
                arg_configs.at(i).at(j).at(k).c().args.resize(n_kernel_args);
                /* Initialize arguments and kernel indices for this
                 * configuration */
                kernel_arguments(i, j, k);
            }
        }
    }
}



void BiSpectrumDiagram::connecting_lines_factors(
        size_t rearr_idx,
        size_t sign_idx,
        size_t overall_loop_idx,
        const Vec1D<double>& loop_magnitudes,
        const Strided2DVec<double>& dot_products,
        Triple<double>& q_xy1, /* out */
        int& heaviside_theta   /* out */
        ) const
{
    heaviside_theta = 1;

    /* q_ab1 is the momentum of the connecting line ab whose loop momentum
     * is evaluated by the momentum conserving delta function, similarly
     * for q_bc1 and q_ca1.
     * q_ab (or q_bc, q_ca) is the *total* momentum flowing in the connection
     * ab (or bc or ca) between blobs */
    Triple<size_t> q_xy1_label(
        static_cast<size_t>(
            q_xy1_labels.at(rearr_idx).at(sign_idx).at(overall_loop_idx).ab()),
        static_cast<size_t>(
            q_xy1_labels.at(rearr_idx).at(sign_idx).at(overall_loop_idx).bc()),
        static_cast<size_t>(
            q_xy1_labels.at(rearr_idx).at(sign_idx).at(overall_loop_idx).ca())
    );

    q_xy1 = {
        std::sqrt(dot_products(q_xy1_label.ab(), q_xy1_label.ab())),
        std::sqrt(dot_products(q_xy1_label.bc(), q_xy1_label.bc())),
        std::sqrt(dot_products(q_xy1_label.ca(), q_xy1_label.ca()))
    };

    /* How many loop momenta have already been assigned? */
    size_t offset = 0;

    if (overall_loop_) {
        offset = 1;
        /* If the total connecting line momentum with which the overall loop
         * momentum is associated (according to overall_loop_idx) is larger
         * than any of the other total connecting line momenta, the diagram
         * configuration is skipped (heaviside_theta = 0) */
        Triple<size_t> q_xy_label(
            static_cast<size_t>(
                q_xy_labels.at(rearr_idx).at(overall_loop_idx).ab()),
            static_cast<size_t>(
                q_xy_labels.at(rearr_idx).at(overall_loop_idx).bc()),
            static_cast<size_t>(
                q_xy_labels.at(rearr_idx).at(overall_loop_idx).ca())
        );

        Triple<double> q_xy(
            std::sqrt(dot_products(q_xy_label.ab(), q_xy_label.ab())),
            std::sqrt(dot_products(q_xy_label.bc(), q_xy_label.bc())),
            std::sqrt(dot_products(q_xy_label.ca(), q_xy_label.ca()))
        );

        switch (overall_loop_idx % 3) {
            case 0:
                if (q_xy.ab() >= q_xy.bc() || q_xy.ab() >= q_xy.ca()) {
                    heaviside_theta = 0;
                    return;
                }
                break;
            case 1:
                if (q_xy.bc() >= q_xy.ab() || q_xy.bc() >= q_xy.ca()) {
                    heaviside_theta = 0;
                    return;
                }
                break;
            case 2:
                if (q_xy.ca() >= q_xy.ab() || q_xy.ca() >= q_xy.bc()) {
                    heaviside_theta = 0;
                    return;
                }
                break;
        }
    }

    const Vec1D<int>& rearr = rearrangements.at(rearr_idx);

    /* Heaviside-theta (q_ab1 - Q(rearr(offset)) */
    if (n_ab > 1) {
        double q = loop_magnitudes.at(static_cast<size_t>(rearr.at(offset)));
        heaviside_theta *= (q_xy1.ab() > q ? n_ab : 0);
        offset += static_cast<size_t>(n_connecting_loops.ab());
    }
#if DEBUG >= 1
    /* Check that the heaviside-theta (Q2(rearr) - Q3(rearr)) etc. are
     * satisfied by (reparametrized) momenta from CUBA */
    for (size_t i = 1; i < static_cast<size_t>(n_connecting_loops.ab()); ++i) {
      double q_ab_i = loop_magnitudes.at(static_cast<size_t>(rearr.at(offset + i)));
      double q_ab_j = loop_magnitudes.at(static_cast<size_t>(rearr.at(offset + i + 1)));
      if (q_ab_i < q_ab_j) {
          throw(std::logic_error("Heaviside theta(ab): Q" +
                      std::to_string(rearr.at(offset + i) + 1) + " < Q" +
                      std::to_string(rearr.at(offset + i + 1) + 1) + "."));
        }
    }
#endif

    /* Same for n_bc */
    if (n_bc > 1) {
        double q = loop_magnitudes.at(static_cast<size_t>(rearr.at(offset)));
        heaviside_theta *= (q_xy1.bc() > q ? n_bc : 0);
        offset += static_cast<size_t>(n_connecting_loops.bc());
    }
#if DEBUG >= 1
    /* Check that the heaviside-theta (Q2(rearr) - Q3(rearr)) etc. are
     * satisfied by (reparametrized) momenta from CUBA */
    for (size_t i = 1; i < static_cast<size_t>(n_connecting_loops.bc()); ++i) {
      double q_bc_i = loop_magnitudes.at(static_cast<size_t>(rearr.at(offset + i)));
      double q_bc_j = loop_magnitudes.at(static_cast<size_t>(rearr.at(offset + i + 1)));
      if (q_bc_i < q_bc_j) {
          throw(std::logic_error("Heaviside theta(bc): Q" +
                      std::to_string(rearr.at(offset + i) + 1) + " < Q" +
                      std::to_string(rearr.at(offset + i + 1) + 1) + "."));
        }
    }
#endif

    /* Same for n_ca */
    if (n_ca > 1) {
        double q = loop_magnitudes.at(static_cast<size_t>(rearr.at(offset)));
        heaviside_theta *= (q_xy1.ca() > q ? n_ca : 0);
        offset += static_cast<size_t>(n_connecting_loops.ca());
    }
#if DEBUG >= 1
    /* Check that the heaviside-theta (Q2(rearr) - Q3(rearr)) etc. are
     * satisfied by (reparametrized) momenta from CUBA */
    for (size_t i = 1; i < static_cast<size_t>(n_connecting_loops.ca()); ++i) {
      double q_ca_i = loop_magnitudes.at(static_cast<size_t>(rearr.at(offset + i)));
      double q_ca_j = loop_magnitudes.at(static_cast<size_t>(rearr.at(offset + i + 1)));
      if (q_ca_i < q_ca_j) {
          throw(std::logic_error("Heaviside theta(ca): Q" +
                      std::to_string(rearr.at(offset + i) + 1) + " < Q" +
                      std::to_string(rearr.at(offset + i + 1) + 1) + "."));
        }
    }
#endif
}



std::string BiSpectrumDiagram::tags() const
{
    std::ostringstream oss;
    oss << "(n_a, n_b, n_c, n_ab, n_bc, n_ca) = ("
        << n_a << ", "
        << n_b << ", "
        << n_c << ", "
        << n_ab << ", "
        << n_bc << ", "
        << n_ca << ")";
    return oss.str();
}



std::string BiSpectrumDiagram::argument_configuration(
        size_t rearr_idx,
        size_t sign_idx,
        size_t overall_loop_idx
        ) const
{
    std::ostringstream oss;
    oss << labels2string(
        arg_configs.at(rearr_idx).at(sign_idx).at(overall_loop_idx).a().args,
        loop_params.n_coeffs(), loop_params.spectrum())
        << " "
        << labels2string(
            arg_configs.at(rearr_idx).at(sign_idx).at(overall_loop_idx).b().args,
            loop_params.n_coeffs(), loop_params.spectrum())
        << " "
        << labels2string(
            arg_configs.at(rearr_idx).at(sign_idx).at(overall_loop_idx).c().args,
            loop_params.n_coeffs(), loop_params.spectrum());
    return oss.str();
}



std::ostream& operator<<(std::ostream& out, const PowerSpectrumDiagram& diagram)
{
    out << diagram.tags() << std::endl;

    for (size_t a = 0; a < diagram.n_rearrangements(); ++a) {
        for (size_t b = 0; b < diagram.n_sign_configs(); ++b) {
            out << "\t" << diagram.argument_configuration(a, b) << std::endl;
        }
    }
    return out;
}



std::ostream& operator<<(std::ostream& out, const BiSpectrumDiagram& diagram)
{
    out << diagram.tags() << std::endl;

    for (size_t i = 0; i < diagram.n_rearrangements(); ++i) {
        for (size_t j = 0; j < diagram.n_sign_configs(); ++j) {
            if (diagram.overall_loop()) {
                for (size_t k = 0; k < 6; ++k) {
                    out << "\t"
                        << diagram.argument_configuration(i, j, k)
                        << std::endl;
                }
            }
            else {
                out << "\t"
                    << diagram.argument_configuration(i, j)
                    << std::endl;
            }
        }
    }
    return out;
}



Vec1D<PowerSpectrumDiagram> ps::construct_diagrams(const LoopParameters& params) {
    Vec1D<PowerSpectrumDiagram> diagrams;

    int n_loops = params.n_loops();
    if (n_loops < 1) {
        throw std::logic_error(
            "ps::construct_diagrams(): called with n_loops < 1.");
    }

    int m = 0;
    int index = 0;

    // Find (distinct) power spectrum diagrams at L-loop
    // They satisfy: m >= 1; l,r > 0; l + r + m = L + 1
    for (m = 1; m <= n_loops + 1; ++m) {
        int l = n_loops + 1 - m;
        int r = 0;
        while (l >= r) {
            if (index >= 2 * n_loops) {
                throw(std::logic_error(
                    "construct_diagrams(): Index larger than 2 * n_loops."));
            }
            diagrams.emplace_back(params, m, l, r);

            l = n_loops + 1 - m - (++r);
            index++;
        };
    }
    return diagrams;
}



Vec1D<BiSpectrumDiagram> bs::construct_diagrams(const LoopParameters& loop_params)
{
    Vec1D<BiSpectrumDiagram> diagrams;

    int n_loops = loop_params.n_loops();
    if (n_loops < 1) {
        throw std::logic_error(
            "bs::construct_diagrams(): called with n_loops < 1.");
    }

    /* Need to place L + 2 lines */
    /* First, consider diagram where a connecting line is zero, e.g. n_ab = 0 */
    for (int n_bc = 1; n_bc <= n_loops + 1; ++n_bc) {
        for (int n_ca = 1; n_ca <= n_loops - n_bc + 2; ++n_ca) {
            /* Any remaining loops are self-loops */
            int n_self = n_loops + 2 - n_bc - n_ca;
            /* Assining loops to n_a */
            for (int n_a = 0; n_a <= n_self; ++n_a) {
                /* Split remaining loops between n_b and n_c */
                int n_b = n_self - n_a;
                int n_c = 0;
                while (n_b >= n_c) {
                    /* Three diagrams, picking n_ab, n_bc, or n_ca as line being zero */
                    diagrams.emplace_back(loop_params, 0,    n_bc, n_ca, n_a, n_b, n_c);
                    diagrams.emplace_back(loop_params, n_bc, 0,    n_ca, n_a, n_b, n_c);
                    diagrams.emplace_back(loop_params, n_ca, n_bc, 0,    n_a, n_b, n_c);
                    /* If n_b != n_c, there are diagrams with n_b <-> n_c */
                    if (n_b != n_c) {
                        diagrams.emplace_back(loop_params, 0,    n_bc, n_ca, n_a, n_c, n_b);
                        diagrams.emplace_back(loop_params, n_bc, 0,    n_ca, n_a, n_c, n_b);
                        diagrams.emplace_back(loop_params, n_ca, n_bc, 0,    n_a, n_c, n_b);
                    }

                    n_b = n_self - n_a - (++n_c);
                }
            }
        }
    }

    /* Diagrams with n_ab, n_bc, n_ca > 0 */
    for (int n_ab = 1; n_ab <= n_loops + 1; ++n_ab) {
        for (int n_bc = 1; n_bc <= n_loops - n_ab + 2; ++n_bc) {
            for (int n_ca = 1; n_ca <= n_loops - n_ab - n_bc + 2; ++n_ca) {
                /* Any remaining loops are self-loops */
                int n_self = n_loops + 2 - n_ab - n_bc - n_ca;
                /* Assining loops to n_a */
                for (int n_a = 0; n_a <= n_self; ++n_a) {
                    /* Split remaining loops between n_b and n_c */
                    int n_b = n_self - n_a;
                    int n_c = 0;
                    while (n_b >= n_c) {
                        diagrams.emplace_back(loop_params, n_ab, n_bc, n_ca, n_a, n_b, n_c);
                        /* If n_b != n_c, there is a diagram with n_b <-> n_c */
                        if (n_b != n_c) {
                            diagrams.emplace_back(loop_params, n_ab, n_bc, n_ca, n_a, n_c, n_b);
                        }

                        n_b = n_self - n_a - (++n_c);
                    }
                }
            }
        }
    }

    return diagrams;
}
