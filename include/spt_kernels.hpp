#ifndef SPT_KERNELS_HPP
#define SPT_KERNELS_HPP

class IntegrandTables;

int compute_SPT_kernels(
        const int arguments[],
        int kernel_index,
        int n,
        IntegrandTables& tables
        );

/* For kernel computers: Validate that number of non-zero-label arguments are
 * in fact n */
void kernel_computer_validate_n(
        const int arguments[],
        int n,
        IntegrandTables& tables
        );

/* For kernel computers: Validate that kernel_index is indeed the appropriate
 * kernel index for the given arguments */
void kernel_computer_validate_kernel_index(
        const int arguments[],
        int kernel_index,
        IntegrandTables& tables
        );

#endif /* ifndef SPT_KERNELS_HPP */
