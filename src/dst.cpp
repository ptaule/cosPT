
void discrete_sine_transform(
        const Vec1D<double>& input_real,
        Vec1D<double>& output
        )
{
    size_t N = input_real.size();

    /* Following algorithm from CLASS-PT http://arxiv.org/abs/2004.10607 */

    Vec1D<double> tmp_arr(2*N,0);

    for (size_t i = 0; i < N; ++i) {
        tmp_arr[2*i + 1] = pow(-1, i) * input_real.at(i);
    }

    Vec1D<double> dft_real = tmp_arr;

    std::reverse(tmp_arr.begin(), tmp_arr.end());
    dft_real.insert(dft_real.end(), tmp_arr.begin(), tmp_arr.end());

    gsl_fft_real_radix2_transform(dft_real.data(), 1, 4*N);

    /* Real coefficients from DFT are the first half, moreover, we care about the first N */
    output.assign(dft_real.begin(), dft_real.begin() + N);
    std::reverse(output.begin(), output.end());
}
