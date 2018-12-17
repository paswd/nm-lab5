#include <cmath>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "types.hpp"
#include "parabolic_pde.hpp"

using namespace std;

static constexpr ldbl phi_0_t(const ldbl &, const ldbl &, const ldbl &,
    const ldbl &) noexcept;
static constexpr ldbl phi_l_t(const ldbl &, const ldbl &, const ldbl &,
    const ldbl &) noexcept;
static constexpr ldbl psi_x(const ldbl &, const ldbl &, const ldbl &,
    const ldbl &) noexcept;
static constexpr ldbl u_exact(const ldbl &, const ldbl &, const ldbl &,
    const ldbl &, const ldbl &) noexcept;

static constexpr ldbl L = PI_LDBL;
static constexpr ldbl ALPHA = 1.0,
                      BETA  = 1.0,
                      GAMMA = 1.0,
                      DELTA = 1.0;

int main(int argc, const char *argv[])
{
    try
    {
        ldbl a, b, c, t;
        size_t n, k;
        uint boundary;
        cin >> a >> b >> c >> t >> n >> k >> boundary;

        const ldbl h = L / n;
        const ldbl tau = t / k;

        vector<ldbl> explicit_fdm_error,
            implicit_fdm_error,
            crank_nicolson_error;

        const ublas::vector<ldbl>
            explicit_fdm_u = explicit_fdm<ldbl>(a, b, c,
                L, n, t, k, ALPHA, BETA, GAMMA, DELTA,
                phi_0_t, phi_l_t, psi_x, static_cast<NumDiff>(boundary),
                [&] (const ublas::vector<ldbl> &u_k) -> void
                {
                    ldbl error = 0.0;
                    const size_t size = u_k.size(),
                        j = explicit_fdm_error.size();
                    for (size_t i = 0; i < size; ++i)
                    {
                        error = max(error, abs(u_k[i] - u_exact(a, b, c,
                            i * h, j * tau)));
                    }
                    explicit_fdm_error.push_back(error);
                }),
            implicit_fdm_u = implicit_fdm<ldbl>(a, b, c,
                L, n, t, k, ALPHA, BETA, GAMMA, DELTA,
                phi_0_t, phi_l_t, psi_x, static_cast<NumDiff>(boundary),
                [&] (const ublas::vector<ldbl> &u_k) -> void
                {
                    ldbl error = 0.0;
                    const size_t size = u_k.size(),
                        j = implicit_fdm_error.size();
                    for (size_t i = 0; i < size; ++i)
                    {
                        error = max(error, abs(u_k[i] - u_exact(a, b, c,
                            i * h, j * tau)));
                    }
                    implicit_fdm_error.push_back(error);
                }),
            crank_nicolson_u = crank_nicolson<ldbl>(a, b, c,
                L, n, t, k, ALPHA, BETA, GAMMA, DELTA,
                phi_0_t, phi_l_t, psi_x, static_cast<NumDiff>(boundary),
                [&] (const ublas::vector<ldbl> &u_k) -> void
                {
                    ldbl error = 0.0;
                    const size_t size = u_k.size(),
                        j = crank_nicolson_error.size();
                    for (size_t i = 0; i < size; ++i)
                    {
                        error = max(error, abs(u_k[i] - u_exact(a, b, c,
                            i * h, j * tau)));
                    }
                    crank_nicolson_error.push_back(error);
                });

        for (size_t i = 0; i < explicit_fdm_u.size(); ++i)
        {
            cout << fixed << setprecision(18) << i * h << ',' <<
                u_exact(a, b, c, i * h, t) << ',' << explicit_fdm_u[i] << ',' <<
                implicit_fdm_u[i] << ',' << crank_nicolson_u[i] << '\n';
        }
        if (argc > 1)
        {
            fstream error_stream(argv[1], ios_base::out | ios_base::trunc);
            for (size_t i = 0; i < explicit_fdm_error.size(); ++i)
            {
                error_stream << fixed << setprecision(18) << i * tau << ',' <<
                    explicit_fdm_error[i] << ',' << implicit_fdm_error[i] << ',' <<
                    crank_nicolson_error[i] << '\n';
            }
        }
    } catch (const exception &ex)
    {
        cout << ex.what() << '\n';
    }

    return 0;
}

static constexpr ldbl phi_0_t(const ldbl &a, const ldbl &b, const ldbl &c,
    const ldbl &t) noexcept
{
    return exp((c - a) * t);
}

static constexpr ldbl phi_l_t(const ldbl &a, const ldbl &b, const ldbl &c,
    const ldbl &t) noexcept
{
    return -exp((c - a) * t);
}

static constexpr ldbl psi_x(const ldbl &, const ldbl &, const ldbl &,
    const ldbl &x) noexcept
{
    return sin(x);
}

static constexpr ldbl u_exact(const ldbl &a, const ldbl &b, const ldbl &c,
    const ldbl &x, const ldbl &t) noexcept
{
    return exp((c - a) * t) * sin(x);
}
