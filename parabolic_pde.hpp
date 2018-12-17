#ifndef __PARABOLIC_PDE_HPP__
#define __PARABOLIC_PDE_HPP__

#include <cstddef>

#include <array>
#include <limits>
#include <type_traits>

#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include "thomas_algorithm.hpp"

enum NumDiff : uint
{
    TWO_POINT_FIRST_ORDER,
    TWO_POINT_SECOND_ORDER,
    THREE_POINT_SECOND_ORDER,
};

static constexpr std::size_t TWO = 2U;

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
ublas::vector<T> explicit_fdm(const T, const T, const T,
    const T, std::size_t, const T, std::size_t,
    const T, const T, const T, const T,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    NumDiff,
    const std::function<void (const ublas::vector<T> &)> &get_error);

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
ublas::vector<T> implicit_fdm(const T, const T, const T,
    const T, std::size_t, const T, std::size_t,
    const T, const T, const T, const T,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    NumDiff,
    const std::function<void (const ublas::vector<T> &)> &);

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
ublas::vector<T> crank_nicolson(const T, const T, const T,
    const T, std::size_t, const T, std::size_t,
    const T, const T, const T, const T,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    const std::function<T (const T &, const T &, const T &, const T &)> &,
    NumDiff,
    const std::function<void (const ublas::vector<T> &)> &);

static bool is_num_diff_includes(NumDiff) noexcept;

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
static void to_three_diagonal(ublas::vector<T> &, ublas::vector<T> &, ublas::vector<T> &,
    ublas::vector<T> &, const T, const T) noexcept;

template<typename T, typename>
ublas::vector<T> explicit_fdm(const T a, const T b, const T c,
    const T l, const std::size_t n_upper,
    const T t, const std::size_t k_upper,
    const T alpha, const T beta,
    const T gamma, const T delta,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_0_t,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_l_t,
    const std::function<T (const T &, const T &, const T &, const T &)> &psi_x,
    const NumDiff boundary,
    const std::function<void (const ublas::vector<T> &)> &get_error)
{
    static constexpr T EPSILON = std::numeric_limits<T>::epsilon();
    const T h = l / n_upper, tau = t / k_upper, sigma = a * a * tau / (h * h),
        b_0 = alpha * (2.0 * a * a / h + h / tau - c * h) -
            beta * (2.0 * a * a - b * h),
        c_0 = -2.0 * alpha * a * a / h,
        a_n = -2.0 * gamma * a * a / h,
        b_n = gamma * (2.0 * a * a / h + h / tau - c * h) +
            delta * (2.0 * a * a + b * h);
    if (n_upper < 2U || !k_upper || h < EPSILON || tau < EPSILON || sigma > 0.5 ||
        !is_num_diff_includes(boundary))
    {
        throw std::logic_error("n < 2 || !k || h < epsilon || tau < epsilon || sigma > 0.5");
    }

    std::array<ublas::vector<T>, TWO> w_h_tau;
    w_h_tau.fill(ublas::vector<T>(n_upper + 1));
    for (std::size_t j = 0; j <= n_upper; ++j)
    {
        w_h_tau[0][j] = psi_x(a, b, c, j * h);
    }
    get_error(w_h_tau[0]);

    for (std::size_t k = 1; k <= k_upper; ++k)
    {
        ublas::vector<T> &u_k = w_h_tau[k % TWO],
            &u_k_minus_1 = w_h_tau[(k - 1) % TWO];
        for (std::size_t j = 1; j < n_upper; ++j)
        {
            u_k[j] = (sigma + b * tau / (2.0 * h)) * u_k_minus_1[j + 1] +
                (1.0 - 2.0 * sigma + c * tau) * u_k_minus_1[j] +
                (sigma - b * tau / (2.0 * h)) * u_k_minus_1[j - 1];
        }
        switch (boundary)
        {
            default:
            case TWO_POINT_FIRST_ORDER:
                u_k[0] =
                    (-alpha * u_k[1] + phi_0_t(a, b, c, k * tau) * h) /
                    (beta * h - alpha);
                u_k[n_upper] =
                    (gamma * u_k[n_upper - 1] + phi_l_t(a, b, c, k * tau) * h) /
                    (delta * h + gamma);
                break;
            case TWO_POINT_SECOND_ORDER:
                u_k[0] = (alpha * h / tau * u_k_minus_1[0] -
                    phi_0_t(a, b, c, k * tau) * (2.0 * a * a - b * h) -
                    c_0 * u_k[1]) / b_0;
                u_k[n_upper] = (gamma * h / tau * u_k_minus_1[n_upper] +
                    phi_l_t(a, b, c, k * tau) * (2.0 * a * a + b * h) -
                    a_n * u_k[n_upper - 1]) / b_n;
                break;
            case THREE_POINT_SECOND_ORDER:
                u_k[0] =
                    (-4.0 * alpha * u_k[1] + alpha * u_k[2U] +
                        2.0 * phi_0_t(a, b, c, k * tau) * h)
                    / (2.0 * beta * h - 3.0 * alpha);
                u_k[n_upper] =
                    (4.0 * gamma * u_k[n_upper - 1] - gamma * u_k[n_upper - 2U] +
                        2.0 * phi_l_t(a, b, c, k * tau) * h) /
                    (2.0 * delta * h + 3.0 * gamma);
                break;
        }
        get_error(u_k);
    }

    return w_h_tau[k_upper % TWO];
}

template<typename T, typename>
ublas::vector<T> implicit_fdm(const T a, const T b, const T c,
    const T l, const std::size_t n_upper,
    const T t, const std::size_t k_upper,
    const T alpha, const T beta,
    const T gamma, const T delta,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_0_t,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_l_t,
    const std::function<T (const T &, const T &, const T &, const T &)> &psi_x,
    const NumDiff boundary,
    const std::function<void (const ublas::vector<T> &)> &get_error)
{
    static constexpr T EPSILON = std::numeric_limits<T>::epsilon();
    const T h = l / n_upper, tau = t / k_upper, sigma = a * a * tau / (h * h);
    if (n_upper < 4U || !k_upper || h < EPSILON || tau < EPSILON ||
        !is_num_diff_includes(boundary))
    {
        throw std::logic_error("n < 4 || !k || h < epsilon || tau < epsilon");
    }

    ublas::vector<T> x, d_j(n_upper + 1),
        a_j(n_upper + 1, sigma - b * tau / (2.0 * h)),
        b_j(n_upper + 1, c * tau - 1.0 - 2.0 * sigma),
        c_j(n_upper + 1, sigma + b * tau / (2.0 * h));
    a_j[0] = c_j[n_upper] = 0.0;
    switch (boundary)
    {
        default:
        case TWO_POINT_FIRST_ORDER:
            b_j[0] = beta * h - alpha;
            c_j[0] = alpha;
            a_j[n_upper] = -gamma;
            b_j[n_upper] = gamma + delta * h;
            break;
        case TWO_POINT_SECOND_ORDER:
            b_j[0] = alpha * (2.0 * a * a / h + h / tau - c * h) -
                beta * (2.0 * a * a - b * h);
            c_j[0] = -2.0 * alpha * a * a / h;
            a_j[n_upper] = -2.0 * gamma * a * a / h;
            b_j[n_upper] = gamma * (2.0 * a * a / h + h / tau - c * h) +
                delta * (2.0 * a * a + b * h);
            break;
        case THREE_POINT_SECOND_ORDER:
            break;
    }
    std::array<ublas::vector<T>, TWO> w_h_tau;
    w_h_tau.fill(ublas::vector<T>(n_upper + 1));
    for (std::size_t j = 0; j <= n_upper; ++j)
    {
        w_h_tau[0][j] = psi_x(a, b, c, j * h);
    }
    get_error(w_h_tau[0]);

    for (std::size_t k = 1; k <= k_upper; ++k)
    {
        ublas::vector<T> &u_k = w_h_tau[k % TWO],
            &u_k_minus_1 = w_h_tau[(k - 1) % TWO];
        for (std::size_t j = 1; j < n_upper; ++j)
        {
            d_j[j] = -u_k_minus_1[j];
        }
        switch (boundary)
        {
            default:
            case TWO_POINT_FIRST_ORDER:
                d_j[0] = phi_0_t(a, b, c, k * tau) * h;
                d_j[n_upper] = phi_l_t(a, b, c, k * tau) * h;
                break;
            case TWO_POINT_SECOND_ORDER:
                d_j[0] = alpha * h / tau * u_k_minus_1[0] -
                    phi_0_t(a, b, c, k * tau) * (2.0 * a * a - b * h);
                d_j[n_upper] = gamma * h / tau * u_k_minus_1[n_upper] +
                    phi_l_t(a, b, c, k * tau) * (2.0 * a * a + b * h);
                break;
            case THREE_POINT_SECOND_ORDER:
                b_j[0] = 2.0 * beta * h - 3.0 * alpha;
                c_j[0] = 4.0 * alpha;
                d_j[0] = 2.0 * phi_0_t(a, b, c, k * tau) * h;
                a_j[n_upper] = -4.0 * gamma;
                b_j[n_upper] = 2.0 * delta * h + 3.0 * gamma;
                d_j[n_upper] = 2.0 * phi_l_t(a, b, c, k * tau) * h;
                to_three_diagonal(a_j, b_j, c_j, d_j, -alpha, gamma);
                break;
        }
        u_k = thomas_algorithm(a_j, b_j, c_j, d_j);
        get_error(u_k);
    }

    return w_h_tau[k_upper % TWO];
}

template<typename T, typename>
ublas::vector<T> crank_nicolson(const T a, const T b, const T c,
    const T l, const std::size_t n_upper,
    const T t, const std::size_t k_upper,
    const T alpha, const T beta,
    const T gamma, const T delta,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_0_t,
    const std::function<T (const T &, const T &, const T &, const T &)> &phi_l_t,
    const std::function<T (const T &, const T &, const T &, const T &)> &psi_x,
    const NumDiff boundary,
    const std::function<void (const ublas::vector<T> &)> &get_error)
{
    static constexpr T THETA = 0.5, EPSILON = std::numeric_limits<T>::epsilon();
    const T h = l / n_upper, tau = t / k_upper, sigma = a * a * tau / (h * h);
    if (n_upper < 4U || !k_upper || h < EPSILON || tau < EPSILON ||
        !is_num_diff_includes(boundary))
    {
        throw std::logic_error("n < 4 || !k || h < epsilon || tau < epsilon");
    }

    ublas::vector<T> x, d_j(n_upper + 1),
        a_j(n_upper + 1, THETA * (sigma - b * tau / (2.0 * h))),
        b_j(n_upper + 1, -1.0 + THETA * (c * tau - 2.0 * sigma)),
        c_j(n_upper + 1, THETA * (sigma + b * tau / (2.0 * h)));
    a_j[0] = c_j[n_upper] = 0.0;
    switch (boundary)
    {
        default:
        case TWO_POINT_FIRST_ORDER:
            b_j[0] = beta * h - alpha;
            c_j[0] = alpha;
            a_j[n_upper] = -gamma;
            b_j[n_upper] = gamma + delta * h;
            break;
        case TWO_POINT_SECOND_ORDER:
            b_j[0] = alpha * (2.0 * a * a / h + h / tau - c * h) -
                beta * (2.0 * a * a - b * h);
            c_j[0] = -2.0 * alpha * a * a / h;
            a_j[n_upper] = -2.0 * gamma * a * a / h;
            b_j[n_upper] = gamma * (2.0 * a * a / h + h / tau - c * h) +
                delta * (2.0 * a * a + b * h);
            break;
        case THREE_POINT_SECOND_ORDER:
            break;
    }
    std::array<ublas::vector<T>, TWO> w_h_tau;
    w_h_tau.fill(ublas::vector<T>(n_upper + 1));
    for (std::size_t j = 0; j <= n_upper; ++j)
    {
        w_h_tau[0][j] = psi_x(a, b, c, j * h);
    }
    get_error(w_h_tau[0]);

    for (std::size_t k = 1; k <= k_upper; ++k)
    {
        ublas::vector<T> &u_k = w_h_tau[k % TWO],
            &u_k_minus_1 = w_h_tau[(k - 1) % TWO];
        for (std::size_t j = 1; j < n_upper; ++j)
        {
            d_j[j] = (THETA - 1.0) * (sigma + b * tau / (2.0 * h)) * u_k_minus_1[j + 1] +
                ((THETA - 1.0) * (c * tau - 2.0 * sigma) - 1.0) * u_k_minus_1[j] +
                (THETA - 1.0) * (sigma - b * tau / (2.0 * h)) * u_k_minus_1[j - 1];
        }
        switch (boundary)
        {
            default:
            case TWO_POINT_FIRST_ORDER:
                d_j[0] = phi_0_t(a, b, c, k * tau) * h;
                d_j[n_upper] = phi_l_t(a, b, c, k * tau) * h;
                break;
            case TWO_POINT_SECOND_ORDER:
                d_j[0] = alpha * h / tau * u_k_minus_1[0] -
                    phi_0_t(a, b, c, k * tau) * (2.0 * a * a - b * h);
                d_j[n_upper] = gamma * h / tau * u_k_minus_1[n_upper] +
                    phi_l_t(a, b, c, k * tau) * (2.0 * a * a + b * h);
                break;
            case THREE_POINT_SECOND_ORDER:
                b_j[0] = 2.0 * beta * h - 3.0 * alpha;
                c_j[0] = 4.0 * alpha;
                d_j[0] = 2.0 * phi_0_t(a, b, c, k * tau) * h;
                a_j[n_upper] = -4.0 * gamma;
                b_j[n_upper] = 2.0 * delta * h + 3.0 * gamma;
                d_j[n_upper] = 2.0 * phi_l_t(a, b, c, k * tau) * h;
                to_three_diagonal(a_j, b_j, c_j, d_j, -alpha, gamma);
                break;
        }
        u_k = thomas_algorithm(a_j, b_j, c_j, d_j);
        get_error(u_k);
    }

    return w_h_tau[k_upper % TWO];
}

static bool is_num_diff_includes(const NumDiff boundary) noexcept
{
    switch (boundary)
    {
        case TWO_POINT_FIRST_ORDER:
        case TWO_POINT_SECOND_ORDER:
        case THREE_POINT_SECOND_ORDER:
            return true;
        default:
            return false;
    }
}

template<typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
static void to_three_diagonal(ublas::vector<T> &a, ublas::vector<T> &b, ublas::vector<T> &c,
    ublas::vector<T> &d, const T m_0_2, const T m_n_n_minus_2) noexcept
{
    const std::size_t n = d.size() - 1;
    const T alpha_0 = -m_0_2 / c[1], alpha_n = -m_n_n_minus_2 / a[n - 1];
    b[0] = std::fma(alpha_0, a[1], b[0]);
    c[0] = std::fma(alpha_0, b[1], c[0]);
    d[0] = std::fma(alpha_0, d[1], d[0]);
    a[n] = std::fma(alpha_n, b[n - 1], a[n]);
    b[n] = std::fma(alpha_n, c[n - 1], b[n]);
    d[n] = std::fma(alpha_n, d[n - 1], d[n]);
}

#endif
