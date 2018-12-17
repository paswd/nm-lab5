#ifndef __THOMAS_ALGORITHM_HPP__
#define __THOMAS_ALGORITHM_HPP__

#include <cstddef>

#include <exception>
#include <stdexcept>
#include <type_traits>

#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

using namespace boost::numeric;

template <typename T,
    typename = std::enable_if<std::is_floating_point<T>::value>>
ublas::vector<T> thomas_algorithm(const ublas::vector<T> &,
    const ublas::vector<T> &,
    const ublas::vector<T> &,
    const ublas::vector<T> &);

template <typename T, typename>
ublas::vector<T> thomas_algorithm(const ublas::vector<T> &a,
    const ublas::vector<T> &b,
    const ublas::vector<T> &c,
    const ublas::vector<T> &d)
{
    if (a.size() != d.size() || b.size() != d.size() || c.size() != d.size() ||
        d.size() < 3U)
    {
        throw std::domain_error("Can't be solved");
    }

    const std::size_t size = d.size();
    ublas::vector<T> p(size), q(size), x(size);

    p(0) = -c(0) / b(0);
    q(0) = d(0) / b(0);

    for (std::size_t i = 1; i < size; ++i)
    {
        p(i) = -c(i) / (b(i) + a(i) * p(i - 1));
        q(i) = (d(i) - a(i) * q(i - 1)) / (b(i) + a(i) * p(i - 1));
    }

    x(size - 1) = q(size - 1);

    for (std::ptrdiff_t i = size - 2U; i >= 0; --i)
    {
        x(i) = p(i) * x(i + 1) + q(i);
    }

    return x;
}

#endif
