#pragma once

#include <cmath>
#include <cstddef>
#include <functional>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>

namespace jw {

    template<typename T>
    struct vec2 {
        T x{};
        T y{};

        constexpr bool operator==(const vec2<T>& rhs) const = default;

        constexpr vec2<T>& operator+=(const vec2<T>& rhs) {
            x += rhs.x;
            y += rhs.y;
            return *this;
        }

        constexpr vec2<T>& operator-=(const vec2<T>& rhs) {
            x -= rhs.x;
            y -= rhs.y;
            return *this;
        }

        constexpr vec2<T>& operator*=(const T& scalar) {
            x *= scalar;
            y *= scalar;
            return *this;
        }

        constexpr vec2<T>& operator/=(const T& scalar) {
            x /= scalar;
            y /= scalar;
            return *this;
        }
    };

    template<typename T>
    [[nodiscard]]
    constexpr vec2<T> operator+(const vec2<T>& v) {
        return v;
    }

    template<typename T>
    [[nodiscard]]
    constexpr vec2<T> operator-(const vec2<T>& v) {
        return {
            -v.x,
            -v.y
        };
    }

    template<typename T>
    [[nodiscard]]
    constexpr vec2<T> operator+(vec2<T> lhs, const vec2<T>& rhs) {
        lhs += rhs;
        return lhs;
    }

    template<typename T>
    [[nodiscard]]
    constexpr vec2<T> operator-(vec2<T> lhs, const vec2<T>& rhs) {
        lhs -= rhs;
        return lhs;
    }

    template<typename T>
    [[nodiscard]]
    constexpr vec2<T> operator*(vec2<T> lhs, const T& scalar) {
        lhs *= scalar;
        return lhs;
    }

    template<typename T>
    [[nodiscard]]
    constexpr vec2<T> operator*(const T& scalar, vec2<T> rhs) {
        rhs *= scalar;
        return rhs;
    }

    template<typename T>
    [[nodiscard]]
    constexpr vec2<T> operator/(vec2<T> lhs, const T& scalar) {
        lhs /= scalar;
        return lhs;
    }

    template<typename T>
    [[nodiscard]]
    constexpr T dot(const vec2<T>& lhs, const vec2<T>& rhs) {
        return lhs.x * rhs.x + lhs.y * rhs.y;
    }

    // 2D wedge product / scalar cross product.
    template<typename T>
    [[nodiscard]]
    constexpr T cross(const vec2<T>& lhs, const vec2<T>& rhs) {
        return lhs.x * rhs.y - lhs.y * rhs.x;
    }

    template<typename T>
    [[nodiscard]]
    constexpr T length_squared(const vec2<T>& v) {
        return dot(v, v);
    }

    template<typename T>
    [[nodiscard]]
    double length(const vec2<T>& v) {
        return std::sqrt(static_cast<double>(length_squared(v)));
    }

    template<typename T>
    [[nodiscard]]
    constexpr T manhattan_length(const vec2<T>& v) {
        if constexpr (std::is_signed_v<T>) {
            const T ax = v.x < T{} ? -v.x : v.x;
            const T ay = v.y < T{} ? -v.y : v.y;
            return ax + ay;
        } else {
            return v.x + v.y;
        }
    }

    template<typename T>
    [[nodiscard]]
    constexpr T manhattan_distance(const vec2<T>& lhs, const vec2<T>& rhs) {
        return manhattan_length(lhs - rhs);
    }

    template<typename T>
    [[nodiscard]]
    constexpr vec2<T> min_components(const vec2<T>& lhs, const vec2<T>& rhs) {
        return {
            lhs.x < rhs.x ? lhs.x : rhs.x,
            lhs.y < rhs.y ? lhs.y : rhs.y
        };
    }

    template<typename T>
    [[nodiscard]]
    constexpr vec2<T> max_components(const vec2<T>& lhs, const vec2<T>& rhs) {
        return {
            lhs.x > rhs.x ? lhs.x : rhs.x,
            lhs.y > rhs.y ? lhs.y : rhs.y
        };
    }

    template<typename T>
    [[nodiscard]]
    constexpr vec2<T> abs_components(const vec2<T>& v) {
        if constexpr (std::is_signed_v<T>) {
            return {
                v.x < T{} ? -v.x : v.x,
                v.y < T{} ? -v.y : v.y
            };
        } else {
            return v;
        }
    }

    template<typename T>
    [[nodiscard]]
    constexpr vec2<T> sign_components(const vec2<T>& v) {
        return {
            (T{} < v.x) - (v.x < T{}),
            (T{} < v.y) - (v.y < T{})
        };
    }

    template<typename T>
    [[nodiscard]]
    constexpr vec2<T> clockwise_90(const vec2<T>& v) {
        return {
            v.y,
            -v.x
        };
    }

    template<typename T>
    [[nodiscard]]
    constexpr vec2<T> counterclockwise_90(const vec2<T>& v) {
        return {
            -v.y,
            v.x
        };
    }

    namespace detail {

        inline void hash_combine(std::size_t& seed, std::size_t value) {
            // Same basic mixing constant traditionally used by boost::hash_combine,
            // but with no Boost dependency.
            seed ^= value + 0x9e3779b97f4a7c15ull + (seed << 6) + (seed >> 2);
        }

    }

    template<typename T>
    struct hash_vec2 {
        [[nodiscard]]
        std::size_t operator()(const vec2<T>& pt) const {
            std::size_t seed = 0;
            detail::hash_combine(seed, std::hash<T>{}(pt.x));
            detail::hash_combine(seed, std::hash<T>{}(pt.y));
            return seed;
        }
    };

    template<typename T>
    using vec2_set = std::unordered_set<vec2<T>, hash_vec2<T>>;

    template<typename T, typename U>
    using vec2_map = std::unordered_map<vec2<T>, U, hash_vec2<T>>;

}