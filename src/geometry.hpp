#pragma once

#include "vec2.hpp"
#include <Eigen/Dense>
#include <vector>
#include <tuple>

/*------------------------------------------------------------------------------------------------*/

namespace jw {

    using point = vec2<double>;
    using points = std::vector<point>;
    using matrix = Eigen::Matrix<double, 3, 3>;

    struct line_segment {
        point u;
        point v;
    };

    auto rotation_matrix(double cos_theta, double sin_theta) -> matrix;
    auto translation_matrix(double x, double y) -> matrix;
    auto scale_matrix(double x_scale, double y_scale) -> matrix;
    auto scale_matrix(double scale)-> matrix;
    auto identity_matrix() -> matrix;
    auto translation_matrix(point pt) -> matrix;
    auto distance(point u, point v) -> double;
    auto apply_matrix(const matrix& matrix, const point& pt) -> point;
    auto apply_matrix(const matrix& matrix, const points& pt) -> points;
    auto apply_matrix(const matrix& matrix, const line_segment& seg) -> line_segment;
    auto map_line_seg_on_to_line_seg_matrix(
        const line_segment& from, const line_segment& to) -> matrix;
    auto map_pt_to_pt_and_rotate(point from, point to, double theta) -> matrix;
    auto inset_line_seg_matrix(const line_segment& seg, double scale_len, double distance) -> matrix;

    std::vector<point> equilateral_triangle_vertices(double base);
    std::vector<point> equilateral_trapezoid_vertices(double base, double side);
}