#include "geometry.hpp"
#include <ranges>
#include <cmath>

namespace r = std::ranges;
namespace rv = std::ranges::views;

/*------------------------------------------------------------------------------------------------*/

auto jw::rotation_matrix(double cos_theta, double sin_theta) -> matrix {

	matrix rotation;
	rotation <<
		cos_theta, -sin_theta, 0,
		sin_theta, cos_theta, 0,
		0, 0, 1;

	return rotation;
}

auto jw::translation_matrix(double x, double y) -> matrix {
	matrix translation;
	translation <<
		1, 0, x,
		0, 1, y,
		0, 0, 1;
	return translation;
}

auto jw::scale_matrix(double x_scale, double y_scale) -> matrix {
	
	matrix scale;
	scale <<
		x_scale, 0, 0,
		0, y_scale, 0,
		0, 0, 1;

	return scale;
}

auto jw::scale_matrix(double scale) -> matrix {
	
	matrix mat;
	mat <<
		scale, 0, 0,
		0, scale, 0,
		0, 0, 1;

	return mat;
}

auto jw::identity_matrix() -> matrix {

	return scale_matrix( 1.0 );

}

auto jw::translation_matrix(point pt) -> matrix {
	return translation_matrix( pt.x, pt.y );
}

auto jw::distance(point u, point v) -> double {
	auto [x1, y1] = u;
	auto [x2, y2] = v;
	auto x_diff = x2 - x1;
	auto y_diff = y2 - y1;
	return std::sqrt(x_diff * x_diff + y_diff * y_diff);
}

auto jw::apply_matrix(const matrix& matrix, const point& pt) -> point {
    Eigen::Matrix<double, 3, 1> v;
    v << pt.x, pt.y, 1.0;
    v = matrix * v;
    return { v[0], v[1] };
}

auto jw::apply_matrix(const matrix& matrix, const points& poly) -> points {
    return poly | rv::transform(
            [&matrix](auto pt) {
                return apply_matrix(matrix, pt);
            }
        ) | r::to<std::vector>();
}

auto jw::apply_matrix(const matrix& matrix, const line_segment& seg) -> line_segment {
    return {
        apply_matrix(matrix, seg.u),
        apply_matrix(matrix, seg.v)
    };
}

auto jw::map_pt_to_pt_and_rotate(point from, point to, double theta) -> matrix {
    const double cos_theta = std::cos(theta);
    const double sin_theta = std::sin(theta);

    return translation_matrix(to)
        * rotation_matrix(cos_theta, sin_theta)
        * translation_matrix(-from.x, -from.y);
}

auto jw::inset_line_seg_matrix(const line_segment& seg, double len_marg, double distance) -> matrix {

    if (len_marg < 0.0) {
        throw std::runtime_error(
            "inset_line_seg_matrix: len_marg must be non-negative."
        );
    }

    const point delta = seg.v - seg.u;
    const double len = std::sqrt(delta.x * delta.x + delta.y * delta.y);

    if (len <= 0.0) {
        throw std::runtime_error(
            "inset_line_seg_matrix: segment must have non-zero length."
        );
    }

    const double inset_len = len - 2.0 * len_marg;

    if (inset_len <= 0.0) {
        throw std::runtime_error(
            "inset_line_seg_matrix: len_marg is too large for segment."
        );
    }

    const double scale = inset_len / len;

    // Unit direction from u to v.
    const point dir{
        delta.x / len,
        delta.y / len
    };

    // Left-hand unit normal when walking from u to v.
    const point normal{
        -dir.y,
        dir.x
    };

    /*
        Work in the local basis:

            p = seg.u + x * dir + y * normal

        Then transform as:

            x' = len_marg + scale * x
            y' = distance + y

        So endpoints become:

            seg.u -> seg.u + len_marg * dir + distance * normal
            seg.v -> seg.v - len_marg * dir + distance * normal
    */

    const double a00 = scale * dir.x * dir.x + normal.x * normal.x;
    const double a01 = scale * dir.x * dir.y + normal.x * normal.y;
    const double a10 = scale * dir.y * dir.x + normal.y * normal.x;
    const double a11 = scale * dir.y * dir.y + normal.y * normal.y;

    const point target_u =
        seg.u +
        len_marg * dir +
        distance * normal;

    const point mapped_u{
        a00 * seg.u.x + a01 * seg.u.y,
        a10 * seg.u.x + a11 * seg.u.y
    };

    const point translation = target_u - mapped_u;

    matrix result;
    result <<
        a00, a01, translation.x,
        a10, a11, translation.y,
        0.0, 0.0, 1.0;

    return result;

}

auto jw::map_line_seg_on_to_line_seg_matrix( 
        const line_segment& from, const line_segment& to ) -> matrix {

    const point from_dir = from.v - from.u;
    const point to_dir = to.v - to.u;

    const double from_len = distance(from.u, from.v);
    const double to_len = distance(to.u, to.v);

    if (from_len == 0.0) {
        throw std::runtime_error(
            "map_line_seg_on_to_line_seg_matrix: from segment has zero length."
        );
    }

    if (to_len == 0.0) {
        throw std::runtime_error(
            "map_line_seg_on_to_line_seg_matrix: to segment has zero length."
        );
    }

    const double from_x = from_dir.x / from_len;
    const double from_y = from_dir.y / from_len;

    const double to_x = to_dir.x / to_len;
    const double to_y = to_dir.y / to_len;

    const double cos_theta =  from_x * to_x + from_y * to_y;
    const double sin_theta =  from_x * to_y - from_y * to_x;

    return translation_matrix(to.u)
        * rotation_matrix(cos_theta, sin_theta)
        * translation_matrix(-from.u.x, -from.u.y);
}

std::vector<jw::point> jw::equilateral_triangle_vertices(double base) {

    const auto half_root_3 = std::sqrt(3.0) / 2.0;
    return {
        {0.0, 0.0},
        {0.5 * base , half_root_3 * base},
        {base, 0.0}
    };
}

std::vector<jw::point> jw::equilateral_trapezoid_vertices(double base, double side) {
    // side is the length of the two non-parallel sides. The left side rises
    // at +60 degrees from the x-axis; the right side falls at -60 degrees.
    //
    // Vertices are returned counter-clockwise, starting at the lower-left corner.

    const double half_side = 0.5 * side;
    const double height = std::sqrt(3.0) * half_side;

    return {
        { 0.0, 0.0 },
        { half_side, height },
        { base - half_side, height },
        { base, 0.0 }
    };
}