// Natural Neighbor Interpolation following the CGAL User Manual:
// https://doc.cgal.org/latest/Interpolation/index.html

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;

#include <cstdint>
#include <stdint.h>
#include <ctime>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K>    Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>      Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>        Delaunay_triangulation;
typedef Delaunay_triangulation::Vertex_handle         Vertex_handle;
typedef K::FT                                         Coord_type;
typedef K::Point_2                                    Point;
typedef std::vector< std::pair<Point, Coord_type> >   Coordinate_vector;
typedef std::map<Point, Coord_type, K::Less_xy_2>     Point_value_map;
typedef std::vector<std::vector<double>>              CVdouble2d;
using Vdouble = nb::ndarray<double, nb::ndim<1>, nb::device::cpu>;
using Vdouble2d = nb::ndarray<double, nb::ndim<2>, nb::device::cpu>;

template<typename T>
std::vector<std::vector<T>> make_2d_vector(std::size_t rows, std::size_t cols)
{
    return std::vector<std::vector<T>>(rows, std::vector<T>(cols));
}

Delaunay_triangulation triangulate(Vdouble2d candidates) {

    int ntrainpoints = candidates.shape(0);
    std::vector< std::pair<Point,unsigned> > points;
    points.reserve(ntrainpoints);

    Delaunay_triangulation T;
    for(int i=0; i<ntrainpoints ; i++){
        Point p(candidates(i,0), candidates(i,1));
        points.push_back(std::make_pair(p,i));
    }

    T.insert( points.begin(),points.end() );
    return T;
}

CVdouble2d nninterpol(Delaunay_triangulation T, Vdouble z, Vdouble2d xx, Vdouble2d yy, double fill_value) {
    int dst_x_size = xx.shape(0);
    int dst_y_size = xx.shape(1);
    CVdouble2d raster = make_2d_vector<double>(dst_x_size, dst_y_size);
    Point_value_map values;

    for (Delaunay_triangulation::Vertex_handle v : T.finite_vertex_handles()) {
          values.insert(std::make_pair(v->point(), z(v->info())));
    }

    Delaunay_triangulation::Face_handle fh;
    for(int i=0; i<dst_x_size; i++) {
        for(int j=0; j<dst_y_size; j++) {
            Point p(xx(i,j), yy(i,j));
            fh = T.locate(p, fh);
            std::vector< std::pair< Point, Coord_type > > coords;
            CGAL::Triple<std::back_insert_iterator<Coordinate_vector>, K::FT, bool> natneighbor = CGAL::natural_neighbor_coordinates_2(T, p, std::back_inserter(coords),fh);

            if(!natneighbor.third)
            {
                raster[i][j] = fill_value;
            }
            else
            {
                raster[i][j] = CGAL::linear_interpolation(coords.begin(), coords.end(), natneighbor.second, CGAL::Data_access<Point_value_map>(values));
            }
        }
    }
    return raster;
}

NB_MODULE(nninterpol, m) {
    m.doc() = "Natural Neighbor Interpolation using CGAL";
    nb::class_<Delaunay_triangulation>(m, "Delaunay_triangulation");
    m.def("triangulate", &triangulate, "Triangulate the source grid using CGAL");
    m.def("nninterpol", &nninterpol,nb::call_guard<nb::gil_scoped_release>(), "Linear Interpolation using CGAL");
}
