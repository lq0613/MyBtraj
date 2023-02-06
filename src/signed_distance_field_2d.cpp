/**
 * Copyright (C) 2022, RAM-LAB, Hong Kong University of Science and Technology
 * This file is part of GPIR (https://github.com/jchengai/gpir).
 * If you find this repo helpful, please cite the respective publication as
 * listed on the above website.
 */

#include "signed_distance_field_2d.h"

#include <omp.h>

#include <future>
#include <iostream>
// #include <opencv2/core.hpp>
// #include <opencv2/core/eigen.hpp>
// #include <opencv2/highgui.hpp>
// #include <opencv2/imgproc.hpp>
//#include "common/matplotlib-cpp/matplotlibcpp.h"


using namespace planning;
using namespace std;
using namespace Eigen;

namespace planning {
  
SignedDistanceField2D::SignedDistanceField2D(){};
SignedDistanceField2D::SignedDistanceField2D(std::array<double, 2> origin,
                                             std::array<int, 2> dim,
                                             const double map_resolution)
    : map_resolution_(map_resolution) {
  occupancy_map_.set_origin(origin);
  occupancy_map_.set_cell_number(dim);
  occupancy_map_.set_resolution(
      std::array<double, 2>{map_resolution, map_resolution});
  esdf_.ResizeFrom(occupancy_map_);
}

SignedDistanceField2D::SignedDistanceField2D(OccupancyMap&& occupancy_map)
    : occupancy_map_(std::move(occupancy_map)) {
  esdf_.ResizeFrom(occupancy_map_);
  map_resolution_ = occupancy_map_.resolution()[0];
}

void SignedDistanceField2D::EuclideanDistanceTransform(
    std::array<int, 2> dim,
    std::function<bool(const int x, const int y)> is_occupied,
    DistanceMap* output_map) {
  int inf = dim[0] + dim[1] + 10;

  std::vector<std::vector<int>> g(dim[0], std::vector<int>(dim[1], 0));
  omp_set_num_threads(4);
  {
#pragma omp parallel for
    // column scan
    for (int x = 0; x < dim[0]; ++x) {
      g[x][0] = is_occupied(x, 0) ? 0 : inf;

      for (int y = 1; y < dim[1]; ++y) {
        g[x][y] = is_occupied(x, y) ? 0 : 1 + g[x][y - 1];
      }

      for (int y = dim[1] - 2; y >= 0; --y) {
        if (g[x][y + 1] < g[x][y]) g[x][y] = 1 + g[x][y + 1];
      }
    }
  }

  // row scan
  omp_set_num_threads(4);
  {
#pragma omp parallel for
    for (int y = 0; y < dim[1]; ++y) {
      int q = 0, w;
      std::vector<int> s(dim[0], 0);
      std::vector<int> t(dim[0], 0);

      auto f = [&g, &y](int x, int i) -> double {
        return (x - i) * (x - i) + g[i][y] * g[i][y];
      };

      for (int u = 1; u < dim[0]; ++u) {
        while (q >= 0 && f(t[q], s[q]) > f(t[q], u)) {
          --q;
        }

        if (q < 0) {
          q = 0;
          s[0] = u;
        } else {
          w = 1 + std::floor((u * u - s[q] * s[q] + g[u][y] * g[u][y] -
                              g[s[q]][y] * g[s[q]][y]) /
                             (2 * (u - s[q])));
          if (w < dim[0]) {
            ++q;
            s[q] = u;
            t[q] = w;
          }
        }
      }

      for (int u = dim[0] - 1; u >= 0; --u) {
        output_map->SetValue(u, y, map_resolution_ * std::sqrt(f(u, s[q])));
        if (u == t[q]) --q;
      }
    }
  }
}

void SignedDistanceField2D::VerticalEuclideanDistanceTransform(
    std::array<int, 2> dim,
    std::function<bool(const int x, const int y)> is_occupied,
    DistanceMap* output_map) {
  int inf = 1e9;

  std::vector<std::vector<int>> g(dim[0], std::vector<int>(dim[1], 0));
  omp_set_num_threads(4);
  {
#pragma omp parallel for
    // column scan
    for (int x = 0; x < dim[0]; ++x) {
      g[x][0] = is_occupied(x, 0) ? 0 : inf;

      for (int y = 1; y < dim[1]; ++y) {
        g[x][y] = is_occupied(x, y) ? 0 : 1 + g[x][y - 1];
      }

      for (int y = dim[1] - 2; y >= 0; --y) {
        if (g[x][y + 1] < g[x][y]) g[x][y] = 1 + g[x][y + 1];
      }

      for (int y = 0; y < dim[1]; ++y) {
        output_map->SetValue(x, y, map_resolution_ * g[x][y]);
      }
    }
  }
}

void SignedDistanceField2D::UpdateSDF() {
  DistanceMap distance_map, inv_distance_map;
  distance_map.ResizeFrom(occupancy_map_);
  inv_distance_map.ResizeFrom(occupancy_map_);
  OccupancyMap& map = occupancy_map_;

  auto dim = occupancy_map_.cell_num();
  EuclideanDistanceTransform(
      dim,
      [&map](const int x, const int y) -> bool { return map.IsOccupied(x, y); },
      &distance_map);
  EuclideanDistanceTransform(
      dim,
      [&map](const int x, const int y) -> bool {
        return !map.IsOccupied(x, y);
      },
      &inv_distance_map);

  const auto& dis_map_data = distance_map.data();
  const auto& inv_dis_map_data = inv_distance_map.data();
  auto& esdf_data = *esdf_.mutable_data();

  omp_set_num_threads(4);
  {
#pragma omp parallel for
    for (int x = 0; x < dim[0]; ++x) {
      for (int y = 0; y < dim[1]; ++y) {
        int address = occupancy_map_.Index2Address(x, y);
        esdf_data[address] = dis_map_data[address];
        if (inv_dis_map_data[address] > 0) {
          esdf_data[address] += (-inv_dis_map_data[address] + map_resolution_);
        }
      }
    }
  }
}

void SignedDistanceField2D::UpdateVerticalSDF() {
  DistanceMap distance_map, inv_distance_map;
  distance_map.ResizeFrom(occupancy_map_);
  inv_distance_map.ResizeFrom(occupancy_map_);
  OccupancyMap& map = occupancy_map_;

  auto dim = occupancy_map_.cell_num();
  VerticalEuclideanDistanceTransform(
      dim,
      [&map](const int x, const int y) -> bool { return map.IsOccupied(x, y); },
      &distance_map);
  VerticalEuclideanDistanceTransform(
      dim,
      [&map](const int x, const int y) -> bool {
        return !map.IsOccupied(x, y);
      },
      &inv_distance_map);

  const auto& dis_map_data = distance_map.data();
  const auto& inv_dis_map_data = inv_distance_map.data();
  auto& esdf_data = *esdf_.mutable_data();

  omp_set_num_threads(4);
  {
#pragma omp parallel for
    for (int x = 0; x < dim[0]; ++x) {
      for (int y = 0; y < dim[1]; ++y) {
        int address = occupancy_map_.Index2Address(x, y);
        esdf_data[address] = dis_map_data[address];
        if (inv_dis_map_data[address] > 0) {
          esdf_data[address] += (-inv_dis_map_data[address] + map_resolution_);
        }
      }
    }
  }
}

void SignedDistanceField2D::TestBasic(){
  int dim = 10;
  std::array<int, 2> size{dim, dim};

  GridMap2D<uint8_t> grid_map;
  grid_map.set_cell_number(size);
  grid_map.set_origin(std::array<double, 2>{0.0, 0.0});
  grid_map.set_resolution(std::array<double, 2>{1, 1});

  for (int i = 0; i < dim; ++i) {
    grid_map.SetValue(Eigen::Vector2i(i, i), 1);
  }
   //cv::imshow("grid map", grid_map.BinaryImage());
  SignedDistanceField2D sdf(std::move(grid_map));
  sdf.UpdateSDF();

  auto esdf = sdf.esdf();

  std::cout << esdf.Matrix() << std::endl;
   //cv::imshow("sdf", sdf.esdf().ImageSec());
   //cv::waitKey(0);
  // EXPECT_EQ(esdf.GetValue(Vector2i(0, 2)), 0.1 * std::sqrt(2));
  // EXPECT_EQ(esdf.GetValue(Vector2i(0, 1)), 0.1 * 1);
}

void SignedDistanceField2D::TestSignedDistance(){
  int dim = 500;
  std::array<int, 2> size{800, 200};

  GridMap2D<uint8_t> grid_map;
  grid_map.set_cell_number(size);
  grid_map.set_origin(std::array<double, 2>{0.0, 0.0});
  grid_map.set_resolution(std::array<double, 2>{0.1, 0.1});

  grid_map.FillCircle(Eigen::Vector2d(15, 15), 10);
  grid_map.FillConvexPoly(
      vector_Eigen<Vector2d>{Vector2d(30, 15), Vector2d(40, 15),
                       Vector2d(40, 20), Vector2d(30, 20)});
  grid_map.FillPoly(vector_Eigen<Vector2d>{Vector2d(0, 20), Vector2d(20, 20),
                                     Vector2d(25, 25), Vector2d(15, 22),
                                     Vector2d(0, 30)});
 // cv::imshow("grid map", grid_map.BinaryImage());

  SignedDistanceField2D sdf(std::move(grid_map));
  auto t0 = chrono::high_resolution_clock::now();

  sdf.UpdateSDF();
  auto t1 = chrono::high_resolution_clock::now();

  double total_ms =
      chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000.0;

  cout << "time for 500x500 sdf: " << total_ms << " ms" << endl;
  // cv::imshow("sdf", sdf.esdf().ImageSec());
  // cv::waitKey(0);
}
void SignedDistanceField2D::TestSignedPolyLines(){
  std::array<int, 2> size{800, 200};
  GridMap2D<uint8_t> grid_map;
  grid_map.set_cell_number(size);
  grid_map.set_origin(std::array<double, 2>{0.0, 0.0});
  grid_map.set_resolution(std::array<double, 2>{0.1, 0.1});

  vector_Eigen<Vector2d> points;
  for (int i = 0; i < 140; ++i) {
    points.emplace_back(Eigen::Vector2d(0.5 * i, 10 + 3 * std::sin(0.5 * i)));
  }

  grid_map.PolyLine(points);
  // cv::imshow("f(x) = 3*sin(x)", grid_map.BinaryImage());
  // cv::waitKey(0);
}
// void SignedDistanceField2D::TestSignedGradient(){
//   int dim = 50;
//   std::array<int, 2> size{dim, dim};

//   GridMap2D<uint8_t> grid_map;
//   grid_map.set_cell_number(size);
//   grid_map.set_origin(std::array<double, 2>{0.0, 0.0});
//   grid_map.set_resolution(std::array<double, 2>{1.0, 1.0});
//   grid_map.FillCircle(Eigen::Vector2d(25, 25), 10);
//   grid_map.FillPoly(vector_Eigen<Vector2d>{Vector2d(10, 10), Vector2d(20, 10),
//                                      Vector2d(20, 20), Vector2d(10, 20)});

//   SignedDistanceField2D sdf(std::move(grid_map));
//   sdf.UpdateSDF();

//   auto esdf = sdf.esdf();

//   std::vector<double> x, y, u, v;
//   for (int i = 0; i < dim; ++i) {
//     for (int j = 0; j < dim; ++j) {
//       x.emplace_back(i + 0.5);
//       y.emplace_back(j + 0.5);
//       Eigen::Vector2d grad;
//       esdf.GetValueBilinear(Eigen::Vector2d(i + 0.5, j + 0.5), &grad);
//       u.emplace_back(grad.x());
//       v.emplace_back(grad.y());
//     }
//   }

//   auto image = esdf.ImageSec();
//   cv::Mat scaled_image;
//   cv::resize(image, scaled_image, {0, 0}, 10.0, 10.0);
//   cv::imshow("sdf", scaled_image);
//   cv::waitKey(0);

//   // namespace plt = matplotlibcpp;
//   // plt::quiver(x, y, u, v);
//   // plt::axis("equal");
//   // plt::show();
// }
}