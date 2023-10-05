// kmeans_clustering.h
#ifndef KMEANS_CLUSTERING_H
#define KMEANS_CLUSTERING_H
#include "generate_path.h"
#include <vector>
#include <opencv2/opencv.hpp>

// 声明 kMeansClustering 函数
std::vector<std::vector<Point3D>> kMeansClustering(const std::vector<Point3D>& points, int k);
std::vector<std::vector<Point3D>> generateZigzagTrajectories(const std::vector<std::vector<Point3D>>& clusters, double d, double h);
#endif // KMEANS_CLUSTERING_H
