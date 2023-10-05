#ifndef GENERATE_PATH_H
#define GENERATE_PATH_H

#include <iostream>
#include <opencv2/opencv.hpp>
#include <gdal.h>
#include <ogr_geometry.h>
#include <ogr_feature.h>
#include <ogrsf_frmts.h>
#include <vector>
#include <constants.h>
#include <Cluster_and_direct_planning.h>
#include <sstream> // 用于 std::ostringstream

// Struct declarations
struct Point3D {
    double x, y, z;
public:
    bool operator==(const Point3D& other) const {
        return this->x == other.x && this->y == other.y && this->z == other.z;
    }
    static double distance(const Point3D& a, const Point3D& b) {
        return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) +
                         (a.z - b.z) * (a.z - b.z));
    }
};

struct Vector3D {
    double x, y, z;
};
struct Polygon3D {
    Point3D points[4];

    // 支持迭代的begin和end方法
    auto begin() -> Point3D* {
        return std::begin(points);
    }

    auto end() -> Point3D* {
        return std::end(points);
    }

    // const版本的begin和end方法
    auto begin() const -> const Point3D* {
        return std::begin(points);
    }

    auto end() const -> const Point3D* {
        return std::end(points);
    }

    // size函数，返回多边形的点数
    size_t size() const {
        return std::size(points);
    }
};
inline Vector3D crossProduct(Vector3D vec1, Vector3D vec2) {
    Vector3D crossVec;
    crossVec.x = vec1.y * vec2.z - vec1.z * vec2.y;
    crossVec.y = vec1.z * vec2.x - vec1.x * vec2.z;
    crossVec.z = vec1.x * vec2.y - vec1.y * vec2.x;
    return crossVec;
}
inline Vector3D makeZPositive(Vector3D input){
    if(input.z < 0){
        input.x = -input.x;
        input.y = -input.y;
        input.z = -input.z;
    }
    return input;
}

inline Vector3D averageNormalVector(Polygon3D polygon) {
    // Calculate vectors for two triangles (A, B, C) and (A, C, D)
    Vector3D vec1, vec2, vec3, vec4;
    vec1.x = polygon.points[1].x - polygon.points[0].x;
    vec1.y = polygon.points[1].y - polygon.points[0].y;
    vec1.z = polygon.points[1].z - polygon.points[0].z;

    vec2.x = polygon.points[2].x - polygon.points[0].x;
    vec2.y = polygon.points[2].y - polygon.points[0].y;
    vec2.z = polygon.points[2].z - polygon.points[0].z;

    vec3.x = polygon.points[2].x - polygon.points[3].x;
    vec3.y = polygon.points[2].y - polygon.points[3].y;
    vec3.z = polygon.points[2].z - polygon.points[3].z;

    vec4.x = polygon.points[0].x - polygon.points[3].x;
    vec4.y = polygon.points[0].y - polygon.points[3].y;
    vec4.z = polygon.points[0].z - polygon.points[3].z;

    // Calculate the cross products to get the normal vectors
    Vector3D normalVec1 = makeZPositive(crossProduct(vec1, vec2));
    Vector3D normalVec2 = makeZPositive(crossProduct(vec3, vec4));

    // Calculate the average normal vector
    Vector3D averageNormalVec;
    averageNormalVec.x = (normalVec1.x + normalVec2.x) / 2.0;
    averageNormalVec.y = (normalVec1.y + normalVec2.y) / 2.0;
    averageNormalVec.z = (normalVec1.z + normalVec2.z) / 2.0;

    // Optionally, normalize the average normal vector
    double magnitude = std::sqrt(averageNormalVec.x * averageNormalVec.x +
                                 averageNormalVec.y * averageNormalVec.y +
                                 averageNormalVec.z * averageNormalVec.z);
    averageNormalVec.x /= magnitude;
    averageNormalVec.y /= magnitude;
    averageNormalVec.z /= magnitude;


    return averageNormalVec;
}
// Function prototypes
double distance(const Point3D& a, const Point3D& b);
std::vector<Point3D> greedyTSP(const Point3D& start, const std::vector<Point3D>& points);
std::vector<std::vector<Point3D>> interpolatePoints(const std::vector<std::vector<Point3D>>& originalResult, double assumedZ, double interval);
std::vector<Point3D> interpolatePoints(const std::vector<Point3D>& originalResult, double assumedZ, double interval);
bool generateKMLsFromWaypoints(const std::vector<std::vector<Waypoint>>& waypointsGroups, const std::string& baseFilename);
#endif // GENERATE_PATH_H
