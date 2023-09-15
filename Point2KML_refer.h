// WaypointKMLGenerator.h

#ifndef _WAYPOINT_KML_GENERATOR_H
#define _WAYPOINT_KML_GENERATOR_H

#include <vector>
#include <string>
#include "tinyxml2.h"

// 声明Waypoint结构体
typedef struct tagWaypoint
{
    double x, y, z;
    double lon, lat, alt;
    double pitch, heading;
    double baseLine;
    int nStripID;
} Waypoint;

// 声明generateKMLFromWaypoints函数
bool generateKMLFromWaypoints(const std::vector<Waypoint>& M300KML, const std::string& outputFilename);

#endif // _WAYPOINT_KML_GENERATOR_H
