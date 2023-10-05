#pragma once

#include <gdal.h>
#include <gdal_priv.h>	
#include <string>
#include <vector>
#include <tuple>
#include <ogrsf_frmts.h>
#include <Eigen/Dense>
#include <geos_c.h>
#include <iostream>

#include "POINT.hpp"


using namespace std;



class VectorReader
{
public:
	VectorReader() :poDS(nullptr) {}
	~VectorReader() {};
	bool Open(string filepath);
	bool Close();
	bool Search(double MinX, double MinY, double MaxX, double MaxY, vector<POINT3D>&);
	bool SearchWithArea(double UpPercent, double DownPercent, vector<POINT2D>& QueryResult); // 用于读取原始影像光伏提取结果的shp，在上下分位数范围内筛选矢量并输出中心点坐标
private:
	GDALDataset* poDS;

};

void ComputePlaneSVD(Eigen::MatrixXd& m_cloud, Eigen::Vector4d& ABCD);

double ComputeZ(double x, double y, Eigen::Vector4d& ABCD);