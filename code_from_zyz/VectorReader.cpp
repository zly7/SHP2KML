#include "VectorReader.h"

//VectorReader::~VectorReader()
//{
//	if (poDS != nullptr)
//		Close();
//}

bool VectorReader::Open(string filepath)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");  //支持中文路径
	CPLSetConfigOption("SHAPE_ENCODING", "");

	poDS = (GDALDataset*)GDALOpenEx(filepath.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL);
	if (!poDS)
	{
		cerr << "Can't Open Vector File!" << endl;
		return false;
	}
	if (poDS->GetLayerCount() < 1)
	{
		cerr << "Layer Error" << endl;
		return false;
	}

	return true;
}

bool VectorReader::Close()
{
	if (poDS != nullptr)
	{
		GDALClose(poDS);
		return true;
	}
	else
		return false;

}

bool VectorReader::Search(double MinX, double MinY, double MaxX, double MaxY, vector<POINT3D>& QueryResult)
{
	OGRLayer* poLayer = poDS->GetLayer(0);
	poLayer->ResetReading();

	OGRFeature* poFeature;

	int i = 0;
	poLayer->SetSpatialFilterRect(MinX, MinY, MaxX, MaxY);
	while ((poFeature = poLayer->GetNextFeature()) != NULL)
	{
		OGRGeometry* pGeo = poFeature->GetGeometryRef();
		OGRwkbGeometryType pGeoType = pGeo->getGeometryType();

		//cout << poFeature->GetFieldAsString("CENTER") << endl;
		string center = poFeature->GetFieldAsString("CENTER");
		int s1 = center.find_first_of(',');
		int s2 = center.find_last_of(',');
		double x = stod(center.substr(0, s1 - 1));
		double y = stod(center.substr(s1 + 1, s2 - 1));
		double z = stod(center.substr(s2 + 1, center.size()));
		QueryResult.push_back(POINT3D(x, y, z));
	}

	return true;
}
bool VectorReader::SearchWithArea(double UpPercent, double DownPercent, vector<POINT2D>& QueryResult)
{

	int im_height = 512;
	int im_width = 640;
	double u = 12e-6;
	double x, y;
	OGRLayer* poLayer = poDS->GetLayer(0);
	poLayer->ResetReading();
	OGRFeature* poFeature;

	vector<int> Areas;
	while ((poFeature = poLayer->GetNextFeature()) != NULL)
	{
		OGRGeometry* pGeo = poFeature->GetGeometryRef();
		OGRPolygon* pPoly = pGeo->toPolygon();
		//cout << pPoly->get_Area() << endl;
		Areas.push_back(pPoly->get_Area());
	}
	std::sort(Areas.begin(), Areas.end());
	int UpN = ceil(UpPercent * Areas.size()-1);
	int DownN = floor(DownPercent * Areas.size()-1);
	int UpArea = Areas[UpN]; //面积上限
	int DownArea = Areas[DownN]; //面积下限
	poLayer->ResetReading();
	OGRPoint poPoint; //临时保存中心信息
	while ((poFeature = poLayer->GetNextFeature()) != NULL)
	{
		OGRGeometry* pGeo = poFeature->GetGeometryRef();
		OGRPolygon* pPoly = pGeo->toPolygon();
		if (pPoly->get_Area() >= DownArea && pPoly->get_Area() <= UpArea)
		{
			pPoly->Centroid(&poPoint); //计算中心位置
			x = (poPoint.getX() - im_width / 2) * u; 
			y = (-poPoint.getY() + im_height / 2) * u;//跑出来的矢量坐标系上下颠倒，在此处纠正。。。
			QueryResult.push_back(POINT2D(x, y));
		}
	}
	return true;
}

void ComputePlaneSVD(Eigen::MatrixXd& m_cloud, Eigen::Vector4d& ABCD)
{
	// 1、计算质心
	Eigen::RowVector3d centroid = m_cloud.colwise().mean();
	// 2、去质心
	Eigen::MatrixXd demean = m_cloud;
	demean.rowwise() -= centroid;
	// 3、SVD分解求解协方差矩阵的特征值特征向量
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(demean, Eigen::ComputeThinU | Eigen::ComputeThinV);
	Eigen::Matrix3d V = svd.matrixV();
	Eigen::MatrixXd U = svd.matrixU();
	//Eigen::Matrix3d S = U.inverse() * demean * V.transpose().inverse();
	Eigen::Matrix3d S = U.transpose() * demean * V.transpose().inverse();
	// 5、平面的法向量a,b,c
	Eigen::RowVector3d normal;
	normal << V(0, 2), V(1, 2), V(2, 2);
	// 6、原点到平面的距离d
	double d = -normal * centroid.transpose();
	// 7、获取拟合平面的参数a,b,c,d和质心x,y,z。
	ABCD << normal.transpose(), d;
}

double ComputeZ(double x, double y, Eigen::Vector4d& ABCD)
{
	double z = (ABCD(0) * x + ABCD(1) * y + ABCD(3)) / -ABCD(2);
	return z;
}

//int main()
//{
//	VectorReader* test = new VectorReader();
//	test->Open(".\\OriginShp\\DJI_20230227125935_0005_T.shp");
//	vector<POINT2D> QueryResult;
//	test->SearchWithArea(0.95, 0.2, QueryResult);
//	//LARGE_INTEGER t1, t2, tc;
//	//QueryPerformanceFrequency(&tc);
//	//QueryPerformanceCounter(&t1);
//
//	//test->Search(0, 0, 0, 0, QueryResult);
//	//Eigen::MatrixXd Result_Mat(QueryResult.size(), 3);
//	//Eigen::Vector4d PlanePara; //ax+by+cz+d = 0 	
//	//int i = 0;
//	//for (auto iter : QueryResult)
//	//{
//	//	Result_Mat(i, 0) = get<0>(iter);
//	//	Result_Mat(i, 1) = get<1>(iter);
//	//	Result_Mat(i, 2) = get<2>(iter);
//	//	i += 1;
//	//}
//	//QueryPerformanceCounter(&t2);
//	//double time = (double)(t2.QuadPart - t1.QuadPart) / (double)tc.QuadPart;
//	//cout << "time = " << time << endl;  //输出时间（单位：ｓ）
//	////cout << Result_Mat;
//	//ComputePlaneSVD(Result_Mat, PlanePara);
//	return 0;
//}