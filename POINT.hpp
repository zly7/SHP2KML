#pragma once

struct POINT3D
{
	POINT3D() { X = 0; Y = 0; Z = 0; }
	POINT3D(double x, double y, double z) :X(x), Y(y), Z(z) {}
	~POINT3D() {}
	double X, Y, Z;
};

struct POINT2D
{
	POINT2D() { X = 0; Y = 0; }
	POINT2D(double x, double y) :X(x), Y(y) {}
	~POINT2D() {}
	double X, Y;
};

struct POINTLINK //关联像片点与地面点的结构体
{
	POINTLINK() { PhotoPoint = POINT2D(0, 0); VectorPoint = POINT3D(0, 0, 0); }
	POINTLINK(POINT2D& a, POINT3D& b) :PhotoPoint(a), VectorPoint(b) {}
	~POINTLINK() {}
	POINT2D PhotoPoint;
	POINT3D VectorPoint;
};