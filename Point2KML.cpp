﻿// Point2KML.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "afxwin.h"
#include "afxdialogex.h"
#include "atlconv.h"

#include "tinyxml2.h"



using namespace std;
using namespace tinyxml2;

#define PI 3.141592653589793
#define MATH_PI         3.14159265358979323846264338327950288   // pi
#define MATH_PI_2       1.57079632679489661923132169163975144   // pi/2
#define MATH_PI_4       0.785398163397448309615660845819875721  // pi/4
#define MATH_R2D        57.295779513082320876798154814105e0
#define MATH_D2R        0.017453292519943295769236907684886e0

#ifndef _WAY_POINT
#define _WAY_POINT
typedef struct tagWaypoint
{
	double x, y, z;
	double lon, lat, alt;
	double pitch, heading;
	double baseLine;
	int nStripID;
}Waypoint;
#endif




//字符串分割函数
CStringArray* Split(CString str, char ch)
{
	// 去掉两边的ch<br />
	str.TrimLeft(ch);
	str.TrimRight(ch);
	int nStart = 0;
	int nLastStart = 0;
	CStringArray* pStrArray = new CStringArray();
	while (-1 != (nStart = str.Find(ch, nStart))) {

		if (nLastStart != nStart) {// 不是连续ch,<br />
			pStrArray->Add(str.Mid(nLastStart, nStart - nLastStart));
		}
		nStart++;
		nLastStart = nStart;
	}
	pStrArray->Add(str.Mid(nLastStart));	// 最后一个子串<br />
	return pStrArray;
}

template<class TYPE>
tinyxml2::XMLElement* InsertChildElement(tinyxml2::XMLDocument& doc,
	tinyxml2::XMLElement* pItem, const char* lpstrTag, TYPE val, bool bVal = true)
{
	tinyxml2::XMLElement* pEle = doc.NewElement(lpstrTag);
	if (bVal) pEle->SetText(val);
	pItem->InsertEndChild(pEle);
	return pEle;
}


int main(int argc, char* argv[])
{
	CStringArray* strMatchNameArray = Split(argv[1], '.');
	string szXmlFileName = strMatchNameArray->GetAt(0) + ".kml";

	ifstream inFile5(argv[1]);
	vector<string> v_str;

	if (inFile5)
	{
		string strLine5;
		while (getline(inFile5, strLine5)) // line中不包括每行的换行符  
		{
			v_str.push_back(strLine5);
		}
	}


	vector<Waypoint> M300KML;


	for (int i = 0; i < v_str.size(); i = i + 1)
	{
		CStringArray* strMatchFileArray = Split(v_str[i].c_str(), ',');

		Waypoint test_Aero;

		test_Aero.lon = atof(strMatchFileArray->GetAt(0));
		test_Aero.lat = atof(strMatchFileArray->GetAt(1));
		test_Aero.alt = atof(strMatchFileArray->GetAt(2));

		test_Aero.heading = atof(strMatchFileArray->GetAt(3));
		test_Aero.pitch = atof(strMatchFileArray->GetAt(4));

		M300KML.push_back(test_Aero);
	}

	int nWaypoint = M300KML.size();



	tinyxml2::XMLDocument doc;
	tinyxml2::XMLDeclaration* declaration = doc.NewDeclaration();
	doc.InsertFirstChild(declaration);

	tinyxml2::XMLElement* root = doc.NewElement("kml");
	if (root == NULL) {
		return false;
	}
	root->SetAttribute("xmlns", "http://www.opengis.net/kml/2.2");
	doc.InsertEndChild(root);

	tinyxml2::XMLElement* pDoc = doc.NewElement("Document");
	pDoc->SetAttribute("xmlns", "");
	root->InsertEndChild(pDoc);

	InsertChildElement(doc, pDoc, "name", strMatchNameArray->GetAt(0));
	InsertChildElement(doc, pDoc, "open", 1);

	tinyxml2::XMLElement* pExt = doc.NewElement("ExtendedData");
	pExt->SetAttribute("xmlns:mis", "www.dji.com");
	pDoc->InsertEndChild(pExt);
	InsertChildElement(doc, pExt, "mis:type", "Waypoint");
	InsertChildElement(doc, pExt, "mis:stationType", "0");

	// waylineGreenPoly
	tinyxml2::XMLElement* pLineStyle = doc.NewElement("Style");
	pLineStyle->SetAttribute("id", "waylineGreenPoly");
	pDoc->InsertEndChild(pLineStyle);
	tinyxml2::XMLElement* pLine = doc.NewElement("LineStyle");
	pLineStyle->InsertEndChild(pLine);
	InsertChildElement(doc, pLine, "color", "FF0AEE8B");
	InsertChildElement(doc, pLine, "width", 6);

	// waypointStyle
	tinyxml2::XMLElement* pLinePoint = doc.NewElement("Style");
	pLinePoint->SetAttribute("id", "waypointStyle");
	pDoc->InsertEndChild(pLinePoint);
	tinyxml2::XMLElement* pIconStyle = doc.NewElement("IconStyle");
	pLinePoint->InsertEndChild(pIconStyle);
	tinyxml2::XMLElement* pIcon = doc.NewElement("Icon");
	pIconStyle->InsertEndChild(pIcon);
	InsertChildElement(doc, pIcon, "href", "https://cdnen.dji-flighthub.com/static/app/images/point.png");

	// Folder
	tinyxml2::XMLElement* pFolder = doc.NewElement("Folder");
	pDoc->InsertEndChild(pFolder);
	InsertChildElement(doc, pFolder, "name", "Waypoints");
	InsertChildElement(doc, pFolder, "description", "Waypoints in the Mission.");



	for (int i = 0; i < nWaypoint; i++) {
		tinyxml2::XMLElement* pPlacemark = doc.NewElement("Placemark");
		pFolder->InsertEndChild(pPlacemark);

		char strValue[128] = ""; sprintf(strValue, "Waypoint%d", i + 1);
		InsertChildElement(doc, pPlacemark, "name", strValue);
		InsertChildElement(doc, pPlacemark, "visibility", 1);
		InsertChildElement(doc, pPlacemark, "description", "Waypoint");
		InsertChildElement(doc, pPlacemark, "styleUrl", "#waypointStyle");

		tinyxml2::XMLElement* pExtendedData = doc.NewElement("ExtendedData");
		pExtendedData->SetAttribute("xmlns:mis", "www.dji.com");
		pPlacemark->InsertEndChild(pExtendedData);
		InsertChildElement(doc, pExtendedData, "mis:useWaylineAltitude", "false");
		//InsertChildElement(doc, pExtendedData, "mis:heading", int(M300KML[i].heading * MATH_R2D + 0.5));
		InsertChildElement(doc, pExtendedData, "mis:heading", int(M300KML[i].heading));
		InsertChildElement(doc, pExtendedData, "mis:turnMode", "Auto");
		//InsertChildElement(doc, pExtendedData, "mis:gimbalPitch", int(M300KML[i].pitch * MATH_R2D + 0.5));
		InsertChildElement(doc, pExtendedData, "mis:gimbalPitch", int(M300KML[i].pitch));
		InsertChildElement(doc, pExtendedData, "mis:useWaylineSpeed", "true");
		InsertChildElement(doc, pExtendedData, "mis:speed", "3.0");
		InsertChildElement(doc, pExtendedData, "mis:useWaylineHeadingMode", "false");
		InsertChildElement(doc, pExtendedData, "mis:useWaylinePointType", "false");
		InsertChildElement(doc, pExtendedData, "mis:pointType", "LineStop");
		InsertChildElement(doc, pExtendedData, "mis:headingMode", "UsePointSetting");
		InsertChildElement(doc, pExtendedData, "mis:cornerRadius", "0.2");

		tinyxml2::XMLElement* pMisAction = doc.NewElement("mis:actions");


		pMisAction->SetAttribute("param", "0");
		pMisAction->SetAttribute("accuracy", "0");
		pMisAction->SetAttribute("cameraIndex", "0");
		pMisAction->SetAttribute("payloadType", "0");
		pMisAction->SetAttribute("payloadIndex", "0");
		pMisAction->SetText("ShootPhoto");

		pExtendedData->InsertEndChild(pMisAction);


		tinyxml2::XMLElement* pPoint = doc.NewElement("Point");
		pPlacemark->InsertEndChild(pPoint);
		InsertChildElement(doc, pPoint, "altitudeMode", "relativeToGround");
		//sprintf(strValue, "%.8lf,%.8lf,%.3lf", pWaypoint[i].lon, pWaypoint[i].lat, pWaypoint[i].alt - heiOfHome);
		sprintf(strValue, "%.8lf,%.8lf,%.3lf", M300KML[i].lon, M300KML[i].lat, M300KML[i].alt);
		InsertChildElement(doc, pPoint, "coordinates", strValue);
	}

	// Placemark
	tinyxml2::XMLElement* pPlacemark = doc.NewElement("Placemark");
	pDoc->InsertEndChild(pPlacemark);
	InsertChildElement(doc, pPlacemark, "name", "Wayline");
	InsertChildElement(doc, pPlacemark, "description", "Wayline");
	InsertChildElement(doc, pPlacemark, "visibility", 1);

	tinyxml2::XMLElement* pExtendedData = doc.NewElement("ExtendedData");
	pExtendedData->SetAttribute("xmlns:mis", "www.dji.com");
	pPlacemark->InsertEndChild(pExtendedData);
	InsertChildElement(doc, pExtendedData, "mis:altitude", 50);
	InsertChildElement(doc, pExtendedData, "mis:autoFlightSpeed", "4.0");
	InsertChildElement(doc, pExtendedData, "mis:actionOnFinish", "GoHome");
	InsertChildElement(doc, pExtendedData, "mis:headingMode", "Auto");
	InsertChildElement(doc, pExtendedData, "mis:gimbalPitchMode", "UsePointSetting");
	InsertChildElement(doc, pExtendedData, "mis:powerSaveMode", "false");
	InsertChildElement(doc, pExtendedData, "mis:waypointType", "LineStop");



	tinyxml2::XMLElement* pdroneInfo = doc.NewElement("mis:droneInfo");
	pExtendedData->InsertEndChild(pdroneInfo);
	InsertChildElement(doc, pdroneInfo, "mis:droneType", "PM430");
	InsertChildElement(doc, pdroneInfo, "mis:advanceSettings", "true");

	tinyxml2::XMLElement* pdroneCameras = doc.NewElement("mis:droneCameras");
	pdroneInfo->InsertEndChild(pdroneCameras);

	tinyxml2::XMLElement* pdroneCamera = doc.NewElement("mis:camera");
	pdroneCameras->InsertEndChild(pdroneCamera);
	InsertChildElement(doc, pdroneCamera, "mis:cameraIndex", "0");
	InsertChildElement(doc, pdroneCamera, "mis:cameraName", "Zenmuse P1");
	InsertChildElement(doc, pdroneCamera, "mis:cameraType", "31");
	InsertChildElement(doc, pdroneCamera, "mis:payloadCameraType", "2");


	tinyxml2::XMLElement* pdroneHeight = doc.NewElement("mis:droneHeight");
	pdroneInfo->InsertEndChild(pdroneHeight);
	InsertChildElement(doc, pdroneHeight, "mis:useAbsolute", "true");
	InsertChildElement(doc, pdroneHeight, "mis:hasTakeoffHeight", "false");
	InsertChildElement(doc, pdroneHeight, "mis:takeoffHeight", "0.0");


	InsertChildElement(doc, pPlacemark, "styleUrl", "#waylineGreenPoly");
	tinyxml2::XMLElement* pLineString = doc.NewElement("LineString");
	pPlacemark->InsertEndChild(pLineString);
	InsertChildElement(doc, pLineString, "tessellate", "1");
	InsertChildElement(doc, pLineString, "altitudeMode", "relativeToGround");

	std::string coordinatesStr;
	for (int i = 0; i < nWaypoint; i++) {
		char strValue[128] = "";
		//sprintf(strValue, "%.8lf,%.8lf,%.3lf ", pWaypoint[i].lon, pWaypoint[i].lat, pWaypoint[i].alt - heiOfHome);
		sprintf(strValue, "%.8lf,%.8lf,%.3lf ", M300KML[i].lon, M300KML[i].lat, M300KML[i].alt);
		coordinatesStr += strValue;
	}

	InsertChildElement(doc, pLineString, "coordinates", coordinatesStr.c_str());
	return (tinyxml2::XML_SUCCESS == doc.SaveFile(szXmlFileName.c_str())) ? true : false;

}

