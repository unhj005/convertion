#include "WebCoordSystemTransform.h"
#include <iostream>
#include <cmath>
using namespace std;

const double pi = 3.14159265358979324;

//
// Krasovsky 1940
//
// a = 6378245.0, 1/f = 298.3
// b = a * (1 - f)
// ee = (a^2 - b^2) / a^2;
const double a = 6378245.0;
const double ee = 0.00669342162296594323;



void WEBCSTRANSFORM::wgs2gcj( double wgLat, double wgLon, double& mgLat, double& mgLon )
{
	if (outOfChina(wgLat, wgLon))
	{
		mgLat = wgLat;
		mgLon = wgLon;
		return;
	}
	double dLat = transformLat(wgLon - 105.0, wgLat - 35.0);
	double dLon = transformLon(wgLon - 105.0, wgLat - 35.0);
	double radLat = wgLat / 180.0 * pi;
	double magic = sin(radLat);
	magic = 1 - ee * magic * magic;
	double sqrtMagic = sqrt(magic);
	dLat = (dLat * 180.0) / ((a * (1 - ee)) / (magic * sqrtMagic) * pi);
	dLon = (dLon * 180.0) / (a / sqrtMagic * cos(radLat) * pi);
	mgLat = wgLat + dLat;
	mgLon = wgLon + dLon;
}

bool WEBCSTRANSFORM::outOfChina( double lat, double lon )
{
	if (lon < 72.004 || lon > 137.8347)
		return true;
	if (lat < 0.8293 || lat > 55.8271)
		return true;
	return false;
}

double WEBCSTRANSFORM::transformLat( double x, double y )
{
	double ret = -100.0 + 2.0 * x + 3.0 * y + 0.2 * y * y + 0.1 * x * y + 0.2 * sqrt(abs(x));
	ret += (20.0 * sin(6.0 * x * pi) + 20.0 * sin(2.0 * x * pi)) * 2.0 / 3.0;
	ret += (20.0 * sin(y * pi) + 40.0 * sin(y / 3.0 * pi)) * 2.0 / 3.0;
	ret += (160.0 * sin(y / 12.0 * pi) + 320 * sin(y * pi / 30.0)) * 2.0 / 3.0;
	return ret;
}

double WEBCSTRANSFORM::transformLon( double x, double y )
{
	double ret = 300.0 + x + 2.0 * y + 0.1 * x * x + 0.1 * x * y + 0.1 * sqrt(abs(x));
	ret += (20.0 * sin(6.0 * x * pi) + 20.0 * sin(2.0 * x * pi)) * 2.0 / 3.0;
	ret += (20.0 * sin(x * pi) + 40.0 * sin(x / 3.0 * pi)) * 2.0 / 3.0;
	ret += (150.0 * sin(x / 12.0 * pi) + 300.0 * sin(x / 30.0 * pi)) * 2.0 / 3.0;
	return ret;
}

void WEBCSTRANSFORM::gcj2wgs( double mglat, double mglon, double& wglat, double& wglon )
{
	//¶þ·Ö·¨
	double glat, glon;
	wgs2gcj(mglat, mglon, glat, glon);
	double dlat = glat - mglat;
	double dlon = glon - mglon;
	wglat = mglat - dlat;
	wglon = mglon - dlon;
}

void WEBCSTRANSFORM::gcj2bd09( double mglat, double mglon, double& bdlat, double& bdlon )
{
	double x = mglon, y = mglat;  
	double z = sqrt(x * x + y * y) + 0.00002 * sin(y * pi);  
	double theta = atan2(y, x) + 0.000003 * cos(x * pi);  
	bdlon = z * cos(theta) + 0.0065;  
	bdlat = z * sin(theta) + 0.006;  
}

void WEBCSTRANSFORM::bd092gcj( double bdlat, double bdlon, double& mglat, double& mglon )
{
	double x = bdlon - 0.0065, y = bdlat - 0.006;  
	double z = sqrt(x * x + y * y) - 0.00002 * sin(y * pi);  
	double theta = atan2(y, x) - 0.000003 * cos(x * pi);  
	mglon = z * cos(theta);  
	mglat = z * sin(theta);  
}

void WEBCSTRANSFORM::bd092wgs( double bdlat, double bdlon, double& wglat, double& wglon )
{
	double mglat, mglon;
	bd092gcj(bdlat, bdlon, mglat, mglon);
	gcj2wgs(mglat, mglon, wglat, wglon);
}

void WEBCSTRANSFORM::wgs2bd09( double wglat, double wglon, double& bdlat, double& bdlon )
{
	double mglat, mglon;
	wgs2gcj(wglat, wglon, mglat, mglon);
	gcj2bd09(mglat, mglon, bdlat, bdlon);
}
