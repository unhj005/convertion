#ifndef GCJ02_WGS84_TRANSFORM
#define GCJ02_WGS84_TRANSFORM

namespace WEBCSTRANSFORM
{
	//
	// World Geodetic System ==> Mars Geodetic System
	void wgs2gcj(double wgLat, double wgLon, double& mgLat, double& mgLon);
	// World Geodetic System <== Mars Geodetic System
	void gcj2wgs(double mglat, double mglon, double& wglat, double& wglon);

	// Mars Geodetic System ==> Baidu Geodetic System
	void gcj2bd09(double mglat, double mglon, double& bdlat, double& bdlon);
	// Mars Geodetic System <== Baidu Geodetic System
	void bd092gcj(double bdlat, double bdlon, double& mglat, double& mglon);

	// Baidu Geodetic System ==> World Geodetic System
	void bd092wgs(double bdlat, double bdlon, double& wglat, double& wglon);
	// Baidu Geodetic System <== World Geodetic System
	void wgs2bd09(double wglat, double wglon, double& bdlat, double& bdlon);

	/////////////////////////////////////////////////////////////////////
	bool outOfChina(double lat, double lon);
	double transformLat(double x, double y);
	double transformLon(double x, double y);
}



#endif