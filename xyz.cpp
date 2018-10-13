#include <iostream>
#include <math.h>
#include <cmath>
using namespace std;

const double PI = 3.1415926535897932384; ///< PI
const double D2R = 0.017453292519943; ///< Coefficient of convert Degree to Radian
const double R2D = 57.29577951308232; ///< Coefficient of convert Radian to Degree
const double a = 6378245.0;
const double ee = 0.00669342162296594323;


//把WGS84的空间直角坐标系XYZ，转化为大地坐标系BLH
static int xyz2BLH(const double& x, const double& y, const double& z,double& B, double& L, double& H)
{
	const double a = 6378137;
	const double b = 6356752.3142;
	const double e2 = 0.0066943799013;
	const double ee2 = 0.00673949674227;

	double sqrtxy = sqrt(x*x + y * y);
	double theta = atan(a*z / (b*sqrtxy));
	double tanB0 = (z + ee2 * b*pow(sin(theta), 3)) / (sqrtxy - e2 * a*pow(cos(theta), 3));
	double B0 = atan(tanB0);
	double DB = 1;
	while (DB > 1e-10)
	{
		double B1 = z / sqrtxy * (1 + a * e2*sin(B0) / (z*sqrt(1 - e2 * sin(B0)*sin(B0))));
		B1 = atan(B1);
		DB = fabs(B1 - B0);
		B0 = B1;
	}
	B = B0;

	double N = a / sqrt(1 - e2 * pow(sin(B), 2));
	if (x > 1e-10)
	{
		L = atan(y / x);
	}
	else if (x < -1e-10)
	{
		L = atan(y / x) + PI;
	}
	else
	{
		if (y > 0)
			L = PI * 0.5;
		else
			L = PI * 1.5;
	}
	H = sqrtxy * cos(B) + z * sin(B) - N * (1 - e2 * pow(sin(B), 2));

	B = B * R2D;
	L = L * R2D;
	if (L > 180.0)
		L = L - 360.0;
	return 0;
}

//把大地坐标系BLH,转化为WGS84的空间直角坐标系XYZ
static int BLH2xyz(double B, double L, double H, double& x, double& y, double& z)
{
	B = B * D2R; L = L * D2R;
	const double a = 6378137;
	const double e2 = 0.0066943799013;
	double N = a / sqrt(1 - e2 * pow(sin(B), 2));
	x = (N + H)*cos(B)*cos(L);
	y = (N + H)*cos(B)*sin(L);
	z = (N*(1 - e2) + H) * sin(B);
	return 0;
}

//是否越出中国国界
bool outOfChina( double lat, double lon )
{
	if (lon < 72.004 || lon > 137.8347)
		return true;
	if (lat < 0.8293 || lat > 55.8271)
		return true;
	return false;
}

//WGS84转火星坐标系的纬度转换
double transformLat( double x, double y )
{
	double ret = -100.0 + 2.0 * x + 3.0 * y + 0.2 * y * y + 0.1 * x * y + 0.2 * sqrt(abs(x));
	ret += (20.0 * sin(6.0 * x * PI) + 20.0 * sin(2.0 * x * PI)) * 2.0 / 3.0;
	ret += (20.0 * sin(y * PI) + 40.0 * sin(y / 3.0 * PI)) * 2.0 / 3.0;
	ret += (160.0 * sin(y / 12.0 * PI) + 320 * sin(y * PI / 30.0)) * 2.0 / 3.0;
	return ret;
}

//WGS84转火星坐标系的经度转换
double transformLon( double x, double y )
{
	double ret = 300.0 + x + 2.0 * y + 0.1 * x * x + 0.1 * x * y + 0.1 * sqrt(abs(x));
	ret += (20.0 * sin(6.0 * x * PI) + 20.0 * sin(2.0 * x * PI)) * 2.0 / 3.0;
	ret += (20.0 * sin(x * PI) + 40.0 * sin(x / 3.0 * PI)) * 2.0 / 3.0;
	ret += (150.0 * sin(x / 12.0 * PI) + 300.0 * sin(x / 30.0 * PI)) * 2.0 / 3.0;
	return ret;
}

//把WGS84的BLH转火星坐标系的BLH
static int wgs2gcj( double wgLat, double wgLon, double wH, double& mgLat, double& mgLon, double& mH )
{
	if (outOfChina(wgLat, wgLon))
	{
		mgLat = wgLat;
		mgLon = wgLon;
		mH = wH;
		return 0;
	}
	double dLat = transformLat(wgLon - 105.0, wgLat - 35.0);
	double dLon = transformLon(wgLon - 105.0, wgLat - 35.0);
	double radLat = wgLat / 180.0 * PI;
	double magic = sin(radLat);
	magic = 1 - ee * magic * magic;
	double sqrtMagic = sqrt(magic);
	dLat = (dLat * 180.0) / ((a * (1 - ee)) / (magic * sqrtMagic) * PI);
	dLon = (dLon * 180.0) / (a / sqrtMagic * cos(radLat) * PI);
	mgLat = wgLat + dLat;
	mgLon = wgLon + dLon;
	mH = wH;
    return 0;
}

//把火星坐标系的BLH转百度坐标系的BLH
static int gcj2bd09( double mglat, double mglon, double mH, double& bdlat, double& bdlon,double& bdH )
{
	double x = mglon, y = mglat;
	double z = sqrt(x * x + y * y) + 0.00002 * sin(y * PI);
	double theta = atan2(y, x) + 0.000003 * cos(x * PI);
	bdlon = z * cos(theta) + 0.0065;
	bdlat = z * sin(theta) + 0.006;
	bdH = mH;
	return 0;

}

//把WGS84坐标系的BLH转百度坐标系的BLH
static int wgs2bd09( double wglat, double wglon, double wH, double& bdlat, double& bdlon, double &bdH)
{
	double mglat, mglon,mH;

	wgs2gcj(wglat, wglon, wH,mglat, mglon,mH);
	gcj2bd09(mglat, mglon, mH,bdlat, bdlon,bdH);
    return 0;
}



int main()
{

	double wx = -953076.900000;
	double wy = -6542517.500000;
	double wz = 2453130.200000;
	double wB, wL, wH,bdB,bdL,bdH;
	xyz2BLH(wx, wy, wz, wB, wL, wH);

		cout <<"wx= "<<wx << endl;
	cout <<"wy= "<<wy << endl;
	cout <<"wz= "<<wz << endl<<endl;
	cout <<"wB= "<<wB << endl;
	cout <<"wL= "<<wL << endl;
	cout <<"wH= "<<wH << endl<<endl;

    wgs2bd09(wB,wL,wH,bdB,bdL,bdH);

    cout <<"translation successful！"<<endl;
    cout <<"bdB= "<<bdB << endl;
	cout <<"bdL= "<<bdL << endl;
	cout <<"bdH= "<<bdH << endl<<endl;


	double bdx,bdy,bdz;
    BLH2xyz(bdB,bdL,bdH,bdx,bdy,bdz);

    cout<<"bdx= "<<bdx<<endl;
    cout<<"bdy= "<<bdy<<endl;
    cout<<"bdz= "<<bdz<<endl;


} 