#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>
#include "Eigen/Dense"
#include "boost/math/distributions.hpp"

using namespace std;
using namespace Eigen;
using namespace boost::math;

struct  Point
{
	string benchmark_name;
	double height;
};

struct Meas
{
	string from_station;
	string to_station;
	double elevation;
	double distance;
	string leveling_class;
};
