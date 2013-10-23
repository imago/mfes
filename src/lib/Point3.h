#include <string>
#include <iostream>
#include <ostream>

using namespace std;

class Point3 {
public:
	Point3(float _x, float  _y, float _z):
		x(_x), y(_y), z(_z)
	{}

	float x;
	float y;
	float z;

	friend ostream& operator<< (std::ostream& stream, Point3& p){
		stream << "(" << p.x << ", " << p.y << ", " << p.z << ")";
		return stream;
	}
};
