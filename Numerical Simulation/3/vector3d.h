#include <math.h>
using namespace std;
struct v3{
    double x, y, z;

    // default + parameterized constructor
    v3(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z){}

    // assignment operator modifies object, therefore non-const
    v3& operator=(const v3& a){
        x = a.x;
        y = a.y;
        z = a.z;
        return *this;
    }

    // addop. doesn't modify object. therefore const.
    v3 operator+(const v3& a) const{
        return v3(a.x + x, a.y + y, a.z + z);
    }

    // equality comparison. doesn't modify object. therefore const.
    bool operator==(const v3& a) const{
        return (x == a.x && y == a.y && z == a.z);
    }
    
};

double dot(v3 a, v3 b){
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

v3 cross(v3 a, v3 b){
	return v3(a.y*b.z - a.z*b.y, -a.x*b.z + a.z*b.x, a.x*b.y - a.y*b.x);
}

double magnitude(v3 a){
	return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

v3 neg(v3 a){
	return v3(-a.x, -a.y, -a.z);
}

v3 vadd(v3 a, v3 b){
	return v3(a.x + b.x, a.y + b.y, a.z + b.z);
}
