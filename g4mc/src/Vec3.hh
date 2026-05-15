#ifndef Vec3_HH
#define Vec3_HH

#include <cmath> 

/// @brief very basic struct representing a 3-vector
struct Vec3 {
    double x, y, z; 
    
    static double norm(const Vec3& v) { return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z); }
};

#endif