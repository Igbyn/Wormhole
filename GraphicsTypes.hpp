#pragma once
#include <cmath>

namespace Graphics {
    struct Vec3 {
        double x, y, z;
        Vec3(double x=0, double y=0, double z=0) : x(x), y(y), z(z) {}
    
        Vec3 operator+(const Vec3& v) const {
            return Vec3(x + v.x, y + v.y, z + v.z);
        }
    
        Vec3 operator-(const Vec3& v) const {
            return Vec3(x - v.x, y - v.y, z - v.z);
        }
    
        Vec3 operator*(double scalar) const {
            return Vec3(x * scalar, y * scalar, z * scalar);
        }
    
        Vec3 operator/(double scalar) const {
            double inv = 1.0 / scalar;
            return Vec3(x * inv, y * inv, z * inv);
        }
    
        Vec3 operator*(const Vec3& v) const {
            return Vec3(x * v.x, y * v.y, z * v.z);
        }
    
        Vec3 normalized() const {
            double magSq = x * x + y * y + z * z;
            if (magSq > 0.0) {
                double invMag = 1.0 / std::sqrt(magSq);
                return Vec3(x * invMag, y * invMag, z * invMag);
            }
            return Vec3(0.0, 0.0, 0.0);
        }
    };
    
    inline Vec3 operator*(double factor, const Vec3& v) { return v * factor; }
    
    inline double dot(const Vec3& a, const Vec3& b) {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }
    
    inline Vec3 cross(const Vec3& a, const Vec3& b) {
        return Vec3(
            a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x
        );
    }
    
    
    struct ColorRGB {
        double r = 0.0, g = 0.0, b = 0.0;
        
        ColorRGB(double r=0, double g=0, double b=0) : r(r), g(g), b(b) {}
        
        ColorRGB operator*(const ColorRGB& other) const {
            return {r * other.r, g * other.g, b * other.b};
        }
    
        ColorRGB operator+(const ColorRGB& other) const {
            return {r + other.r, g + other.g, b + other.b};
        }
        
        ColorRGB operator*(double factor) const {
            return {r * factor, g * factor, b * factor};
        }
    };
    
    inline ColorRGB operator*(double factor, const ColorRGB& c) { return c * factor; }
    
    
    
    struct Skybox {
        float* data;
        int width;
        int height;
    };

    struct Frame {
        int width = 0;
        int height = 0;

        std::vector<ColorRGB> buffer;

        void resize(int w, int h) {
            width = w;
            height = h;
            buffer.resize(width * height);
        }

        inline ColorRGB& operator()(int x, int y) {
            return buffer[y * width + x];
        }

        inline const ColorRGB& operator()(int x, int y) const {
            return buffer[y * width + x];
        }
    };
}