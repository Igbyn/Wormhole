#include <cmath>
#include <cstdint>

namespace Noise {

    double hash3D(int32_t x, int32_t y, int32_t z) {
        uint32_t h = static_cast<uint32_t>(x) * 1664525u + static_cast<uint32_t>(y) * 1013904223u + static_cast<uint32_t>(z) * 33241243u;
        h = (h ^ (h >> 16)) * 0x45d9f3bu;
        h = (h ^ (h >> 16)) * 0x45d9f3bu;
        h = h ^ (h >> 16);
        return static_cast<double>(h) / 4294967295.0;
    }
    
    double fade(double t) {
        return t * t * t * (t * (t * 6.0 - 15.0) + 10.0);
    }
    
    double lerp(double t, double a, double b) {
        return a + t * (b - a);
    }

    double sample(double x, double y, double z) {
        int32_t X = static_cast<int32_t>(std::floor(x));
        int32_t Y = static_cast<int32_t>(std::floor(y));
        int32_t Z = static_cast<int32_t>(std::floor(z));
        
        double fx = x - std::floor(x);
        double fy = y - std::floor(y);
        double fz = z - std::floor(z);
        
        double u = fade(fx);
        double v = fade(fy);
        double w = fade(fz);
        
        double c000 = hash3D(X,     Y,     Z);
        double c100 = hash3D(X + 1, Y,     Z);
        double c010 = hash3D(X,     Y + 1, Z);
        double c110 = hash3D(X + 1, Y + 1, Z);

        double c001 = hash3D(X,     Y,     Z + 1);
        double c101 = hash3D(X + 1, Y,     Z + 1);
        double c011 = hash3D(X,     Y + 1, Z + 1);
        double c111 = hash3D(X + 1, Y + 1, Z + 1);
        
        double x1 = lerp(u, c000, c100);
        double x2 = lerp(u, c010, c110);
        double y1 = lerp(v, x1, x2);
    
        double x3 = lerp(u, c001, c101);
        double x4 = lerp(u, c011, c111);
        double y2 = lerp(v, x3, x4);
    
        return lerp(w, y1, y2);
    }
    
    double fBm(double x, double y, double z, int octaves = 8, double lacunarity = 2.0, double gain = 0.5) {
        double total = 0.0;
        double amplitude = 1.0;
        double maxValue = 0.0;
    
        for (int i = 0; i < octaves; ++i) {
            total += amplitude * sample(x, y, z);
            maxValue += amplitude;

            x *= lacunarity;
            y *= lacunarity;
            z *= lacunarity;
            amplitude *= gain;
        }
        
        return total / maxValue;
    }

}