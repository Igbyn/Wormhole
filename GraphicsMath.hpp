#pragma once
#include <cmath>
#include "GraphicsTypes.hpp"
#include "WormholeTypes.hpp"

namespace Graphics {
    Vec3 getCartesianPosition(double l, double theta, double phi, double r0) {
        double r = std::sqrt(l * l + r0 * r0);
    
        return Vec3(
            r * std::sin(theta) * std::cos(phi),
            r * std::sin(theta) * std::sin(phi),
            r * std::cos(theta)
        );
    }
    
    Vec3 getCartesianDirection(const Wormhole::RayState& ray, double r0) {
        double r = std::sqrt(ray.l * ray.l + r0 * r0);
        double sinTheta = std::sin(ray.theta);
        double cosTheta = std::cos(ray.theta);
        double sinPhi = std::sin(ray.phi);
        double cosPhi = std::cos(ray.phi);
    
        double v_r = ray.p_l;
        double v_theta = r * ray.p_theta;
        double v_phi = r * sinTheta * ray.p_phi;
    
        Vec3 dir(
            v_r * sinTheta * cosPhi + v_theta * cosTheta * cosPhi - v_phi * sinPhi,
            v_r * sinTheta * sinPhi + v_theta * cosTheta * sinPhi + v_phi * cosPhi,
            v_r * cosTheta - v_theta * sinTheta
        );
    
        return dir.normalized();
    }

    ColorRGB sampleSkyBox(const Vec3& direction, const Skybox& skybox) {
        Vec3 normDir = direction.normalized();
    
        double sky_phi = std::atan2(normDir.y, normDir.x);
        double sky_theta = std::acos(normDir.z); 
    
        double sky_u = (sky_phi + M_PI) / (2.0 * M_PI);
        double sky_v = sky_theta / M_PI;
    
        sky_u += (sky_u < 0.0) ? 1.0 : (sky_u >= 1.0 ? -1.0 : 0.0);
        sky_v = (sky_v < 0.0) ? 0.0 : (sky_v > 1.0 ? 1.0 : sky_v);
    
        int sky_x = static_cast<int>(sky_u * (skybox.width - 1));
        int sky_y = static_cast<int>(sky_v * (skybox.height - 1));
    
        int pixelIndex = (sky_y * skybox.width + sky_x) * 3;
    
        double r = skybox.data[pixelIndex];
        double g = skybox.data[pixelIndex + 1];
        double b = skybox.data[pixelIndex + 2];
    
        return ColorRGB(r, g, b);
    }
}