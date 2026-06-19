#pragma once
#include <cmath>
#include "GraphicsTypes.hpp"





namespace Wormhole
{
    struct QuadraticLUTField {
        double min;
        double mid;
        double max;
    };

    struct QuadraticLUT {
        std::vector<QuadraticLUTField> table;
        int resolution = 0;

        inline void finalize() {
            resolution = static_cast<int>(table.size());
        }

        inline double sample(double t, int i) const {
            int cycle = i / resolution;
            int j = i % resolution;

            const QuadraticLUTField& field = table[j];
            const QuadraticLUTField& lastField = table.back();

            double wMin = (1.0 - 2.0 * t) * (1.0 - t);
            double wMid = 4.0 * t * (1.0 - t);
            double wMax = t * (2.0 * t - 1.0);

            return (lastField.min * cycle + field.min) * wMin
                 + (lastField.mid * cycle + field.mid) * wMid
                 + (lastField.max * cycle + field.max) * wMax;
        }
    };

    struct DustPoint {
        Graphics::ColorRGB color;
        double density;
    };
    
    
    struct UserInputState {
        bool w = false;
        bool a = false;
        bool s = false;
        bool d = false;
        bool shift = false;
        bool space = false;
        double mouseDx = 0.0;
        double mouseDy = 0.0;
    };
    
    
    struct Camera {
        double l, theta, phi;
        double theta_dir, phi_dir;

        double movementSpeed = 0.5;
        double rotationMultiplier = 0.1;
        
        double fov = M_PI / 3.0;
        int height = 360;
        double aspectRatio = 16.0 / 9.0;

        void updateAspectRatio(int width, int height) {
            aspectRatio = static_cast<double>(width) / static_cast<double>(height);
        }

        void updateCamera(const UserInputState& inputState, double dt) {
            double forwardMultiplier = dt * movementSpeed * (static_cast<double>(inputState.w) - static_cast<double>(inputState.s));
            double rightMultiplier = dt * movementSpeed * (static_cast<double>(inputState.d) - static_cast<double>(inputState.a));
            double upMultiplier = dt * movementSpeed * (static_cast<double>(inputState.space) - static_cast<double>(inputState.shift));

            double dphi = rotationMultiplier * inputState.mouseDx;
            double dtheta = rotationMultiplier * inputState.mouseDy;

            phi_dir += dphi;
            theta_dir -= dtheta;

            Graphics::Vec3 forward(std::cos(phi_dir) * std::sin(theta_dir), std::cos(theta_dir), std::sin(phi_dir) * std::sin(theta_dir));
            Graphics::Vec3 right(-std::sin(phi_dir), 0.0, std::cos(phi_dir));
            Graphics::Vec3 up = Graphics::cross(forward, right);

            forward = forward * forwardMultiplier;
            right = right * rightMultiplier;
            up = up * upMultiplier;

            Graphics::Vec3 dir = forward + right + up;

            Graphics::Vec3 e_l(std::cos(phi) * std::sin(theta), std::cos(theta), std::sin(phi) * std::sin(theta));
            Graphics::Vec3 e_theta(std::cos(phi) * std::cos(theta), -std::sin(theta), std::sin(phi) * std::cos(theta));
            Graphics::Vec3 e_phi(-std::sin(phi), 0.0, std::cos(phi));

            l += Graphics::dot(dir, e_l);
            theta += Graphics::dot(dir, e_theta) / l;
            phi += Graphics::dot(dir, e_phi) / (l * std::sin(theta));
        }
    };

    
    
    
    struct RingSettings {
        double r_multiplier_min = 3;
        double r_multiplier_max = 5;
        double epsilon = std::sqrt(1.0/2.0);
        int a = 3;
        int b = 5;
        double p_phi = 0.5;
        double d_theta_max = M_PI / 18.0;
        double baseIntensity = 1.0;
    };
    
    
    
    struct Settings {
        double r0 = 5.0;
        double smallerStepSize = 0.05;
        double biggerStepSize = 0.1;
        int maxSteps = 200;
    };
    
    
    
    struct RayState {
        double l, theta, phi;
        double p_l, p_theta, p_phi;
        
        RayState operator+(const RayState& a) const {
            return {l+a.l, theta+a.theta, phi+a.phi, p_l+a.p_l, p_theta+a.p_theta, p_phi+a.p_phi};
        }
        
        RayState operator-(const RayState& a) const {
            return {l-a.l, theta-a.theta, phi-a.phi, p_l-a.p_l, p_theta-a.p_theta, p_phi-a.p_phi};
        }
        
        RayState operator*(double factor) const {
            return {l*factor, theta*factor, phi*factor, p_l*factor, p_theta*factor, p_phi*factor};
        }
    };
    
    inline RayState operator*(double factor, const RayState& r) { return r * factor; }
}
