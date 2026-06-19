#pragma once
#include <cmath>
#include <vector>
#include "WormholeTypes.hpp"
#include "GraphicsTypes.hpp"
#include "GraphicsMath.hpp"
#include "Noise.hpp"

namespace Wormhole {
    QuadraticLUT getLUT(const Settings& settings, const RingSettings& ringSettings, int resolution) {
        QuadraticLUT LUT;
        LUT.table.reserve(resolution);
        
        double accumulatedFlow_min = 0.0;
        double accumulatedFlow_mid = 0.0;
        double accumulatedFlow_max = 0.0;
        
        double k = static_cast<double>(ringSettings.a) / static_cast<double>(ringSettings.b);
        double period = M_PI / k;

        double epsilon_sq = ringSettings.epsilon * ringSettings.epsilon;
        
        double p_term = settings.r0 * (1.0 - epsilon_sq) / 2.0;
        
        double p_min = p_term * std::sqrt(ringSettings.r_multiplier_min * ringSettings.r_multiplier_min - 1);
        double p_max = p_term * std::sqrt(ringSettings.r_multiplier_max * ringSettings.r_multiplier_max - 1);
        
        double gamma = 1 / std::sqrt(1 - ringSettings.p_phi * ringSettings.p_phi);
        
        double stepSize = period / resolution;
        
        double r0_sq = settings.r0 * settings.r0;
        double w_term = ringSettings.p_phi / gamma;

        double flow_term = stepSize / w_term;
        

        for (int i = 0; i < resolution; ++i) {
            double phi = i * stepSize;

            double cos_term = std::cos(k * phi);
            double l_term = (2.0 * cos_term) / (1.0 - epsilon_sq * cos_term * cos_term);
                
            double l_min = l_term * p_min;
            double l_max = l_term * p_max;

            double l_mid = (l_min + l_max) / 2.0;

            double r_min_sq = r0_sq + l_min * l_min;
            double r_mid_sq = r0_sq + l_mid * l_mid;
            double r_max_sq = r0_sq + l_max * l_max;
            
            double flow_min = flow_term * r_min_sq;
            double flow_mid = flow_term * r_mid_sq;
            double flow_max = flow_term * r_max_sq;
            
            accumulatedFlow_min += flow_min;
            accumulatedFlow_mid += flow_mid;
            accumulatedFlow_max += flow_max;

            QuadraticLUTField field;
            field.min = accumulatedFlow_min;
            field.mid = accumulatedFlow_mid;
            field.max = accumulatedFlow_max;

            LUT.table.push_back(field);
        }

        LUT.finalize();
        return LUT;
    }

    double getNoiseDensity(const QuadraticLUT& LUT, double period, double phiTotal, double t_theta, double t_l, double time) {
        int i = phiTotal * LUT.resolution / period;

        double t_phi = LUT.sample(t_l, i) + time;

        return Noise::fBm(t_l, t_theta, t_phi);
    }


    RayState getDerivatives(const RayState& s, double r0) {
        RayState d;
    
        double r2 = s.l * s.l + r0 * r0;
        
        double sinTheta = std::sin(s.theta);
        double cosTheta = std::cos(s.theta);
            
        if (std::abs(sinTheta) < 1e-6) sinTheta = 1e-6;
    
    
        d.l = s.p_l;
        d.theta = s.p_theta / r2;
        d.phi = s.p_phi / (r2 * sinTheta * sinTheta);
    
        d.p_l = s.l / (r2 * r2) * (s.p_theta * s.p_theta + s.p_phi * s.p_phi / (sinTheta * sinTheta));
        d.p_theta = cosTheta * s.p_phi * s.p_phi / (r2 * sinTheta * sinTheta * sinTheta);
        d.p_phi = 0.0;
            
        return d;
    }
    
    RayState rk4Step(const RayState& currentState, double stepSize, double r0) {
    
        RayState k1 = getDerivatives(currentState, r0);
    
        RayState k2 = getDerivatives(currentState + k1 * (stepSize * 0.5), r0);
    
        RayState k3 = getDerivatives(currentState + k2 * (stepSize * 0.5), r0);
    
        RayState k4 = getDerivatives(currentState + k3 * stepSize, r0);
    
        RayState finalDerivative = (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (1.0 / 6.0);
    
        return currentState + finalDerivative * stepSize;
    }
    
    
    bool rayMissesBoundarySphere(const RayState& ray, double r2, double boundaryRadiusSq) {
    
        double sinTheta = std::sin(ray.theta);
        double angularMomentumSq = (r2 * r2 * ray.p_theta * ray.p_theta) + 
                                   (r2 * r2 * sinTheta * sinTheta * ray.p_phi * ray.p_phi);
    
        if (angularMomentumSq > boundaryRadiusSq) {
            return true; 
        }
    
        if (r2 > boundaryRadiusSq && ray.p_l * ray.l >= 0.0) {
            return true; 
        }
    
        return false;
    }
    
    DustPoint sampleDustPoint(const RayState& currentState, const Settings& settings, const RingSettings& ringSettings, const QuadraticLUT& LUT, double time) {
        DustPoint dustPoint;
        double theta_mean = M_PI / 2;
        double d_theta = currentState.theta - theta_mean;
    
    
        if (std::abs(d_theta) < ringSettings.d_theta_max) {
    
            double delta_theta = d_theta / ringSettings.d_theta_max;

            double t_theta = (d_theta + ringSettings.d_theta_max) / (2 * ringSettings.d_theta_max);
    
    
            double r = std::sqrt(currentState.l * currentState.l + settings.r0 * settings.r0);
            double sinTheta = std::sin(currentState.theta);
            double gamma = 1 / std::sqrt(1 - ringSettings.p_phi * ringSettings.p_phi);
            double w_phi = ringSettings.p_phi / (gamma * r * r * sinTheta * sinTheta);
    
    
    
            double k = static_cast<double>(ringSettings.a) / static_cast<double>(ringSettings.b);
            double period = M_PI / k;
    
            double p_term = settings.r0 * (1 - ringSettings.epsilon * ringSettings.epsilon) / 2;
    
            double p_min = p_term * std::sqrt(ringSettings.r_multiplier_min * ringSettings.r_multiplier_min - 1);
            double p_max = p_term * std::sqrt(ringSettings.r_multiplier_max * ringSettings.r_multiplier_max - 1);
    
    
    
            for (int i = 0; i < ringSettings.b; ++i) {
                double phiTotal = currentState.phi - 2 * M_PI * i;
                double cos_term = std::cos(k * phiTotal);
                double l_term = (2.0 * cos_term) / (1.0 - ringSettings.epsilon * ringSettings.epsilon * cos_term * cos_term);
                
                double l_min = l_term * p_min;
                double l_max = l_term * p_max;
    
                double l_mean = (l_min + l_max) / 2;
                double l_spread = (l_max - l_min) / 2;
    
                double delta_l = (currentState.l - l_mean) / l_spread;
                double delta2 = (delta_l * delta_l) + (delta_theta * delta_theta);
                if (delta2 < 1.0) {
    
                        
    
                    double dustArea = std::abs(l_spread) * r * d_theta;
                    double v_phi = r * sinTheta * w_phi;
    
                    double fluxDensity = ringSettings.baseIntensity / (dustArea * v_phi);
                    double falloff = 1.0 - delta2;
                    double density = fluxDensity * falloff * falloff * falloff;

                    double t_l = (currentState.l - l_min) / (2 * l_spread);

                    double noiseDensity = getNoiseDensity(LUT, period, phiTotal, t_theta, t_l, time);
                }
                
            }
        }
    
        return dustPoint;
    }
    

    
    
    
    
    
    Graphics::ColorRGB castRay(RayState& ray, const Settings& settings, const RingSettings& ringSettings, const Graphics::Skybox& skybox, const QuadraticLUT& LUT, double time) {
    
        double boundaryRadiusSq = settings.r0 * settings.r0 * 15.0 * 15.0;
        double r2 = ray.l * ray.l + settings.r0 * settings.r0;
    
    
        //Ray missed wormhole and dust
        if (rayMissesBoundarySphere(ray, r2, boundaryRadiusSq)) {
            Graphics::Vec3 rayDir = Graphics::getCartesianDirection(ray, settings.r0);
            return Graphics::sampleSkyBox(rayDir, skybox);
        }
    
    
    
        Graphics::ColorRGB pixelColor = Graphics::ColorRGB(0.0, 0.0, 0.0);
        double remainingTransmittance = 1.0;
    
        double stepSize;
        bool rayStuck = true;
    
    
        //Casting ray loop
        for (int i = 0; i < settings.maxSteps; ++i) {
    
            //Breaking if pixel completely colored
            if (remainingTransmittance < 0.001) {
                remainingTransmittance = 0.0;
                rayStuck = false;
                break;
            }
    
            r2 = ray.l * ray.l + settings.r0 * settings.r0;
    
    
            //Breaking if ray is out of wormhole and dust
            if (r2 > boundaryRadiusSq) {
                if (ray.p_l * ray.l >= 0.0) {
                    rayStuck = false;
                    break;
                }
                stepSize = settings.biggerStepSize;
            }
            else {
                stepSize = settings.smallerStepSize;
            }
    
    
            DustPoint dustPoint = sampleDustPoint(ray, settings, ringSettings, LUT, time);
    
    
            //Updating pixel color if inside dust
            if (dustPoint.density > 0.0) {
                double stepOpacity = dustPoint.density * stepSize;
    
                pixelColor = pixelColor + dustPoint.color * stepOpacity * remainingTransmittance;
    
                remainingTransmittance *= 1.0 - stepOpacity;
            }
    
    
            //Updating ray
            ray = rk4Step(ray, stepSize, settings.r0);
        }
    
    
        //Adding skybox color
        if (remainingTransmittance > 0.001) {
            Graphics::ColorRGB skyBoxColor;
            if (rayStuck) {
                skyBoxColor = Graphics::ColorRGB(0.0, 0.0, 0.0);
            }
            else {
                Graphics::Vec3 rayDir = Graphics::getCartesianDirection(ray, settings.r0);
                skyBoxColor = Graphics::sampleSkyBox(rayDir, skybox);
            }
            pixelColor = pixelColor + skyBoxColor * remainingTransmittance;
        }
    
        return pixelColor;
    }

    Graphics::Frame renderFrame(const Camera& camera, const Settings& settings, const RingSettings& ringSettings, const Graphics::Skybox& skybox, const QuadraticLUT& LUT, double time) {
        Graphics::Frame frame;
        
        double width = 2 * std::tan(camera.fov / 2);
        double height = width / camera.aspectRatio;
        
        int widthResolution = camera.aspectRatio * camera.height;

        frame.width = widthResolution;
        frame.height = camera.height;

        double widthStep = width / widthResolution;
        double heightStep = height / camera.height;

        double sinTheta = std::sin(camera.theta_dir);
        double sinPhi = std::sin(camera.phi_dir);
        double cosPhi = std::cos(camera.phi_dir);

        Graphics::Vec3 forward(cosPhi * sinTheta, std::cos(camera.theta_dir), sinPhi * sinTheta);
        Graphics::Vec3 right(-sinPhi, 0.0, cosPhi);
        Graphics::Vec3 up = Graphics::cross(forward, right);

        Graphics::Vec3 start = forward - right * width / 2 - up * height / 2;

        double cosTheta_pos = std::cos(camera.theta);
        double sinTheta_pos = std::sin(camera.theta);
        double cosPhi_pos = std::cos(camera.phi);
        double sinPhi_pos = std::sin(camera.phi);

        Graphics::Vec3 e_l(cosPhi_pos * sinTheta_pos, cosTheta_pos, sinPhi_pos * sinTheta_pos);
        Graphics::Vec3 e_theta(cosPhi_pos * cosTheta_pos, -sinTheta_pos, sinPhi_pos * cosTheta_pos);
        Graphics::Vec3 e_phi(-sinPhi_pos, 0.0, cosPhi_pos);

        double p_phiTerm = (std::abs(sinTheta_pos) > 1e-6) ? (1.0 / (camera.l * sinTheta_pos)) : 0.0;

        RayState ray;

        for (int y = 0; y < camera.height; ++y) {
            for (int x = 0; x < widthResolution; ++x) {
                Graphics::Vec3 rayDir = start + right * widthStep * x + up * heightStep * y;
                rayDir = rayDir.normalized();

                ray.l = camera.l;
                ray.theta = camera.theta;
                ray.phi = camera.phi;

                ray.p_l = Graphics::dot(rayDir, e_l);
                ray.p_theta = Graphics::dot(rayDir, e_theta) / camera.l;
                ray.p_phi = Graphics::dot(rayDir, e_phi) * p_phiTerm;
                
                frame.buffer.push_back(castRay(ray, settings, ringSettings, skybox, LUT, time));
            }
        }
    }

    void updateCamera(Camera& camera) {
        
    }
}