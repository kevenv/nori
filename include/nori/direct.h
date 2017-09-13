#pragma once

#include <nori/integrator.h>

class Scene;

NORI_NAMESPACE_BEGIN

class DirectIntegrator : public Integrator {
public:
    DirectIntegrator(const PropertyList& props);
    // used via another integrator
    DirectIntegrator(const std::string& samplingMethod, int emitterSamples, int brdfSamples);

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray, const Intersection* _its=nullptr) const override;

    std::string toString() const override;

private:
    const std::string m_samplingMethod;
    const int m_emitterSamples;
    const int m_brdfSamples;

    float balanceHeuristic(int n1, float pdf1, int n2, float pdf2) const;
    float powerHeuristic(int n1, float pdf1, int n2, float pdf2, float power) const;
};

NORI_NAMESPACE_END