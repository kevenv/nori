
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class AmbientOcclusionIntegrator : public Integrator {
public:
    AmbientOcclusionIntegrator(const PropertyList &props) :
        m_sampleCount(props.getInteger("sampleCount", 1)),
        m_albedo(props.getFloat("albedo", 1.0f)),
        m_samplingMethod(props.getString("samplingMethod", "cosine-weighted"))
    {

    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray, const Intersection* _its) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        Normal3f n = its.shFrame.n;
        float maxt = scene->getBoundingBox().getExtents().norm();

        // monte carlo!
        float Lr = 0.0f;
        if (m_samplingMethod == "uniform") {

            for (int i = 0; i < m_sampleCount; ++i) {
                Vector3f d = Warp::squareToUniformHemisphere(sampler->next2D());
                d = its.toWorld(d); // transform to world space so it aligns with the its
                d.normalize();
                Ray3f shadowRay(its.p, d, Epsilon, maxt);
                float V = (float)!scene->rayIntersect(shadowRay);
                float cosTheta = std::max(0.0f, d.dot(n));
                Lr += V * cosTheta;
            }
            Lr *= 2.0f * m_albedo / m_sampleCount;

        }
        else {

            for (int i = 0; i < m_sampleCount; ++i) {
                Vector3f d = Warp::squareToCosineHemisphere(sampler->next2D());
                d = its.toWorld(d); // transform to world space so it aligns with the its
                d.normalize();
                Ray3f shadowRay(its.p, d, Epsilon, maxt);
                float V = (float)!scene->rayIntersect(shadowRay);
                Lr += V;
            }
            Lr *= m_albedo / m_sampleCount;

        }

        return Color3f(Lr);
    }

    std::string toString() const {
        return tfm::format(
            "AmbientOcclusionIntegrator[\n"
            " samplingMethod = %s\n"
            " sampleCount = %d\n"
            " albedo = %d\n"
            "]",
            m_samplingMethod,
            m_sampleCount,
            m_albedo
        );
    }

private:
    const std::string m_samplingMethod;
    const int m_sampleCount;
    const float m_albedo;
};

NORI_REGISTER_CLASS(AmbientOcclusionIntegrator, "ao");
NORI_NAMESPACE_END