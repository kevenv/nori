
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/shape.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/direct.h>

NORI_NAMESPACE_BEGIN

class VolumePT : public Integrator {
public:
    VolumePT(const PropertyList &props) :
        m_tracerType(props.getString("tracerType", "explicit")),
        m_termination(props.getString("termination", "russian-roulette")),
        m_terminationProb(props.getFloat("terminationProb", 0.2f)),
        m_terminationBounds(props.getInteger("terminationBounds", 15)),
        m_distanceSampling(props.getString("distanceSampling", "transmittance")),
        m_sigmaA(props.getFloat("sigmaA", 0.1f)),
        m_sigmaS(props.getFloat("sigmaS", 0.5f)),
        m_sigmaT(m_sigmaA + m_sigmaS)
    {

    }

    float Tr(const nori::Point3f& x, const nori::Point3f& y) const {
        return std::exp( -m_sigmaT * (x - y).norm() );
    }

    float distSample(const std::string& method, float tmax, Sampler* sampler) const {
        if (method == "uniform") {
            return sampler->next1D() * tmax;
        }
        else { // method == "transmittance"
            return -std::log(1 - sampler->next1D()) / m_sigmaT;
        }
    }

    float distPdf(float t, const std::string& method, float tmax) const {
        if (method == "uniform") {
            return 1.0f / tmax;
        }
        else { // method == "transmittance"
            return m_sigmaT * std::exp(-m_sigmaT*t);
        }
    }

    float distPdf_failure(float s, const std::string& method, float tmax) const {
        if (method == "uniform") {
            return 1.0f / s;
        }
        else { // method == "transmittance"
            return std::exp(-m_sigmaT*s);
        }
    }

    /*
    - infinite
    - homogeneous
        Tr(x,y) = e^(-sigmaS * norm(x-y))
    - isotropic
        PDF = 1 / 4pi
    - monochromatic
    - no surfaces
    - single scattering
        w NEE
        1 x area light (area sampling)
    > in-scattering (3rd part of VRE)
    */
    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray, const Intersection* _its) const {
        if (m_tracerType == "explicit") {
            return Li_explicit(scene, sampler, ray);
        }
        else { //if (m_tracerType == "implicit") {
            return 0.0f;
        }
    }

    Color3f Li_explicit(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        //float maxt = scene->getBoundingBox().getExtents().norm();
        float maxt = 100.0f;
        
        Point3f x = ray.o;
        Vector3f wi = ray.d;

        // sample distance (t)
        float t = distSample(m_distanceSampling, maxt, sampler);
        float pdfT = distPdf(t, m_distanceSampling, maxt);

        Intersection its;
        if (scene->rayIntersect(ray, its) && its.shape->isEmitter()) {
            float s = (x-its.p).norm();
            if(t >= s) {
                float pdfT = distPdf_failure(s, m_distanceSampling, maxt);
                Color3f Le = Tr(x,its.p) * its.shape->getEmitter()->eval(its, -ray.d) / pdfT;
                return Le;
            }
        }

        // calc point xt
        Point3f xt = x + t*wi;

        // sample light
        Emitter* em = scene->getEmitters()[0];
        Color3f Le_(0.0f);
        /*
        // - area
        Normal3f yN;
        Point3f xe = em->sample(sampler, yN);
        Vector3f wo = (xe - xt).normalized();

        Ray3f lightRay(xt, wo, Epsilon, maxt);
        Intersection itsLight;
        bool intersects = scene->rayIntersect(lightRay, itsLight);
        if (intersects && (itsLight.shape->isEmitter())) {
            float cosThetaY = std::max(0.0f, (-wo).dot(yN));
            if (cosThetaY > 0.0f) { // check for division by zero
                float pA = 1.0f / em->getShape()->getArea();
                float d2 = (xe - xt).squaredNorm();
                float pdf = d2 / cosThetaY * pA;

                Le_ = em->eval(itsLight, wo) / pdf;
            }
        }
        */

        // - solidangle
        Point3f xe; Normal3f yN; float pWi;
        Vector3f wo = em->sampleSolidAngle(sampler, xt, yN, pWi, xe);

        Ray3f lightRay(xt, wo, Epsilon, maxt);
        Intersection itsLight;
        bool intersects = scene->rayIntersect(lightRay, itsLight);
        if (intersects && (itsLight.shape->isEmitter())) {
            Le_ = itsLight.shape->getEmitter()->eval(itsLight, wo) / pWi;
        }

        // calc Li w NEE
        Color3f Li = Tr(xt, xe) * Le_;

        // calc Ls
        float fp = 1.0f / (4.0f * M_PI);
        Color3f Ls = fp * Li;

        Color3f Lr = Tr(x, xt) * m_sigmaS * Ls / pdfT;
        return Lr;
    }

    std::string toString() const {
        return tfm::format(
            "VolumePT[\n"
            " tracerType = %s\n"
            " termination = %s\n"
            " terminationProb = %d\n"
            " terminationBounds = %d\n"
            " distanceSampling = %s\n"
            " sigmaA = %d\n"
            " sigmaS = %d\n"
            " sigmaT = %d\n"
            "]",
            m_tracerType,
            m_termination,
            m_terminationProb,
            m_terminationBounds,
            m_distanceSampling,
            m_sigmaA,
            m_sigmaS,
            m_sigmaT
        );
    }

private:
    const std::string m_tracerType;
    const std::string m_termination;
    const float m_terminationProb;
    const int m_terminationBounds;
    const std::string m_distanceSampling;
    const float m_sigmaA;
    const float m_sigmaS;
    const float m_sigmaT;
};

NORI_REGISTER_CLASS(VolumePT, "volumePT");
NORI_NAMESPACE_END