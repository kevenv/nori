
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/shape.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/direct.h>

NORI_NAMESPACE_BEGIN

struct DistSample {
    float t;
    float pdf;
};

class VolumePT : public Integrator {
public:
    VolumePT(const PropertyList &props) :
        m_tracerType(props.getString("tracerType", "explicit")),
        m_termination(props.getString("termination", "russian-roulette")),
        m_terminationProb(props.getFloat("terminationProb", 0.2f)),
        m_terminationBounds(props.getInteger("terminationBounds", 15)),
        m_distanceSampling(props.getString("distanceSampling", "transmittance")),
        m_emitterSampling(props.getString("emitterSampling", "solidangle")),
        m_sigmaA(props.getFloat("sigmaA", 0.1f)),
        m_sigmaS(props.getFloat("sigmaS", 0.5f)),
        m_sigmaT(m_sigmaA + m_sigmaS)
    {

    }

    float Tr(const nori::Point3f& x, const nori::Point3f& y) const {
        return std::exp( -m_sigmaT * (x - y).norm() );
    }

    DistSample distSample(const std::string& method, float tmax, Sampler* sampler, const Ray3f& ray, const Point3f& xe) const {
        DistSample s;

        if (method == "uniform") {
            s.t = sampler->next1D() * tmax;
            s.pdf = 1.0f / tmax;
        }
        else if(method == "transmittance") {
            s.t = -std::log(1 - sampler->next1D()) / m_sigmaT;
            s.pdf = m_sigmaT * std::exp(-m_sigmaT*s.t);
        }
        else { //if(method == "equi-angular")
            float delta = (xe - ray.o).dot(ray.d);
            float D = (xe - (ray.o + delta*ray.d)).norm();
            float thetaA = std::atan((0.0 - delta)/D);
            float thetaB = std::atan((tmax - delta)/D);
            float ei = sampler->next1D();
            s.t = D * std::tan((1-ei)*thetaA + ei*thetaB);
            s.pdf = D / ((thetaB - thetaA)*(D*D+s.t*s.t));
            s.t += delta;
        }

        return s;
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
    - multiple scattering
        w NEE
        1 x area light (area sampling)
    > in-scattering (3rd part of VRE)
    */
    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray, const Intersection* _its) const {
        if (m_tracerType == "explicit") {
            return Li_explicit(scene, sampler, ray, 0);
        }
        else { //if (m_tracerType == "implicit") {
            return Li_implicit(scene, sampler, ray);
        }
    }

    Color3f Li_explicit(const Scene *scene, Sampler *sampler, const Ray3f &ray, int bounds) const {
        if (m_termination == "russian-roulette") {
            if (sampler->next1D() <= m_terminationProb || (bounds >= m_terminationBounds)) return Color3f(0.0f);
        }
        else if(m_termination == "path-depth") {
            if (bounds >= m_terminationBounds) return Color3f(0.0f);
        }

        //float maxt = scene->getBoundingBox().getExtents().norm();
        float maxt = 100.0f;
        
        Point3f x = ray.o;
        Vector3f wi = ray.d;

        // sample light position (xe)
        Emitter* em = scene->getEmitters()[0];
        Normal3f yN;
        // point light = position
        // sphere light = point on surface
        Point3f xe = em->sample(sampler, yN);

        // sample distance (t)
        DistSample dSample = distSample(m_distanceSampling, maxt, sampler, ray, xe);
        float t = dSample.t;
        float pdfT = dSample.pdf;

        // Emission - direct light hit (explicit)
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
        Color3f Le_(0.0f);

        if(!em->isDeltaLight()) {
            if(m_emitterSampling == "area") {
                Vector3f wo = (xe - xt).normalized();

                Ray3f lightRay(xt, wo, Epsilon, maxt);
                Intersection itsLight;
                bool intersects = scene->rayIntersect(lightRay, itsLight);
                if (intersects && (itsLight.shape->isEmitter())) {
                    float cosThetaY = std::max(0.0f, (-wo).dot(yN));
                    if (cosThetaY > 0.0f) { // check for division by zero
                        float pA = 1.0f / itsLight.shape->getArea();
                        float d2 = (xe - xt).squaredNorm();
                        float pdf = d2 / cosThetaY * pA;

                        Le_ = itsLight.shape->getEmitter()->eval(itsLight, wo) / pdf;
                    }
                    //xe = itsLight.p;
                }
            }
            else { //if(m_emitterSampling == "solidangle")
                float pWi;
                Vector3f wo = em->sampleSolidAngle(sampler, xt, yN, pWi, xe);

                Ray3f lightRay(xt, wo, Epsilon, maxt);
                Intersection itsLight;
                bool intersects = scene->rayIntersect(lightRay, itsLight);
                if (intersects && (itsLight.shape->isEmitter())) {
                    xe = itsLight.p;
                    Le_ = itsLight.shape->getEmitter()->eval(itsLight, wo) / pWi;
                }
            }
        }
        else {
            // - point light
            Vector3f wo = (xt - xe).normalized();
            Ray3f lightRay(xt, wo, Epsilon, maxt);
            Intersection itsLight;
            bool intersects = scene->rayIntersect(lightRay, itsLight);
            if (!intersects) {
                Le_ = em->eval(itsLight, wo) / (xt - xe).squaredNorm();
            }
        }

        // calc Li w NEE
        Color3f Li = Tr(xt, xe) * Le_;

        // calc Ldir
        float fp = 1.0f / (4.0f * M_PI);
        Color3f Ldir = fp * Li;

        // calc Lind
        // sample fp
        Vector3f wo_ = Warp::squareToUniformSphere(sampler->next2D());
        float pdfW = Warp::squareToUniformSpherePdf(wo_);
        //wo_ = Frame().toWorld(wo_);
        //wo_.normalize();
        Ray3f traceRay(xt, wo_, Epsilon, maxt);

        Color3f Lind = fp * Li_explicit(scene, sampler, traceRay, ++bounds) / pdfW;

        Color3f Ls = Ldir + Lind;

        Color3f Lr = Tr(x, xt) * m_sigmaS * Ls / pdfT;

        if (m_termination == "russian-roulette") {
            Lr /= (1 - m_terminationProb);
        }

        return Lr;
    }

    Color3f Li_implicit(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        //float maxt = scene->getBoundingBox().getExtents().norm();
        float maxt = 100.0f;

        Point3f x = ray.o;
        Vector3f wi = ray.d;

        // get light position (xe)
        Emitter* em = scene->getEmitters()[0];
        Point3f xe(0.0f);

        // sample distance (t)
        DistSample dSample = distSample(m_distanceSampling, maxt, sampler, ray, xe);
        float t = dSample.t;
        float pdfT = dSample.pdf;

        // Emission - direct light hit (explicit)
        Intersection its;
        if (scene->rayIntersect(ray, its) && its.shape->isEmitter()) {
            float s = (x - its.p).norm();
            if (t >= s) {
                float pdfT = distPdf_failure(s, m_distanceSampling, maxt);
                Color3f Le = Tr(x, its.p) * its.shape->getEmitter()->eval(its, -ray.d) / pdfT;
                return Le;
            }
        }

        // calc point xt
        Point3f xt = x + t*wi;

        // direct illumination
        Color3f Lr(0.0f);

        Vector3f wo = Warp::squareToUniformSphere(sampler->next2D());
        float pdf = Warp::squareToUniformSpherePdf(wo);

        Ray3f lightRay(xt, wo, Epsilon, maxt);
        Intersection itsLight;
        bool intersects = scene->rayIntersect(lightRay, itsLight);
        if (intersects && (itsLight.shape->isEmitter())) {
            Point3f xe = itsLight.p;
            Color3f Le = itsLight.shape->getEmitter()->eval(itsLight, wo) / pdf;

            // calc Li
            Color3f Li = Tr(xt, xe) * Le;

            // calc Ls
            float fp = 1.0f / (4.0f * M_PI);
            Color3f Ls = fp * Li;

            Lr = Tr(x, xt) * m_sigmaS * Ls / pdfT;
        }

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
            " emitterSampling = %s\n"
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
    const std::string m_emitterSampling;
    const float m_sigmaA;
    const float m_sigmaS;
    const float m_sigmaT;
};

NORI_REGISTER_CLASS(VolumePT, "volumePT");
NORI_NAMESPACE_END