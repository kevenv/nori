
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/shape.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/direct.h>

NORI_NAMESPACE_BEGIN

class PathTracer : public Integrator {
public:
    PathTracer(const PropertyList &props) :
        m_tracerType(props.getString("tracerType", "explicit")),
        m_termination(props.getString("termination", "russian-roulette")),
        m_terminationProb(props.getFloat("terminationProb", 0.2f)),
        m_terminationBounds(props.getInteger("terminationBounds", 15)),
        m_directSampling(props.getString("directSampling", "area")),
        m_indirectSampling(props.getString("indirectSampling", "cosine")),
        m_directIntegrator(m_directSampling, 1, 0)
    {

    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray, const Intersection* _its) const {
        if (m_tracerType == "explicit") {
            return Li_explicit(scene, sampler, ray, 0);
        }
        else { //if (m_tracerType == "implicit") {
            return Li_implicit(scene, sampler, ray, 0);
        }
    }

    Color3f Li_explicit(const Scene *scene, Sampler *sampler, const Ray3f &ray, int bounds) const {
        if (m_termination == "russian-roulette") {
            if (sampler->next1D() <= m_terminationProb) return Color3f(0.0f);
        }
        else if(m_termination == "path-depth") {
            if (bounds > m_terminationBounds) return Color3f(0.0f);
        }

        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        if(its.shape->isEmitter()) {
            // stop indirect illumination if hit a light
            if(bounds > 0) {
                return Color3f(0.0f);
            }
            // color lights
            else {
                return Color3f(its.shape->getEmitter()->eval(its, -ray.d));
            }
        }

        Normal3f n = its.shFrame.n;
        float maxt = scene->getBoundingBox().getExtents().norm();

        // direct illumination
        Color3f L_dir = m_directIntegrator.Li(scene, sampler, ray, &its);
        
        // indirect illumination

        // cast a random ray
        Vector3f wo;
        float pdf;
        Ray3f traceRay;
        // avoid double counting DI in ID
        Intersection itsTmp;
        do {
            if (m_indirectSampling == "cosine") {
                wo = Warp::squareToCosineHemisphere(sampler->next2D());
                pdf = Warp::squareToCosineHemispherePdf(wo);
            }
            else if (m_indirectSampling == "uniform") {
                wo = Warp::squareToUniformHemisphere(sampler->next2D());
                pdf = Warp::squareToUniformHemispherePdf(wo);
            }
            wo = its.toWorld(wo); // transform to world space so it aligns with the its
            wo.normalize();
            traceRay = Ray3f(its.p, wo, Epsilon, maxt);
        } while(scene->rayIntersect(traceRay, itsTmp) && itsTmp.shape->isEmitter());

        //wi,wo
        nori::BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(wo), nori::ESolidAngle);
        nori::Color3f fr = its.shape->getBSDF()->eval(bRec); // BRDF * cosTheta
        Color3f L_ind = fr * Li_explicit(scene, sampler, traceRay, ++bounds) / pdf;
        if (m_termination == "russian-roulette") {
            L_ind /= (1 - m_terminationProb);
        }

        Color3f Le(0.0f);
        return Le + L_dir + L_ind;
    }

    Color3f Li_implicit(const Scene *scene, Sampler *sampler, const Ray3f &ray, int bounds) const {
        if (m_termination == "russian-roulette") {
            if (sampler->next1D() <= m_terminationProb) return Color3f(0.0f);
        }
        else if (m_termination == "path-depth") {
            if (bounds > m_terminationBounds) return Color3f(0.0f);
        }

        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        Normal3f n = its.shFrame.n;
        float maxt = scene->getBoundingBox().getExtents().norm();

        // cast a random ray
        Vector3f wo;
        float pdf;
        if (m_indirectSampling == "cosine") {
            wo = Warp::squareToCosineHemisphere(sampler->next2D());
            pdf = Warp::squareToCosineHemispherePdf(wo);
        }
        else if (m_indirectSampling == "uniform") {
            wo = Warp::squareToUniformHemisphere(sampler->next2D());
            pdf = Warp::squareToUniformHemispherePdf(wo);
        }
        wo = its.toWorld(wo); // transform to world space so it aligns with the its
        wo.normalize();
        
        nori::Color3f L(0.0f);
        if (its.shape->isEmitter()) {
            //direct illumination
            const Emitter* emitter = its.shape->getEmitter();
            L = emitter->eval(its, wo);
        }
        else {
            //indirect illumination
            Ray3f traceRay(its.p, wo, Epsilon, maxt);
            L = Li_implicit(scene, sampler, traceRay, ++bounds);
        }

        //wi,wo
        nori::BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(wo), nori::ESolidAngle);
        nori::Color3f fr = its.shape->getBSDF()->eval(bRec); // BRDF * cosTheta
        Color3f Le(0.0f);
        Color3f Lr = fr * L / pdf;
        if (m_termination == "russian-roulette") {
            Lr /= (1 - m_terminationProb);
        }
        return Le + Lr;
    }

    std::string toString() const {
        return tfm::format(
            "PathTracer[\n"
            " tracerType = %s\n"
            " termination = %s\n"
            " terminationProb = %d\n"
            " terminationBounds = %d\n"
            " directSampling = %s\n"
            " indirectSampling = %s\n"
            "]",
            m_tracerType,
            m_termination,
            m_terminationProb,
            m_terminationBounds,
            m_directSampling,
            m_indirectSampling
        );
    }

private:
    const std::string m_tracerType;
    const std::string m_termination;
    const float m_terminationProb;
    const int m_terminationBounds;
    const std::string m_directSampling;
    const std::string m_indirectSampling;

    DirectIntegrator m_directIntegrator;
};

NORI_REGISTER_CLASS(PathTracer, "path_tracer");
NORI_NAMESPACE_END