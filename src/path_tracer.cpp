
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
        else if (m_tracerType == "implicit") {
            return Li_implicit(scene, sampler, ray, 0);
        }
        else if (m_tracerType == "implicit-iter") {
            return Li_implicit_iter(scene, sampler, ray);
        }
        else { //if (m_tracerType == "implicit-exp") {
            return Li_implicit_exp(scene, sampler, ray, 0);
        }
    }

    Color3f Li_explicit(const Scene *scene, Sampler *sampler, const Ray3f &ray, int bounds) const {
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);
        
        // emission - direct light hit (explicit)
        if (bounds == 0 && its.shape->isEmitter()) {
            return Color3f(its.shape->getEmitter()->eval(its, -ray.d));
        }

        // stop indirect illumination if hit a light
        if (bounds > 1 && its.shape->isEmitter()) {
            return Color3f(0.0f);
        }

        if (m_termination == "russian-roulette") {
            if (sampler->next1D() <= m_terminationProb) return Color3f(0.0f);
        }
        else if (m_termination == "path-depth") {
            if (bounds >= m_terminationBounds) return Color3f(0.0f);
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
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        // direct illumination / emission (implicit)
        if (its.shape->isEmitter()) {
            return Color3f(its.shape->getEmitter()->eval(its, -ray.d));
        }

        if (m_termination == "russian-roulette") {
            if (sampler->next1D() <= m_terminationProb) return Color3f(0.0f);
        }
        else if (m_termination == "path-depth") {
            if (bounds >= m_terminationBounds) return Color3f(0.0f);
        }

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
        Ray3f traceRay(its.p, wo, Epsilon, maxt);

        // indirect illumination
        Color3f L = Li_implicit(scene, sampler, traceRay, ++bounds);

        //wi,wo
        nori::BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(wo), nori::ESolidAngle);
        nori::Color3f fr = its.shape->getBSDF()->eval(bRec); // BRDF * cosTheta
        Color3f Lr = fr * L / pdf;
        if (m_termination == "russian-roulette") {
            Lr /= (1 - m_terminationProb);
        }

        Color3f Le(0.0f);
        return Le + Lr;
    }

    Color3f Li_implicit_iter(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Color3f T(1.0f); // throughput
        Color3f Lr(0.0f);
        Ray3f wi(ray);

        Intersection its;
        if(!scene->rayIntersect(ray, its)) return 0.0f;

        // explicit emission
        if(its.shape->isEmitter()) {
            Color3f Le = its.shape->getEmitter()->eval(its, -wi.d);
            return Le;
        }

        if(m_terminationBounds == 0) {
            return 0.0f;
        }

        float maxt = scene->getBoundingBox().getExtents().norm();

        int bounds = 0;
        while(bounds < m_terminationBounds || m_termination == "russian-roulette") {
            /*
            // implicit emission

            if(!intersects) return 0.0f;

            if(its.shape->isEmitter()) {
                Color3f Le = its.shape->getEmitter()->eval(its, -wi.d);
                return Le;
            }

            if(m_terminationBounds == 1 && i == 0) {
                return 0.0f;
            }
            */

            Normal3f n = its.shFrame.n;

            // cast a random ray
            Vector3f woDir;
            float pdf;
            if (m_indirectSampling == "cosine") {
                woDir = Warp::squareToCosineHemisphere(sampler->next2D());
                pdf = Warp::squareToCosineHemispherePdf(woDir);
            }
            else if (m_indirectSampling == "uniform") {
                woDir = Warp::squareToUniformHemisphere(sampler->next2D());
                pdf = Warp::squareToUniformHemispherePdf(woDir);
            }
            woDir = its.toWorld(woDir); // transform to world space so it aligns with the its
            woDir.normalize();
            Ray3f wo(its.p, woDir, Epsilon, maxt);

            //wi,wo
            nori::BSDFQueryRecord bRec(its.toLocal(-wi.d), its.toLocal(wo.d), nori::ESolidAngle);
            nori::Color3f fr = its.shape->getBSDF()->eval(bRec); // BRDF * cosTheta

            if(!scene->rayIntersect(wo, its)) return 0.0f;

            if(its.shape->isEmitter()) {
                // direct illumination
                Color3f Le = its.shape->getEmitter()->eval(its, wo.d);
                Lr += T * fr * Le / pdf;
                break;
            }
            else {
                // indirect illumination
                T *= fr / pdf;
            }

            wi = wo;

            // RR
            if (m_termination == "russian-roulette") {
                T /= (1 - m_terminationProb); // do it on the last bounce, not the full path

                if (sampler->next1D() <= m_terminationProb) return Lr;
            }

            bounds++;
        }

        return Lr;
    }

    Color3f Li_implicit_exp(const Scene *scene, Sampler *sampler, const Ray3f &ray, int bounds) const {
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        // emission - direct light hit (explicit)
        if (bounds == 0 && its.shape->isEmitter()) {
            return Color3f(its.shape->getEmitter()->eval(its, -ray.d));
        }

        if (m_termination == "russian-roulette") {
            if (sampler->next1D() <= m_terminationProb) return Color3f(0.0f);
        }
        if (m_termination == "path-depth") {
            if (bounds >= m_terminationBounds) return Color3f(0.0f);
        }

        Normal3f n = its.shFrame.n;
        float maxt = scene->getBoundingBox().getExtents().norm();

        Color3f Lr(0.0f);

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
        Ray3f traceRay(its.p, wo, Epsilon, maxt);

        Intersection itsTrace;
        if (scene->rayIntersect(traceRay, itsTrace)) {
            Color3f L(0.0f);
            if (itsTrace.shape->isEmitter()) {
                //direct illumination
                const Emitter* emitter = itsTrace.shape->getEmitter();
                L = emitter->eval(itsTrace, wo);
            }
            else {
                //indirect illumination
                L = Li_implicit(scene, sampler, traceRay, ++bounds);
            }

            //wi,wo
            nori::BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(wo), nori::ESolidAngle);
            nori::Color3f fr = its.shape->getBSDF()->eval(bRec); // BRDF * cosTheta
            Lr = fr * L / pdf;
            if (m_termination == "russian-roulette") {
                Lr /= (1 - m_terminationProb);
            }
        }

        Color3f Le(0.0f);
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