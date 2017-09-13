#include <nori/direct.h>

#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

DirectIntegrator::DirectIntegrator(const PropertyList &props) :
    m_samplingMethod(props.getString("samplingMethod", "area")),
    m_emitterSamples(props.getInteger("emitterSamples", 1)),
    m_brdfSamples(props.getInteger("brdfSamples", 1))
{

}

DirectIntegrator::DirectIntegrator(const std::string& samplingMethod, int emitterSamples, int brdfSamples):
    m_samplingMethod(samplingMethod),
    m_emitterSamples(emitterSamples),
    m_brdfSamples(brdfSamples)
{

}

Color3f DirectIntegrator::Li(const Scene *scene, Sampler *sampler, const Ray3f &ray, const Intersection* _its) const {
    /* Find the surface that is visible in the requested direction */
    Intersection its;
    if (!_its) {
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        if (its.shape->isEmitter()) {
            return its.shape->getEmitter()->eval();
        }
    }
    else {
        its = *_its; // DirectIntegrator used by another integrator
    }

    Normal3f n = its.shFrame.n;
    float maxt = scene->getBoundingBox().getExtents().norm();

    // monte carlo!
    Color3f Lr(0.0f);
    if (m_samplingMethod == "uniform") {

        if (m_emitterSamples <= 0) return 0.0f;
            
        for (int i = 0; i < m_emitterSamples; ++i) {
            Vector3f wo = Warp::squareToUniformHemisphere(sampler->next2D());
            float pdf = Warp::squareToUniformHemispherePdf(wo);
            wo = its.toWorld(wo); // transform to world space so it aligns with the its
            wo.normalize();
            Ray3f lightRay(its.p, wo, Epsilon, maxt);

            Intersection itsLight;
            bool intersects = scene->rayIntersect(lightRay, itsLight);
            if (intersects && itsLight.shape->isEmitter()) {
                const Emitter* emitter = itsLight.shape->getEmitter();
                Color3f Le = emitter->eval();
                //wi,wo
                nori::BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(wo), nori::ESolidAngle);
                nori::Color3f brdfValue = its.shape->getBSDF()->eval(bRec); // BRDF * cosTheta

                Lr += brdfValue * Le / pdf;
            }
        }
        Lr *= 1.0f / m_emitterSamples;

    }
    else if (m_samplingMethod == "brdf") {

        if (m_emitterSamples <= 0) return 0.0f;

        for (int i = 0; i < m_emitterSamples; ++i) {
            nori::BSDFQueryRecord bRec(its.toLocal(-ray.d), Vector3f(0.0f), nori::ESolidAngle);
            nori::Color3f brdfValue = its.shape->getBSDF()->sample(bRec, sampler->next2D()); // BRDF / pdf * cosTheta
            Vector3f wo = its.toWorld(bRec.wo); // transform to world space so it aligns with the its
            wo.normalize();

            Ray3f lightRay(its.p, wo, Epsilon, maxt);
            Intersection itsLight;
            bool intersects = scene->rayIntersect(lightRay, itsLight);
            if (intersects && itsLight.shape->isEmitter()) {
                const Emitter* emitter = itsLight.shape->getEmitter();
                Color3f Le = emitter->eval();
                Lr += brdfValue * Le;
            }
        }
        Lr *= 1.0f / m_emitterSamples;

    }
    else if (m_samplingMethod == "area") {

        if (m_emitterSamples <= 0) return 0.0f;

        for (int i = 0; i < m_emitterSamples; ++i) {
            for (const Emitter* emitter : scene->getEmitters()) {
                const Shape* lightShape = emitter->getShape();
                if (lightShape) { // is area light?
                    Color3f Le = emitter->eval();

                    Normal3f yN;
                    Point3f x = its.p;
                    Point3f y = emitter->sample(sampler, yN);
                    Vector3f wo = (y - x).normalized();

                    Ray3f lightRay(x, wo, Epsilon, maxt);
                    Intersection itsLight;
                    bool intersects = scene->rayIntersect(lightRay, itsLight);
                    if (intersects && (itsLight.shape->isEmitter() && itsLight.shape == lightShape)) {
                        float pA = 1.0f / lightShape->getArea();
                        float d2 = (y - x).squaredNorm();
                        float cosThetaY = std::max(0.0f, (-wo).dot(yN)); //todo: division by zero
                        float pdf = d2 / cosThetaY * pA;

                        //wi,wo
                        nori::BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(wo), nori::ESolidAngle);
                        nori::Color3f brdfValue = its.shape->getBSDF()->eval(bRec); // BRDF * cosTheta

                        Lr += brdfValue * Le / pdf;
                    }
                }
            }
        }
        Lr *= 1.0f / m_emitterSamples;

    }
    else if (m_samplingMethod == "solidangle") {
        
        if (m_emitterSamples <= 0) return 0.0f;

        for (int i = 0; i < m_emitterSamples; ++i) {
            for (const Emitter* emitter : scene->getEmitters()) {
                const Shape* lightShape = emitter->getShape();
                if (lightShape) { // is area light?
                    Normal3f yN;
                    float pWi;
                    Vector3f wo = emitter->sampleSolidAngle(sampler, its.p, yN, pWi);
                        
                    Ray3f lightRay(its.p, wo, Epsilon, maxt);
                    Intersection itsLight;
                    bool intersects = scene->rayIntersect(lightRay, itsLight);
                    if (intersects && (itsLight.shape->isEmitter() && itsLight.shape == lightShape)) {
                        Color3f Le = itsLight.shape->getEmitter()->eval();
                        //wi,wo
                        nori::BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(wo), nori::ESolidAngle);
                        nori::Color3f brdfValue = its.shape->getBSDF()->eval(bRec); // BRDF * cosTheta

                        Lr += brdfValue * Le / pWi;
                    }
                }
            }
        }
        Lr *= 1.0f / m_emitterSamples;

    }
    else if (m_samplingMethod == "mis") {

        // light sampling
        if (m_emitterSamples > 0) {
            for (const Emitter* emitter : scene->getEmitters()) {
                const Shape* lightShape = emitter->getShape();
                if (lightShape) { // is area light?

                    for (int i = 0; i < m_emitterSamples; ++i) {
                        Normal3f yN;
                        float pWi;
                        Vector3f wo = emitter->sampleSolidAngle(sampler, its.p, yN, pWi);

                        Ray3f lightRay(its.p, wo, Epsilon, maxt);
                        Intersection itsLight;
                        bool intersects = scene->rayIntersect(lightRay, itsLight);
                        if (intersects && (itsLight.shape->isEmitter() && itsLight.shape == lightShape)) {
                            Color3f Le = itsLight.shape->getEmitter()->eval();
                            //wi,wo
                            nori::BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(wo), nori::ESolidAngle);
                            nori::Color3f brdfValue = its.shape->getBSDF()->eval(bRec); // BRDF * cosTheta

                            float pdfBrdf = its.shape->getBSDF()->pdf(bRec);
                            float weight = balanceHeuristic(m_emitterSamples, pWi, m_brdfSamples, pdfBrdf);

                            Lr += brdfValue * Le * weight / (pWi * m_emitterSamples);
                        }
                    }

                }
            }
        }

        // BRDF sampling
        if (m_brdfSamples > 0) {
            for (int i = 0; i < m_brdfSamples; ++i) {
                nori::BSDFQueryRecord bRec(its.toLocal(-ray.d), Vector3f(0.0f), nori::ESolidAngle);
                nori::Color3f brdfValue = its.shape->getBSDF()->sample(bRec, sampler->next2D()); // BRDF / pdf * cosTheta
                Vector3f wo = its.toWorld(bRec.wo); // transform to world space so it aligns with the its
                wo.normalize();

                Ray3f lightRay(its.p, wo, Epsilon, maxt);
                Intersection itsLight;
                bool intersects = scene->rayIntersect(lightRay, itsLight);
                if (intersects && itsLight.shape->isEmitter()) {
                    const Emitter* emitter = itsLight.shape->getEmitter();
                    Color3f Le = emitter->eval();

                    float pdfEmitter; Normal3f yN;
                    emitter->sampleSolidAngle(sampler, its.p, yN, pdfEmitter);
                    float weight = balanceHeuristic(m_brdfSamples, its.shape->getBSDF()->pdf(bRec), m_emitterSamples, pdfEmitter);

                    Lr += brdfValue * Le * weight / m_brdfSamples;
                }
            }
        }

    }

    return Lr;
}

float DirectIntegrator::balanceHeuristic(int n1, float pdf1, int n2, float pdf2) const {
    return (n1 * pdf1) / (n1 * pdf1 + n2 * pdf2);
}

float DirectIntegrator::powerHeuristic(int n1, float pdf1, int n2, float pdf2, float power) const {
    return powf(n1 * pdf1, power) / pow(n1 * pdf1 + n2 * pdf2, power);
}

std::string DirectIntegrator::toString() const {
    return tfm::format(
        "DirectIntegrator[\n"
        " samplingMethod = %d\n"
        " emitterSamples = %d\n"
        " brdfSamples = %d\n"
        "]",
        m_samplingMethod,
        m_emitterSamples,
        m_brdfSamples
    );
}

NORI_REGISTER_CLASS(DirectIntegrator, "direct");
NORI_NAMESPACE_END