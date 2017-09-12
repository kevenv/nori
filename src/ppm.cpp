
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/kdtree.h>

NORI_NAMESPACE_BEGIN

struct Photon {
    Point3f x;
    Vector3f w;
    Color3f phi;
};

class PPM : public Integrator {
public:
    PPM(const PropertyList &props):
        m_photonCount(props.getInteger("photonCount", 100)),
        m_kPhotons(props.getInteger("kPhotons", 10)),
        m_radius2(props.getFloat("radius2", 10.0f)),
        m_samplesFG(props.getInteger("samplesFG", 100)),
        m_samplesDI(props.getInteger("samplesDI", 50)),
        m_knnMethodStr(props.getString("knnMethod", "radius")),
        m_knnMethod(KNN_METHOD_RADIUS),
        m_currentPhotonCount(0),
        m_emittedPhotonCount(0)
    {
        m_progressive = static_cast<bool>(props.getInteger("progressive",1));
        m_iterations = props.getInteger("iterations", 1);

        if(m_knnMethodStr == "radius" || m_progressive) {
            m_knnMethod = KNN_METHOD_RADIUS;
        } else if(m_knnMethodStr == "photons") {
            m_knnMethod = KNN_METHOD_PHOTONS;
        }

        m_photonMap.resize(m_photonCount);
        m_KDTree.reserve(m_photonCount);

        if(m_kPhotons > m_photonCount) {
            throw NoriException("PPM: kPhotons > photonCount...");
        }
    }

    void preprocess(const Scene *scene) override {
        generatePhotonMap(scene);
    }

    virtual void beforeIteration(const Scene *scene, int iteration) override {
        preprocess(scene);
        updateRadius(iteration);
    }

    void updateRadius(int iteration) {
        float alpha = iteration == 0 ? 1.0f : 2.0f/3.0f;
        m_radius2 *= ((float)iteration + alpha) / (iteration + 1);
    }

    void generatePhotonMap(const Scene* scene) {
        m_KDTree.clear();
        m_currentPhotonCount = 0;
        m_emittedPhotonCount = 0;

        std::unique_ptr<Sampler> samplerPtr = scene->getSampler()->clone();
        Sampler* sampler = samplerPtr.get();

        while(!haveEnoughPhotons()) {
            const Emitter& em = chooseRandomEmitter(scene);
            Photon p = emitPhotonFromEmitter(em, sampler);
            tracePhoton(p, scene, sampler, 0);
        }

        for(int i = 0; i < m_photonCount; ++i) {
            m_photonMap[i].phi /= m_emittedPhotonCount;
        }

        // build KD-Tree
        for(int i = 0; i < m_photonCount; ++i) {
            PhotonKDTreeNode node(m_photonMap[i].x, i);
            m_KDTree.push_back(node);
        }
        m_KDTree.build();
    }

    bool haveEnoughPhotons() {
        return m_currentPhotonCount >= m_photonCount;
    }

    const Emitter& chooseRandomEmitter(const Scene* scene) {
        return *scene->getEmitters()[0]; //todo
    }

    Photon emitPhotonFromEmitter(const Emitter& em, Sampler* sampler) {
        Photon p; // don't store 'starting' photons directly on the emitter

        // sample position
        Normal3f n;
        p.x = em.sample(sampler, n);
        //n = -n; //todo: need -n? w rect light
        float pdfX = 1.0f / em.getShape()->getArea();

        // sample direction
        Vector3f wLoc = Warp::squareToUniformHemisphere(sampler->next2D()); //todo: hemisphere for hemisphere light, cosine for rect light
        float pdfW = Warp::squareToUniformHemispherePdf(wLoc);
        Frame N(n);
        p.w = N.toWorld(wLoc);
        p.w.normalize();

        // compute power
        Color3f Le = em.eval();
        float cosTheta = std::max(0.0f,p.w.dot(n));
        p.phi = Le * cosTheta / (pdfX*pdfW);

        m_emittedPhotonCount++;

        return p;
    }

    void tracePhoton(Photon& p, const Scene* scene, Sampler* sampler, int bounds = 0) {
        if (haveEnoughPhotons()) {
            return; // no more photons todo: will that break things?
        }

        float maxt = scene->getBoundingBox().getExtents().norm();

        Intersection its;
        Ray3f ray(p.x, p.w, Epsilon, maxt);
        if(!scene->rayIntersect(ray, its)) {
            return;
        }

        // store photon if diffuse
        if(its.shape->isEmitter()) {
            return;
        }

        // don't estimate DI with PM
        if (bounds > 0) {
            m_photonMap[m_currentPhotonCount].x = its.p;
            m_photonMap[m_currentPhotonCount].w = p.w;
            m_photonMap[m_currentPhotonCount].phi = p.phi;
            m_currentPhotonCount++;
        }

        // compute new photon power
		//wi,wo
        BSDFQueryRecord bRec(its.toLocal(-p.w), Vector3f(0.0f), ESolidAngle);
        Color3f brdfValue = its.shape->getBSDF()->sample(bRec, sampler->next2D()); // = fr * cos / pdf

        Photon p_;
        p_.x = its.p;
        p_.w = its.toWorld(bRec.wo);
        p_.phi = p.phi * brdfValue;

        // scatter photon
        if(survivedRR(p, p_, sampler->next1D())) {
            tracePhoton(p_, scene, sampler, ++bounds);
        }

        /*
        if(bounds < 12 && !p_.phi.isBlack()) {
            tracePhoton(p_, scene, sampler, ++bounds);
        }
        */
    }

    bool survivedRR(const Photon& p, Photon& p_, float rand) {
        float phi = p.phi.getLuminance();
        float phi_ = p_.phi.getLuminance();

        float rr = 1 - std::min(1.0f, phi_ / phi);
        if(rand < rr) {
            return false; // photon dies, life is cruel for photons... :(
        }
        else {
            p_.phi /= (Color3f(1.0f) - Color3f(rr)); // photon absorption
            return true; // photon survived! :D
        }
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const override {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its)) {
            return 0.0f;
        }
        if (its.shape->isEmitter()) {
            return its.shape->getEmitter()->eval();
        }

        //return photonMapViewer(ray, its);

        Normal3f n = its.shFrame.n;
        float maxt = scene->getBoundingBox().getExtents().norm();

        // DI
        Color3f Ld(0.0f);
        if (m_samplesDI > 0) {
            for (int i = 0; i < m_samplesDI; ++i) {
                for (const Emitter* emitter : scene->getEmitters()) {
                    const Shape* lightShape = emitter->getShape();
                    if (lightShape) { // is area light?
                        Normal3f yN;
                        float pWi;
                        Vector3f wo = emitter->sampleSolidAngle(sampler, its.p, yN, pWi);

                        Ray3f lightRay(its.p, wo, Epsilon, maxt);
                        Intersection itsLight;
                        bool intersects = scene->rayIntersect(lightRay, itsLight);
                        if (intersects && itsLight.shape->isEmitter()) {
                            Color3f Le = itsLight.shape->getEmitter()->eval();
							//wi,wo
                            nori::BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(wo), nori::ESolidAngle, its.toLocal(n));
                            nori::Color3f brdfValue = its.shape->getBSDF()->eval(bRec); // BRDF * cosTheta

                            Ld += brdfValue * Le / pWi;
                        }
                    }
                }
            }
            Ld *= 1.0f / m_samplesDI;
        }

        // GI using PM w final gathering (FG)
        Color3f Li(0.0f);
        if(m_samplesFG > 0) {
            for (int i = 0; i < m_samplesFG; ++i) {
                Vector3f wo = Warp::squareToCosineHemisphere(sampler->next2D());
                float pdf = Warp::squareToCosineHemispherePdf(wo);
                wo = its.toWorld(wo); // transform to world space so it aligns with the its
                wo.normalize();
                Ray3f gatherRay(its.p, wo, Epsilon, maxt);

                Intersection itsGather;
                bool intersects = scene->rayIntersect(gatherRay, itsGather);
                if (intersects && !itsGather.shape->isEmitter()) {
                    Color3f Le = computeLrFromDensityEstimation(gatherRay, itsGather);
					//wi,wo
                    nori::BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(wo), nori::ESolidAngle);
                    nori::Color3f brdfValue = its.shape->getBSDF()->eval(bRec); // BRDF * cosTheta

                    Li += brdfValue * Le / pdf;
                }
            }
            Li *= 1.0f / m_samplesFG;
            return Ld + Li;
        }
        else {
            return Ld + computeLrFromDensityEstimation(ray, its);
        }
    }

    Color3f computeLrFromDensityEstimation(const Ray3f& ray, const Intersection& its) const {
        std::vector<Photon> nearestPhotons;
        if(m_knnMethod == KNN_METHOD_RADIUS) {
            float radius2 = m_radius2;
            if(nnSearch(its.p, m_photonCount, radius2, nearestPhotons)) {
                return densityEstimation(nearestPhotons, nearestPhotons.size(), radius2, ray, its);
            }
        }
        else if(m_knnMethod == KNN_METHOD_PHOTONS) {
            float radius2 = std::numeric_limits<float>::infinity(); // will be changed by nnSearch()
            if(nnSearch(its.p, m_kPhotons, radius2, nearestPhotons)) {
                return densityEstimation(nearestPhotons, nearestPhotons.size()-1 /* ignore the k-th photon */, radius2, ray, its);
            }
        }
        return 0.0f;
    }

    bool nnSearch(const Point3f& p, int kPhotonsMax, float& radius2, std::vector<Photon>& nearestPhotons) const {
        std::vector<PointKDTree<PhotonKDTreeNode>::SearchResult> results(kPhotonsMax+1); // + 1 extra for nnSearch impl
        int kPhotons = m_KDTree.nnSearch(p, radius2, kPhotonsMax, &results[0]); // radius2 will be changed by nnSearch()
        if(kPhotons <= 0) return false;

        nearestPhotons.resize(kPhotons);
        for(int i = 0; i < kPhotons; ++i) {
            unsigned int idxKD = results[i].index;
            PhotonMapIdx idx = m_KDTree[idxKD].getData();
            nearestPhotons[i] = m_photonMap[idx];
        }
        return true;
    }

    Color3f densityEstimation(const std::vector<Photon>& nearestPhotons, int k, float radius2,
                             const Ray3f& ray, const Intersection& its) const {
        // compute radiance estimate from k nearest photons
        Color3f Lr(0.0f);
        for(int i = 0; i < k; ++i) {
			//wi,wo
            BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(-m_photonMap[i].w), ESolidAngle);
            Color3f brdfValue = its.shape->getBSDF()->eval(bRec) /* BRDF * cosTheta */ / Frame::cosTheta(bRec.wo);
            Lr += brdfValue * nearestPhotons[i].phi / (M_PI*radius2);
        }
        return Lr;
    }

    Color3f photonMapViewer(const Ray3f& ray, const Intersection& its) const {
        for(int i = 0; i < m_photonCount; ++i) {
            // ray-point intersection
            Point3f P = m_photonMap[i].x;

            // point is in ray if PO || TO
            Vector3f PO(P - ray.o);
            Vector3f TO(ray.o + 1*ray.d - ray.o);

            // A || B and + <=> A*B == |A|*|B|
            bool intersects = std::abs(PO.dot(TO)-(PO.norm() * TO.norm())) <= 1e-6f;
            //float dot = PO.dot(TO);
            //bool intersects = std::abs(dot*dot - PO.squaredNorm()*TO.squaredNorm()) <= 1e-6f;
            if(intersects) {
                return Color3f(1.0f,0.0f,0.0f);
            }
        }
        return 0.0f;
    }

    std::string toString() const override {
        return tfm::format(
            "PPM[\n"
            "photonCount = %d\n"
            "kPhotons = %d\n"
            "radius2 = %d\n"
            "samplesFG = %d\n"
            "samplesDI = %d\n"
            "knnMethod = %s\n"
            "]",
            m_photonCount,
            m_kPhotons,
            m_radius2,
            m_samplesFG,
            m_samplesDI,
            m_knnMethodStr
        );
    }

private:
    enum KNN_METHOD {
        KNN_METHOD_RADIUS,
        KNN_METHOD_PHOTONS
    };

    const int m_photonCount;
    const int m_kPhotons;
    float m_radius2;
    const int m_samplesFG;
    const int m_samplesDI;
    const std::string m_knnMethodStr;
    KNN_METHOD m_knnMethod;

    std::vector<Photon> m_photonMap;
    int m_currentPhotonCount;
    typedef unsigned int PhotonMapIdx;
    typedef GenericKDTreeNode<Vector3f, PhotonMapIdx> PhotonKDTreeNode;
    PointKDTree<PhotonKDTreeNode> m_KDTree;

    int m_emittedPhotonCount;
};

NORI_REGISTER_CLASS(PPM, "ppm");
NORI_NAMESPACE_END