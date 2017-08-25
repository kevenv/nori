
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/kdtree.h>

//#define USE_NAIVE_KNN

NORI_NAMESPACE_BEGIN

struct Photon {
    Point3f x;
    Vector3f w;
    Color3f phi;
};

struct PhotonDist {
    float dist;
    int idx;

    bool operator<(const PhotonDist& pDist) const {
        return dist < pDist.dist;
    }
};

class PPM : public Integrator {
public:
	PPM(const PropertyList &props):
        m_progressive(static_cast<bool>(props.getInteger("progressive",1))),
        m_photonCount(props.getInteger("photonCount", 100)),
        m_kPhotons(props.getInteger("kPhotons", 10)),
        m_radius2(props.getFloat("radius2", 10.0f)),
        m_samplesFinalGathering(props.getInteger("samplesFinalGathernig",10)),
        m_currentPhotonCount(0)
	{
        m_photonMap.resize(m_photonCount);
        m_KDTree.reserve(m_photonCount);
	}

    void preprocess(const Scene *scene) override {
        generatePhotonMap(scene);
    }

    virtual void beforeIteration(const Scene *scene, int iteration) override {
        preprocess(scene);
        updateRadius(iteration);
    }

    void updateRadius(int iteration) {
        float alpha = 2.0f/3.0f;
        m_radius2 *= ((float)iteration+alpha)/(iteration+1);
    }

    void generatePhotonMap(const Scene* scene) {
        m_KDTree.clear();
        m_currentPhotonCount = 0;

        std::unique_ptr<Sampler> samplerPtr = scene->getSampler()->clone();
        Sampler* sampler = samplerPtr.get();

        while(!haveEnoughPhotons()) {
            const Emitter& em = chooseRandomEmitter(scene);
            Photon p = emitPhotonFromEmitter(em, sampler);
            tracePhoton(p, scene, sampler);
        }

        for(int i = 0; i < m_photonCount; ++i) {
            m_photonMap[i].phi /= m_photonCount;
        }

        // build KD-Tree
        for(int i = 0; i < m_photonCount; ++i) {
            PhotonKDTreeNode node(m_photonMap[i].x, i);
            m_KDTree.push_back(node);
        }
        m_KDTree.build();
    }

    bool haveEnoughPhotons() {
        return m_currentPhotonCount >= m_photonCount-1;
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
        Vector3f wLoc = Warp::squareToUniformSphere(sampler->next2D()); //todo: sphere for sphere light, cosine for rect light
        float pdfW = Warp::squareToUniformSpherePdf(wLoc);
        Frame N(n);
        p.w = N.toWorld(wLoc);
        p.w.normalize();

        // compute power
        Color3f Le = em.eval();
        float cosTheta = std::max(0.0f,p.w.dot(n));
        p.phi = Le * cosTheta / (pdfX*pdfW);

        return p;
    }

    void tracePhoton(const Photon& p, const Scene* scene, Sampler* sampler) {
        float maxt = scene->getBoundingBox().getExtents().norm();

        Intersection its;
        Ray3f ray(p.x, p.w, Epsilon, maxt);
        if(scene->rayIntersect(ray, its)) {
            //store photon if diffuse
            if(!its.shape->isEmitter()) {
                m_photonMap[m_currentPhotonCount].x = its.p;
                m_photonMap[m_currentPhotonCount].w = p.w;
                m_photonMap[m_currentPhotonCount].phi = p.phi;

                if (haveEnoughPhotons()) {
                    return; // no more photons todo: will that break things?
                }
                m_currentPhotonCount++;
            }

            // scatter photon
            BSDFQueryRecord bRec(its.toLocal(p.w), Vector3f(0.0f), ESolidAngle);
            Color3f brdfValue = its.shape->getBSDF()->sample(bRec, sampler->next2D());
            float pdfBRDF = its.shape->getBSDF()->pdf(bRec);

            // compute new photon power
            Photon p_;
            p_.x = its.p;
            p_.w = its.toWorld(bRec.wo);

            Normal3f n(its.shFrame.n);
            float cosTheta = std::max(0.0f,p_.w.dot(n));
            p_.phi = p.phi * cosTheta * brdfValue / pdfBRDF;

            // continue scatter photon
            if(survivedRR(p, p_, sampler->next1D())) {
                tracePhoton(p_, scene, sampler);
            }
        }
    }

    bool survivedRR(const Photon& p, Photon& p_, float rand) {
        float phi = getLuminance(p.phi);
        float phi_ = getLuminance(p_.phi);

        float rr = 1 - std::min(1.0f, phi_ / phi);
        if(rand < rr) {
            return false; // photon dies, life is cruel for photons... :(
        }
        else {
            p_.phi /= (Color3f(1.0f) - Color3f(rr)); // photon absorption
            return true; // photon survived! :D
        }
    }

    /// Return the luminance (assuming the color value is expressed in linear sRGB)
    inline float getLuminance(Color3f c) const {
        return c[0] * 0.212671f + c[1] * 0.715160f + c[2] * 0.072169f;
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

        // final gathering
        if(m_samplesFinalGathering > 0) {
            Normal3f n = its.shFrame.n;
            float maxt = scene->getBoundingBox().getExtents().norm();

            Color3f Lr(0.0f);
            for (int i = 0; i < m_samplesFinalGathering; ++i) {
                Vector3f d = Warp::squareToUniformHemisphere(sampler->next2D());
                float pdf = Warp::squareToUniformHemispherePdf(d);
                d = its.toWorld(d); // transform to world space so it aligns with the its
                d.normalize();
                Ray3f gatherRay(its.p, d, Epsilon, maxt);

                Intersection itsGather;
                bool intersects = scene->rayIntersect(gatherRay, itsGather);
                if (intersects && !itsGather.shape->isEmitter()) {
                    Color3f Le = computeLrFromDensityEstimation(gatherRay, itsGather);
                    float cosTheta = std::max(0.0f, d.dot(n));
                    nori::BSDFQueryRecord bRec(its.toLocal(d), its.toLocal(-ray.d), nori::ESolidAngle);
                    nori::Color3f brdfValue = its.shape->getBSDF()->eval(bRec);

                    Lr += brdfValue * Le * cosTheta / pdf;
                }
            }
            Lr *= 1.0f / m_samplesFinalGathering;

            return Lr;
        }
        else {
            return computeLrFromDensityEstimation(ray, its);
        }

        //photonMapViewer(ray, its);
	}

    Color3f computeLrFromDensityEstimation(const Ray3f& ray, const Intersection& its) const {
        if(m_progressive) {
            // get k nearest photons
            std::vector<unsigned int> results;
            m_KDTree.search(its.p, std::sqrt(m_radius2), results);
            int kPhotons = results.size();
            if(kPhotons <= 0) return 0.0f;

            std::sort(results.begin(), results.end()); //ascending order

            std::vector<Photon> nearestPhotons(kPhotons);
            for(int i = 0; i < kPhotons; ++i) {
                unsigned int idxKD = results[i];
                PhotonMapIdx idx = m_KDTree[idxKD].getData();
                nearestPhotons[i] = m_photonMap[idx];
            }

            return densityEstimation(nearestPhotons, nearestPhotons.size(), m_radius2, ray, its);
        }
        else {
            // get k nearest photons
#ifdef USE_NAIVE_KNN
            std::vector<PhotonDist> dist(m_photonCount);
            for(int i = 0; i < m_photonCount; ++i) {
                Vector3f d = its.p - m_photonMap[i].x;
                dist[i].dist = d.norm();
                dist[i].idx = i;
            }

            std::sort(dist.begin(), dist.end()); //ascending order

            std::vector<Photon> nearestPhotons(m_kPhotons);
            for(int i = 0; i < m_kPhotons; ++i) {
                int idx = dist[i].idx;
                nearestPhotons[i] = m_photonMap[idx];
            }
            float radius2 = dist[m_kPhotons-1].dist; radius2*=radius2; // distance to the k-th photon
#else
            std::vector<PointKDTree<PhotonKDTreeNode>::SearchResult> results(m_kPhotons+1); // + 1 extra for nnSearch impl
            m_KDTree.nnSearch(its.p, m_kPhotons, &results[0]);

            std::sort(results.begin(), results.end()); //ascending order

            std::vector<Photon> nearestPhotons(m_kPhotons);
            for(int i = 0; i < m_kPhotons; ++i) {
                unsigned int idxKD = results[i].index;
                PhotonMapIdx idx = m_KDTree[idxKD].getData();
                nearestPhotons[i] = m_photonMap[idx];
            }
            float radius2 = results[m_kPhotons-1].distSquared; // distance to the k-th photon

            return densityEstimation(nearestPhotons, nearestPhotons.size()-1 /* ignore the k-th photon */, radius2, ray, its);
#endif
        }
    }

    Color3f densityEstimation(const std::vector<Photon>& nearestPhotons, int k, float radius2,
                             const Ray3f& ray, const Intersection& its) const {
        // compute radiance estimate from k nearest photons
        Color3f Lr(0.0f);
        for(int i = 0; i < k; ++i) {
            BSDFQueryRecord bRec(its.toLocal(m_photonMap[i].w), its.toLocal(-ray.d), ESolidAngle);
            Color3f brdfValue = its.shape->getBSDF()->eval(bRec);
            Lr += brdfValue * nearestPhotons[i].phi / (M_PI*radius2);
        }
        return Lr;
    }

    Color3f photonMapViewer(const Ray3f& ray, const Intersection& its) {
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
                BSDFQueryRecord bRec(its.toLocal(m_photonMap[i].w), its.toLocal(-ray.d), ESolidAngle);
                Color3f brdfValue = its.shape->getBSDF()->eval(bRec);
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
			"]",
            m_photonCount,
            m_kPhotons
		);
	}

private:
    const bool m_progressive;
    const int m_photonCount;
    const int m_kPhotons;
    const int m_samplesFinalGathering;

    std::vector<Photon> m_photonMap;
    int m_currentPhotonCount;
    typedef unsigned int PhotonMapIdx;
    typedef GenericKDTreeNode<Vector3f, PhotonMapIdx> PhotonKDTreeNode;
    PointKDTree<PhotonKDTreeNode> m_KDTree;

    float m_radius2;
};

NORI_REGISTER_CLASS(PPM, "ppm");
NORI_NAMESPACE_END