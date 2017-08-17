
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>

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
        m_currentPhotonCount(0)
	{
        m_photonMap.resize(m_photonCount);
	}

    void preprocess(const Scene *scene) override {
        generatePhotonMap(scene);
    }

    void generatePhotonMap(const Scene* scene) {
        std::unique_ptr<Sampler> sampler = scene->getSampler()->clone();

        while(!haveEnoughPhotons()) {
            const Emitter& em = chooseRandomEmitter(scene);
            Photon& p = emitPhotonFromEmitter(em, sampler.get());
            tracePhoton(p, scene, sampler.get());
        }

        for(int i = 0; i < m_photonCount; ++i) {
            m_photonMap[i].phi /= m_photonCount;
        }
    }

    bool haveEnoughPhotons() {
        return m_currentPhotonCount >= m_photonCount-1;
    }

    const Emitter& chooseRandomEmitter(const Scene* scene) {
        return *scene->getEmitters()[0]; //todo
    }

    Photon& emitPhotonFromEmitter(const Emitter& em, Sampler* sampler) {
        Photon p; // don't store 'starting' photons directly on the emitter

        // sample position
        Normal3f n; //todo: need -n?
        p.x = em.sample(sampler, n);
        float pdfX = 1.0f / em.getShape()->getArea();

        // sample direction
        Vector3f wLoc = Warp::squareToCosineHemisphere(sampler->next2D());
        float pdfW = Warp::squareToCosineHemispherePdf(wLoc);
        Frame N(n);
        p.w = N.toWorld(wLoc);

        // compute power
        Color3f Le = em.eval();
        float cosTheta = std::max(0.0f,p.w.dot(n));
        p.phi = 1/m_photonCount * Le * cosTheta / pdfX*pdfW;
    }

    void tracePhoton(const Photon& p, const Scene* scene, Sampler* sampler) {
        Intersection its;
        Ray3f ray(p.x, p.w);
        if(scene->rayIntersect(ray, its)) {
            //store photon if diffuse
            m_photonMap[m_currentPhotonCount].x = its.p;
            m_photonMap[m_currentPhotonCount].w = p.w;
            m_photonMap[m_currentPhotonCount].phi = p.phi;

            if(!haveEnoughPhotons()) {
                return; // no more photons todo: will that break things?
            }
            m_currentPhotonCount++;

            // scatter photon
            BSDFQueryRecord bRec(p.w, Vector3f(0.0f), ESolidAngle);
            Color3f brdfValue = its.shape->getBSDF()->sample(bRec, sampler->next2D());

            // compute new photon power
            Photon p_;
            p_.x = its.p;
            p_.w = bRec.wo;

            Normal3f n(its.shFrame.n);
            float cosTheta = std::max(0.0f,p_.w.dot(n));
            p_.phi = p.phi * cosTheta * brdfValue;

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
            p_.phi /= (Color3f(1.0f) - p_.phi); // photon absorption
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
		if (!scene->rayIntersect(ray, its))
			return Color3f(0.0f);

		if (its.shape->isEmitter()) {
			return its.shape->getEmitter()->eval();
		}

        Color3f Lr(0.0f);
		return Lr;
	}

	std::string toString() const override {
		return tfm::format(
			"PPM[\n"
            "photonCount = %d\n"
			"]",
            m_photonCount
		);
	}

private:
    const int m_photonCount;

    std::vector<Photon> m_photonMap;
    int m_currentPhotonCount;
};

NORI_REGISTER_CLASS(PPM, "ppm");
NORI_NAMESPACE_END