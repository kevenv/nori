
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/shape.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <pcg32.h>

NORI_NAMESPACE_BEGIN

pcg32 s_random;

class PathTracer : public Integrator {
public:
	PathTracer(const PropertyList &props) :
		m_tracerType(props.getString("tracerType", "explicit")),
		m_termination(props.getString("termination", "russian-roulette")),
		m_terminationProb(props.getFloat("terminationProb", 0.2f)),
		m_directSampling(props.getString("directSampling", "area")),
		m_indirectSampling(props.getString("indirectSampling", "cosine"))
	{

	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
		if (m_tracerType == "explicit") {
			return Li_explicit(scene, sampler, ray, 0);
		}
		else if (m_tracerType == "implicit") {
			return Li_implicit(scene, sampler, ray, 0);
		}
	}

	Color3f Li_explicit(const Scene *scene, Sampler *sampler, const Ray3f &ray, int bounds) const {
		if (m_termination == "russian-roulette") {
			if (s_random.nextFloat() <= m_terminationProb) return Color3f(0.0f);
		}
		else if(m_termination == "path-depth") {
			if (bounds > 15) return Color3f(0.0f);
		}

		Intersection its;
		if (!scene->rayIntersect(ray, its))
			return Color3f(0.0f);

		Normal3f n = its.shFrame.n;
		float maxt = scene->getBoundingBox().getExtents().norm();

		// direct illumination
		Color3f L_dir(0.0f);
		for (const Emitter* emitter : scene->getEmitters()) {
			const Shape* lightShape = emitter->getShape();
			if (lightShape) { // is area light?
				Color3f Le = emitter->eval();

				Normal3f yN;
				Point3f x = its.p;
				Point3f y = emitter->sample(sampler, yN);
				Vector3f d = (y - x).normalized();

				Ray3f lightRay(x, d, Epsilon, maxt);
				Intersection itsLight;
				bool intersects = scene->rayIntersect(lightRay, itsLight);
				if (intersects && itsLight.shape->isEmitter()) {
					float cosTheta_i = std::max(0.0f, d.dot(n));
					float cosTheta_o = std::max(0.0f, d.dot(yN));
					float cosTheta = (cosTheta_i * cosTheta_o) / ((x - y).squaredNorm());

					float pA = 1.0f / lightShape->getArea();

					nori::BSDFQueryRecord bRec(its.toLocal(d), its.toLocal(-ray.d), nori::ESolidAngle);
					nori::Color3f brdfValue = its.shape->getBSDF()->eval(bRec);

					L_dir += brdfValue * Le * cosTheta / pA;
				}
			}
		}

		// indirect illumination

		// cast a random ray
		Vector3f d = Warp::squareToCosineHemisphere(sampler->next2D());
		d = its.toWorld(d); // transform to world space so it aligns with the its
		d.normalize();
		Ray3f traceRay(its.p, d, Epsilon, maxt);

		float cosTheta, pWi;
		if (m_indirectSampling == "cosine") {
			cosTheta = 1.0f;
			pWi = INV_PI;
		}
		else if (m_indirectSampling == "uniform") {
			cosTheta = std::max(0.0f, d.dot(n));
			pWi = INV_TWOPI;
		}
		nori::BSDFQueryRecord bRec(d, its.toLocal(-ray.d), nori::ESolidAngle);
		nori::Color3f fr = its.shape->getBSDF()->eval(bRec);
		Color3f L_ind = fr * Li_explicit(scene, sampler, traceRay, bounds++) * cosTheta / pWi;
		if (m_termination == "russian-roulette") {
			L_ind /= (1 - m_terminationProb);
		}

		Color3f Le(0.0f);
		return Le + L_dir + L_ind;
	}

	Color3f Li_implicit(const Scene *scene, Sampler *sampler, const Ray3f &ray, int bounds) const {
		if (m_termination == "russian-roulette") {
			if (s_random.nextFloat() <= m_terminationProb) return Color3f(0.0f);
		}
		else if (m_termination == "path-depth") {
			if (bounds > 15) return Color3f(0.0f);
		}

		Intersection its;
		if (!scene->rayIntersect(ray, its))
			return Color3f(0.0f);

		Normal3f n = its.shFrame.n;
		float maxt = scene->getBoundingBox().getExtents().norm();

		// cast a random ray
		Vector3f d = Warp::squareToCosineHemisphere(sampler->next2D());
		d = its.toWorld(d); // transform to world space so it aligns with the its
		d.normalize();
		Ray3f traceRay(its.p, d, Epsilon, maxt);

		nori::Color3f L(0.0f);
		if (its.shape->isEmitter()) {
			//direct illumination
			const Emitter* emitter = its.shape->getEmitter();
			L = emitter->eval();
		}
		else {
			//indirect illumination
			L = Li_implicit(scene, sampler, traceRay, bounds++);
		}

		float cosTheta, pWi;
		if (m_indirectSampling == "cosine") {
			cosTheta = 1.0f;
			pWi = INV_PI;
		}
		else if (m_indirectSampling == "uniform") {
			cosTheta = std::max(0.0f, d.dot(n));
			pWi = INV_TWOPI;
		}
		nori::BSDFQueryRecord bRec(d, its.toLocal(-ray.d), nori::ESolidAngle);
		nori::Color3f fr = its.shape->getBSDF()->eval(bRec);
		Color3f Le(0.0f);
		Color3f Lr = fr * L * cosTheta / pWi;
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
			" directSampling = %s\n"
			" indirectSampling = %s\n"
			"]",
			m_tracerType,
			m_termination,
			m_terminationProb,
			m_directSampling,
			m_indirectSampling
		);
	}

private:
	const std::string m_tracerType;
	const std::string m_termination;
	const float m_terminationProb;
	const std::string m_directSampling;
	const std::string m_indirectSampling;
};

NORI_REGISTER_CLASS(PathTracer, "path_tracer");
NORI_NAMESPACE_END