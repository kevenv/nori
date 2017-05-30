
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/shape.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

class DirectIntegrator : public Integrator {
public:
	DirectIntegrator(const PropertyList &props) :
		m_sampleCount(props.getInteger("sampleCount", 1))
	{

	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
		/* Find the surface that is visible in the requested direction */
		Intersection its;
		if (!scene->rayIntersect(ray, its))
			return Color3f(0.0f);

		if (its.shape->isEmitter()) {
			return Color3f(1.0f);
		}

		Normal3f n = its.shFrame.n;
		float maxt = scene->getBoundingBox().getExtents().norm();

		// monte carlo!
		/*
		Color3f Lr(0.0f);
		for (int i = 0; i < m_sampleCount; ++i) {
			Vector3f d = Warp::squareToUniformHemisphere(sampler->next2D());
			d = its.toWorld(d); // transform to world space so it aligns with the its
			Ray3f lightRay(its.p, d, Epsilon, maxt);
			
			Intersection itsLight;
			bool intersects = scene->rayIntersect(lightRay, itsLight);
			if (intersects && itsLight.shape->isEmitter()) {
				const Emitter* emitter = itsLight.shape->getEmitter();
				Color3f Le = emitter->eval();
				float cosTheta = std::max(0.0f, d.dot(n));
				nori::BSDFQueryRecord bRec(d, its.toLocal(-ray.d), nori::ESolidAngle);
				nori::Color3f brdfValue = its.shape->getBSDF()->eval(bRec);

				Lr += brdfValue * Le * cosTheta;
			}
		}
		Lr *= 2.0f * M_PI / m_sampleCount;
		*/

		Color3f Lr(0.0f);
		for (int i = 0; i < m_sampleCount; ++i) {
			for (const Emitter* emitter : scene->getEmitters()) {
				const Shape* lightShape = emitter->getShape();
				if (lightShape) { // is area light?
					Color3f Le = emitter->eval();

					Point3f x = its.p;
					Vector3f v = Warp::squareToUniformSphere(sampler->next2D());
					Point3f y = v * 0.1f + Vector3f(1, -1.5, 1.75);
					Normal3f yN = Normal3f(y);
					Vector3f d = (y-x).normalized();
					
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

						Lr += brdfValue * Le * cosTheta / pA;
					}
				}
			}
		}
		
		Lr *= 1.0f / m_sampleCount;

		return Lr;
	}

	std::string toString() const {
		return tfm::format(
			"DirectIntegrator[\n"
			" sampleCount = %d\n"
			"]",
			m_sampleCount
		);
	}

private:
	const int m_sampleCount;
};

NORI_REGISTER_CLASS(DirectIntegrator, "direct");
NORI_NAMESPACE_END