
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
        m_samplingMethod(props.getString("samplingMethod", "area")),
		m_sampleCount(props.getInteger("sampleCount", 1)),
        m_emitterSamples(props.getInteger("emitterSamples", 1)),
        m_brdfSamples(props.getInteger("brdfSamples", 1))
	{
		
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
		/* Find the surface that is visible in the requested direction */
		Intersection its;
		if (!scene->rayIntersect(ray, its))
			return Color3f(0.0f);

		if (its.shape->isEmitter()) {
			return its.shape->getEmitter()->eval();
		}

		Normal3f n = its.shFrame.n;
		float maxt = scene->getBoundingBox().getExtents().norm();

		// monte carlo!
		Color3f Lr(0.0f);
		if (m_samplingMethod == "uniform") {
			
			for (int i = 0; i < m_sampleCount; ++i) {
				Vector3f d = Warp::squareToUniformHemisphere(sampler->next2D());
				float pdf = Warp::squareToUniformHemispherePdf(d);
				d = its.toWorld(d); // transform to world space so it aligns with the its
				d.normalize();
				Ray3f lightRay(its.p, d, Epsilon, maxt);

				Intersection itsLight;
				bool intersects = scene->rayIntersect(lightRay, itsLight);
				if (intersects && itsLight.shape->isEmitter()) {
					const Emitter* emitter = itsLight.shape->getEmitter();
					Color3f Le = emitter->eval();
					float cosTheta = std::max(0.0f, d.dot(n));
					nori::BSDFQueryRecord bRec(its.toLocal(d), its.toLocal(-ray.d), nori::ESolidAngle, its.toLocal(n));
					nori::Color3f brdfValue = its.shape->getBSDF()->eval(bRec);

					Lr += brdfValue * Le * cosTheta / pdf;
				}
			}
			Lr *= 1.0f / m_sampleCount;

		}
		else if (m_samplingMethod == "brdf") {

			for (int i = 0; i < m_sampleCount; ++i) {
				nori::BSDFQueryRecord bRec(its.toLocal(-ray.d), Vector3f(0.0f), nori::ESolidAngle, its.toLocal(n));
				nori::Color3f brdfValue = its.shape->getBSDF()->sample(bRec, sampler->next2D());
				Vector3f d = its.toWorld(bRec.wo); // transform to world space so it aligns with the its
				d.normalize();

				Ray3f lightRay(its.p, d, Epsilon, maxt);
				Intersection itsLight;
				bool intersects = scene->rayIntersect(lightRay, itsLight);
				if (intersects && itsLight.shape->isEmitter()) {
					const Emitter* emitter = itsLight.shape->getEmitter();
					Color3f Le = emitter->eval();
					Lr += brdfValue * Le;
				}
			}
			Lr *= 1.0f / m_sampleCount;

		}
		else if (m_samplingMethod == "area") {

			for (int i = 0; i < m_sampleCount; ++i) {
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

							nori::BSDFQueryRecord bRec(its.toLocal(d), its.toLocal(-ray.d), nori::ESolidAngle, its.toLocal(n));
							nori::Color3f brdfValue = its.shape->getBSDF()->eval(bRec);

							Lr += brdfValue * Le * cosTheta / pA;
						}
					}
				}
			}

			Lr *= 1.0f / m_sampleCount;

		}
		else if (m_samplingMethod == "solidangle") {
			
			for (int i = 0; i < m_sampleCount; ++i) {
				for (const Emitter* emitter : scene->getEmitters()) {
					const Shape* lightShape = emitter->getShape();
					if (lightShape) { // is area light?
						Normal3f yN;
						float pWi;
						Vector3f d = emitter->sampleSolidAngle(sampler, its.p, yN, pWi);
						
						Ray3f lightRay(its.p, d, Epsilon, maxt);
						Intersection itsLight;
						bool intersects = scene->rayIntersect(lightRay, itsLight);
						if (intersects && itsLight.shape->isEmitter()) {
							Color3f Le = itsLight.shape->getEmitter()->eval();
							float cosTheta = std::max(0.0f, d.dot(n));
							nori::BSDFQueryRecord bRec(its.toLocal(d), its.toLocal(-ray.d), nori::ESolidAngle, its.toLocal(n));
							nori::Color3f brdfValue = its.shape->getBSDF()->eval(bRec);

							Lr += brdfValue * Le * cosTheta / pWi;
						}
					}
				}
			}
			Lr *= 1.0f / m_sampleCount;

		}
		else if (m_samplingMethod == "mis") {

            // light sampling
			for (const Emitter* emitter : scene->getEmitters()) {
				const Shape* lightShape = emitter->getShape();
				if (lightShape) { // is area light?

					for (int i = 0; i < m_emitterSamples; ++i) {
						/*
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

							nori::BSDFQueryRecord bRec(its.toLocal(d), its.toLocal(-ray.d), nori::ESolidAngle, its.toLocal(n));
							nori::Color3f brdfValue = its.shape->getBSDF()->eval(bRec);

							float pdfBrdf = its.shape->getBSDF()->pdf(bRec);
							float weight = balanceHeuristic(m_emitterSamples, pA, m_brdfSamples, pdfBrdf);

                            Lr += brdfValue * Le * cosTheta * weight / (pA * m_emitterSamples);
						}
						*/

                        Normal3f yN;
                        float pWi;
                        Vector3f d = emitter->sampleSolidAngle(sampler, its.p, yN, pWi);

                        Ray3f lightRay(its.p, d, Epsilon, maxt);
                        Intersection itsLight;
                        bool intersects = scene->rayIntersect(lightRay, itsLight);
                        if (intersects && itsLight.shape->isEmitter()) {
                            Color3f Le = itsLight.shape->getEmitter()->eval();
                            float cosTheta = std::max(0.0f, d.dot(n));
                            nori::BSDFQueryRecord bRec(its.toLocal(d), its.toLocal(-ray.d), nori::ESolidAngle, its.toLocal(n));
                            nori::Color3f brdfValue = its.shape->getBSDF()->eval(bRec);

                            float pdfBrdf = its.shape->getBSDF()->pdf(bRec);
                            float weight = balanceHeuristic(m_emitterSamples, pWi, m_brdfSamples, pdfBrdf);

                            Lr += brdfValue * Le * cosTheta * weight / (pWi * m_emitterSamples);
                        }

					}

				}
			}

            // BRDF sampling
            for (int i = 0; i < m_brdfSamples; ++i) {
                nori::BSDFQueryRecord bRec(its.toLocal(-ray.d), Vector3f(0.0f), nori::ESolidAngle, its.toLocal(n));
                nori::Color3f brdfValue = its.shape->getBSDF()->sample(bRec, sampler->next2D());
                Vector3f d = its.toWorld(bRec.wo); // transform to world space so it aligns with the its
                d.normalize();

                Ray3f lightRay(its.p, d, Epsilon, maxt);
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

		return Lr;
	}

	float balanceHeuristic(int n1, float pdf1, int n2, float pdf2) const {
		return (n1 * pdf1) / (n1 * pdf1 + n2 * pdf2);
	}

	float powerHeuristic(int n1, float pdf1, int n2, float pdf2, float power) {
		return powf(n1 * pdf1, power) / pow(n1 * pdf1 + n2 * pdf2, power);
	}

	std::string toString() const {
		return tfm::format(
			"DirectIntegrator[\n"
			" samplingMethod = %d\n"
			" sampleCount = %d\n"
            " emitterSamples = %d\n"
            " brdfSamples = %d\n"
			"]",
			m_samplingMethod,
			m_sampleCount,
            m_emitterSamples,
            m_brdfSamples
		);
	}

private:
	const std::string m_samplingMethod;
	const int m_sampleCount;
    const int m_emitterSamples;
    const int m_brdfSamples;
};

NORI_REGISTER_CLASS(DirectIntegrator, "direct");
NORI_NAMESPACE_END