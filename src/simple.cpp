
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class SimpleIntegrator : public Integrator {
public:
	SimpleIntegrator(const PropertyList &props) {
		/* No parameters this time */
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
		/* Find the surface that is visible in the requested direction */
		Intersection its;
		if (!scene->rayIntersect(ray, its))
			return Color3f(0.0f);

		nori::Vector3f light(3.0f, -1.0f, 2.0f);
		light.normalize();
		Ray3f shadowRay(its.p, light);
		bool visible = !scene->rayIntersect(shadowRay);
		if (visible) {
			nori::BSDFQueryRecord bRec(its.toLocal(light), its.toLocal(-ray.d), nori::ESolidAngle);
			nori::Color3f color = its.shape->getBSDF()->eval(bRec);
			Normal3f n = its.shFrame.n;
			float L = light.dot(n);
			L = std::max(L, 0.0f);
			return L * color;
		}
		else {
			return Color3f(0.0f);
		}
	}

	std::string toString() const {
		return "SimpleIntegrator[]";
	}
};

NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END