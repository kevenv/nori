
#include <nori/shape.h>
#include <nori/common.h>

NORI_NAMESPACE_BEGIN

/**
* \brief Sphere
*/
class Sphere : public Shape {
public:
	Sphere(const PropertyList &propList):
		m_center(propList.getPoint("center", Point3f(0.0f))),
		m_radius(propList.getFloat("radius", 1.0f)),
		m_zMin(propList.getFloat("zMin", -m_radius)),
		m_zMax(propList.getFloat("zMax", m_radius)),
		m_thetaMin(propList.getFloat("thetaMin", 0.0f)),
		m_thetaMax(propList.getFloat("thetaMax", M_PI)),
		phiMax(propList.getFloat("phiMax", 2*M_PI))
	{
		m_bbox.expandBy(m_center - Vector3f(m_radius));
		m_bbox.expandBy(m_center + Vector3f(m_radius));
	}

	virtual uint32_t getPrimitiveCount() const override {
		return 1;
	}

	virtual BoundingBox3f getBoundingBox(uint32_t index) const override {
		return m_bbox;
	}

	virtual Point3f getCentroid(uint32_t index) const override {
		return m_center;
	}

	virtual float getArea() const override {
		return 4 * M_PI * m_radius * m_radius;
	}

	virtual bool rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t) const override {
		// transform Ray (world space) to object space

		// compute quadratic sphere coefficients
		float a = ray.d.dot(ray.d);
		float b = 2.0f * ray.d.dot(ray.o - m_center);
		float c = (ray.o - m_center).dot(ray.o - m_center) - m_radius * m_radius;

		// t = solve quadratic equation
		float t0, t1;
		if (!solveQuadratic(a, b, c, t0, t1)) {
			return false;
		}
		// check quadratic t0,t1 for nearest intersection
		const float epsillon = 1e-8f;
		if (t0 + epsillon > ray.maxt || t1 - epsillon <= ray.mint) {
			return false;
		}
		t = t0;
		if (t - epsillon <= 0.0f) {
			t = t1;
			if (t + epsillon > ray.maxt) {
				return false;
			}
		}

		// compute sphere hit position
		nori::Point3f pHit = ray(t);
		if (pHit.x() == 0.0f && pHit.y() == 0.0) {
			pHit.x() = 1e-5f * m_radius;
		}
		float phi = std::atan2(pHit.y(), pHit.x());
		if (phi < 0) {
			phi += 2 * M_PI;
		}

		// test sphere intersection against clipping params
		if ((m_zMin > -m_radius && pHit.z() < m_zMin) || (m_zMax < m_radius && pHit.z() > m_zMax) || (phi > phiMax)) {
			if (t == t1) return false;
			if (t1 + epsillon > ray.maxt) return false;
			t = t1;
			pHit = ray(t);
			if (pHit.x() == 0.0f && pHit.y() == 0.0) {
				pHit.x() = 1e-5f * m_radius;
			}
			phi = std::atan2(pHit.y(), pHit.x());
			if (phi < 0) {
				phi += 2 * M_PI;
			}
			if ((m_zMin > -m_radius && pHit.z() < m_zMin) || (m_zMax < m_radius && pHit.z() > m_zMax) || (phi > phiMax)) {
				return false;
			}
		}

		// compute parametric representation of pHit
		u = phi / phiMax;
		float theta = std::acos(clamp(pHit.z() / m_radius, -1.0f, 1.0f));
		v = (theta - m_thetaMin) / (m_thetaMax - m_thetaMin);

		return true;
	}

	virtual void computeIntersectionInfo(uint32_t index, const Ray3f &ray, Intersection & its) const override {
		float u, v, t;
		rayIntersect(index, ray, u, v, t);
		its.p = ray(t);

		Vector3f l = its.p - m_center;
		Vector3f n = l.normalized();
		its.geoFrame = Frame(n);
		its.shFrame = Frame(n);
	}

	virtual std::string toString() const override {
		return tfm::format(
			"%s\n"
			"Sphere[\n"
			"center = %s"
			"radius = %d"
			"]",
			Shape::toString(),
			m_center.toString(),
			m_radius
		);
	}

protected:
	Point3f m_center;
	float m_radius;
	float m_zMin, m_zMax;
	float m_thetaMin, m_thetaMax, phiMax;
};

NORI_REGISTER_CLASS(Sphere, "sphere");
NORI_NAMESPACE_END
