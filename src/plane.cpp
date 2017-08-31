
#include <nori/shape.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

/**
* \brief Plane
*/
class Plane : public Shape {
public:
	Plane(const PropertyList &propList):
		m_center(propList.getPoint("center", Point3f(0.0f))),
		m_width(propList.getFloat("width", 1.0f)),
		m_height(propList.getFloat("height", 1.0f)),
		m_normalZ(propList.getFloat("normalZ",1.0f))
	{
		m_bbox.expandBy(Vector3f( -m_width/2,  m_height/2, 0.0f) + m_center);
		m_bbox.expandBy(Vector3f(  m_width/2,  m_height/2, 0.0f) + m_center);
		m_bbox.expandBy(Vector3f(  m_width/2, -m_height/2, 0.0f) + m_center);
		m_bbox.expandBy(Vector3f( -m_width/2, -m_height/2, 0.0f) + m_center);

		m_normal = Vector3f(0.0f, 0.0f, m_normalZ); //todo: plane is hardcoded for Z=0
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
		return m_width * m_height;
	}

	virtual Point3f sample(Sampler* sampler, Normal3f& normal) const override {
		Point2f sample = sampler->next2D();
		Vector3f v(sample.x(), sample.y(), 0.0f);

		Point3f y = Vector3f(v.x() * m_width, v.y() * m_height, 0.0f) + m_center;

		normal = m_normal; //todo: normal = y - x; normal.normalize();
		return y;
	}

	virtual Vector3f sampleSolidAngle(Sampler* sampler, Point3f& x, Normal3f& normal, float& pWi) const override {
		throw NoriException("Unimplemented Plane::sampleSolidAngle() !!!");
	}

	virtual bool rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t) const override {
        //see: https://www.cl.cam.ac.uk/teaching/1999/AGraphHCI/SMAG/node2.html#eqn:vectplane

		Vector3f N(m_normal); //plane normal
		Point3f Q(m_center); //point on the plane
		Point3f E(ray.o);
		Vector3f D(ray.d);

		float nDotD = N.dot(D);
		if(nDotD <= 1e-6f) return false;

		t = N.dot(Q-E) / nDotD;
		if(t < 0.0f) return false;

        Point3f P(ray(t)); //in world

        return P.x() >= (-m_width/2 + m_center.x()) && P.x() <= (m_width/2 + m_center.x()) &&
               P.y() >= (-m_height/2 + m_center.y()) && P.y() <= (m_height/2 + m_center.y());
	}

	virtual void computeIntersectionInfo(uint32_t index, const Ray3f &ray, Intersection & its) const override {
		its.p = ray(its.t) - m_center;

		its.geoFrame = Frame(m_normal);
		its.shFrame = Frame(m_normal);
	}

	virtual std::string toString() const override {
		return tfm::format(
			"%s\n"
			"Plane[\n"
			"center = %s\n"
			"width = %d\n"
			"height = %d\n"
			"]",
			Shape::toString(),
			m_center.toString(),
			m_width,
			m_height
		);
	}

protected:
	Point3f m_center;
	float m_width;
	float m_height;

	float m_normalZ; //todo: set normal x,y,z

	Normal3f m_normal;
};

NORI_REGISTER_CLASS(Plane, "plane");
NORI_NAMESPACE_END
