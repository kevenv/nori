
#include <nori/emitter.h>
#include <nori/shape.h>

NORI_NAMESPACE_BEGIN

class AreaLight : public Emitter {
public:
    AreaLight(const PropertyList &propList) :
        m_radiance(propList.getColor("radiance",Color3f(1.0f)))
    {

    }

    virtual Color3f eval(const Intersection& its, const Vector3f& d) const override {
        // careful when modifying this, the its in PPM is fake!
        return m_radiance;
    }

    virtual Point3f sample(Sampler* sampler, Normal3f& normal) const override {
        return m_shape->sample(sampler, normal);
    }

    virtual Vector3f sampleSolidAngle(Sampler* sampler, Point3f& x, Normal3f& normal, float& pWi, Point3f& y) const override {
        return m_shape->sampleSolidAngle(sampler, x, normal, pWi, y);
    }

    /// Return a human-readable summary
    std::string toString() const override {
        return tfm::format(
            "AreaLight[\n");
    }

private:
    Color3f m_radiance;
};

NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END
