
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

class PointLight : public Emitter {
public:
    PointLight(const PropertyList &propList) :
        m_position(propList.getPoint("position", Point3f(0.0f,0.0f,0.0f))),
        m_intensity(propList.getColor("intensity", Color3f(1.0f,1.0f,1.0f)))
    {
        m_deltaLight = true;
    }

    virtual Color3f eval(const Intersection& its, const Vector3f& d) const override {
        return m_intensity;
    }

    virtual Point3f sample(Sampler* sampler, Normal3f& normal) const override {
        return m_position;
    }

    virtual Vector3f sampleSolidAngle(Sampler* sampler, Point3f& x, Normal3f& normal, float& pWi, Point3f& y) const override {
        throw NoriException("Unimplemented PointLight::sampleSolidAngle() !!!");
    }

    /// Return a human-readable summary
    std::string toString() const override {
        return tfm::format(
            "PointLight[\n"
            " position = %s\n",
            " intensity = %s\n",
            m_position.toString(),
            m_intensity.toString()
        );
    }

private:
    const Point3f m_position;
    const Color3f m_intensity;
};

NORI_REGISTER_CLASS(PointLight, "point");
NORI_NAMESPACE_END
