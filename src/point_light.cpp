
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

class PointLight : public Emitter {
public:
    PointLight(const PropertyList &propList) {
        
    }

    virtual Color3f eval(const Intersection& its, const Vector3f& d) const override {
        throw NoriException("Unimplemented PointLight::eval() !!!");
    }

    virtual Color3f evalPosition() const override {
        throw NoriException("Unimplemented PointLight::evalPosition() !!!");
    }

    virtual Point3f sample(Sampler* sampler, Normal3f& normal) const override {
        throw NoriException("Unimplemented PointLight::sample() !!!");
    }

    virtual Vector3f sampleSolidAngle(Sampler* sampler, Point3f& x, Normal3f& normal, float& pWi) const override {
        throw NoriException("Unimplemented PointLight::sampleSolidAngle() !!!");
    }

    /// Return a human-readable summary
    std::string toString() const override {
        return tfm::format(
            "PointLight[\n");
    }

private:
    
};

NORI_REGISTER_CLASS(PointLight, "point");
NORI_NAMESPACE_END
