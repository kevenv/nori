
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

class AreaLight : public Emitter {
public:
	AreaLight(const PropertyList &propList) :
		m_radiance(propList.getColor("radiance",Color3f(1.0f)))
	{

	}

	Color3f eval() const override {
		return m_radiance;
	}

	/// Return a human-readable summary
	std::string toString() const {
		return tfm::format(
			"AreaLight[\n");
	}

private:
	Color3f m_radiance;
};

NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END
