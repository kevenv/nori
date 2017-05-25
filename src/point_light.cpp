
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

class PointLight : public Emitter {
public:
	PointLight(const PropertyList &propList) {
		
	}

	Color3f eval() const override {
		return Color3f(0.0f);
	}

	/// Return a human-readable summary
	std::string toString() const {
		return tfm::format(
			"PointLight[\n");
	}

private:
	
};

NORI_REGISTER_CLASS(PointLight, "point");
NORI_NAMESPACE_END
