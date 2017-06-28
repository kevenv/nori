
#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Phong
 */
class Phong : public BSDF {
public:
    Phong(const PropertyList &propList):
		m_albedo(propList.getColor("albedo", Color3f(1.0f))),
		m_shininess(propList.getFloat("shininess", 1.0f)),
		m_diffuseReflectance(propList.getFloat("diffuseReflectance", 0.25f)),
		m_specularReflectance(propList.getFloat("specularReflectance", 1.0f))
	{

    }

    /// Evaluate the BRDF model
    Color3f eval(const BSDFQueryRecord &bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);
		
		Vector3f wr(2*bRec.N * (bRec.N.dot(bRec.wi)) - bRec.wi);
		float alpha = std::max(wr.dot(bRec.wo), 0.0f);
		float specular = (m_shininess + 2) / (2 * M_PI) * pow(alpha, m_shininess);
		return m_diffuseReflectance * m_albedo/M_PI + m_specularReflectance * specular;
    }

    /// Compute the density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
		return 0.0f; //todo
    }

    /// Draw a a sample from the BRDF model
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
        return Color3f(0.0f); //todo
    }

    bool isDiffuse() const {
        return false;
    }

    /// Return a human-readable summary
    std::string toString() const {
        return tfm::format(
            "Phong[\n"
			"  albedo = %s\n"
            "  shinyness = %d\n"
			"  diffuseReflectance = %d\n"
			"  specularReflectance = %d\n"
            "]",
			m_albedo.toString(),
			m_shininess,
			m_diffuseReflectance,
			m_specularReflectance
			);
    }

    EClassType getClassType() const { return EBSDF; }

private:
	const Color3f m_albedo;
	const float m_shininess;
	const float m_diffuseReflectance;
	const float m_specularReflectance;
};

NORI_REGISTER_CLASS(Phong, "phong");
NORI_NAMESPACE_END
