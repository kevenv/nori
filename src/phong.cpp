
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
		m_shininess(propList.getInteger("shininess", 1)),
		m_diffuseReflectance(propList.getColor("diffuseReflectance", Color3f(1.0f,1.0f,1.0f))),
		m_specularReflectance(propList.getColor("specularReflectance", Color3f(1.0f, 1.0f, 1.0f)))
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
		return (m_diffuseReflectance / M_PI) + m_specularReflectance * specular;
    }

    /// Compute the density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
		/* return zero if the measure is wrong, or when 
		queried for illumination on the backside */
		if (bRec.measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) <= 0
			|| Frame::cosTheta(bRec.wo) <= 0)
			return 0.0f;

		float theta = sphericalCoordinates(bRec.wo).x();
		return (m_shininess + 2)/(2*M_PI) * pow(cos(theta), m_shininess);
    }

    /// Draw a a sample from the BRDF model
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
        if (Frame::cosTheta(bRec.wi) <= 0)
            return Color3f(0.0f);

		bRec.measure = ESolidAngle;

		float theta = acos(pow( (1 - sample.x()), 1.0f/(m_shininess+2) ));
		float phi = 2*M_PI*sample.y();
		bRec.wo = sphericalDirection(theta, phi);

        /* Relative index of refraction: no change */
        bRec.eta = 1.0f;

		return eval(bRec) / pdf(bRec) * Frame::cosTheta(bRec.wo);
    }

    bool isDiffuse() const {
        return false;
    }

    /// Return a human-readable summary
    std::string toString() const {
        return tfm::format(
            "Phong[\n"
            "  shininess = %d\n"
			"  diffuseReflectance = %s\n"
			"  specularReflectance = %s\n"
            "]",
			m_shininess,
			m_diffuseReflectance.toString(),
			m_specularReflectance.toString()
			);
    }

    EClassType getClassType() const { return EBSDF; }

private:
	const int m_shininess;
	const Color3f m_diffuseReflectance;
	const Color3f m_specularReflectance;
};

NORI_REGISTER_CLASS(Phong, "phong");
NORI_NAMESPACE_END
