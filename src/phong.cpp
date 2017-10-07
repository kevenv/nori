
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
        m_shininess(propList.getFloat("shininess", 1.0f)),
        m_diffuseReflectance(propList.getColor("diffuseReflectance", Color3f(1.0f,1.0f,1.0f))),
        m_specularReflectance(propList.getColor("specularReflectance", Color3f(1.0f, 1.0f, 1.0f)))
    {
        // ensure energy conservation
        // kd + ks <= 1
        const float max = 1.0f;
        Color3f tmp(m_diffuseReflectance + m_specularReflectance);
        float actualMax = tmp.max();
        float scale = 1.0f;
        if(actualMax > max) {
            scale = 0.99f * (max / actualMax);
        }
        m_diffuseReflectance *= scale;
        m_specularReflectance *= scale;

        float specAvg = m_specularReflectance.getLuminance();
        float diffAvg = m_diffuseReflectance.getLuminance();
        m_specularSamplingWeight = specAvg / (specAvg + diffAvg);
    }

    /// Reflection in local coordinates
    inline Vector3f reflect(const Vector3f &wi) const {
        //Same as : Vector3f wr(2*bRec.N * (bRec.N.dot(bRec.wi)) - bRec.wi);
        return Vector3f(-wi.x(), -wi.y(), wi.z());
    }

    /// Evaluate the BRDF model
    Color3f eval(const BSDFQueryRecord &bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);

        float alpha = bRec.wo.dot(reflect(bRec.wi));
        float specular = 0.0f;
        if(alpha > 0) {
            specular = (m_shininess + 2) / (2 * M_PI) * pow(alpha, m_shininess);
        }
        return ((m_diffuseReflectance / M_PI) + m_specularReflectance * specular) * Frame::cosTheta(bRec.wo);
    }

    /// Compute the density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        /* return zero if the measure is wrong, or when
        queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        float diffuseProb = Warp::squareToCosineHemispherePdf(bRec.wo);
        float alpha = bRec.wo.dot(reflect(bRec.wi));
        float specularProb = 0.0f;
        if(alpha > 0) {
            specularProb = (m_shininess + 1)/(2*M_PI) * std::pow(alpha, m_shininess);
        }
        return m_specularSamplingWeight * specularProb + (1-m_specularSamplingWeight) * diffuseProb;
    }

    /// Draw a a sample from the BRDF model
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
        if (Frame::cosTheta(bRec.wi) <= 0)
            return Color3f(0.0f);

        // adapt sample according to if it's diffuse or specular
        Point2f sample(_sample);
        bool choseSpecular = true;
        if (sample.x() <= m_specularSamplingWeight) {
            sample = Point2f(sample.x() / m_specularSamplingWeight, sample.y());
        } else {
            sample = Point2f((sample.x() - m_specularSamplingWeight) / (1-m_specularSamplingWeight), sample.y());
            choseSpecular = false;
        }

        // do the correct sampling depending on which part we sample (diffuse or specular)
        if (choseSpecular) {
            // sample from a Phong lobe centered around (0, 0, 1)
            float theta = acos(pow( (1 - sample.x()), 1.0f/(m_shininess+2) ));
            float phi = 2*M_PI*sample.y();
            Vector3f localDir = sphericalDirection(theta, phi);

            // rotate into the correct coordinate system
            Vector3f R = reflect(bRec.wi);
            bRec.wo = Frame(R).toWorld(localDir);

            if (Frame::cosTheta(bRec.wo) <= 0)
                return 0.0f;
        } else {
            // diffuse sampling
            bRec.wo = Warp::squareToCosineHemisphere(sample);
        }

        // Relative index of refraction: no change
        bRec.eta = 1.0f;
        float _pdf = pdf(bRec);
        if (_pdf == 0.0f) {
            return 0.0f;
        }
        return eval(bRec) / _pdf;
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
    const float m_shininess;
    Color3f m_diffuseReflectance;
    Color3f m_specularReflectance;
    float m_specularSamplingWeight;
};

NORI_REGISTER_CLASS(Phong, "phong");
NORI_NAMESPACE_END
