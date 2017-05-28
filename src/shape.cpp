
#include <nori/shape.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

Shape::Shape() {}

Shape::~Shape()
{
	delete m_bsdf;
	delete m_emitter;
}

void Shape::activate() {
	if (!m_bsdf) {
		/* If no material was assigned, instantiate a diffuse BRDF */
		m_bsdf = static_cast<BSDF *>(
			NoriObjectFactory::createInstance("diffuse", PropertyList()));
	}
}

void Shape::addChild(NoriObject *obj) {
	switch (obj->getClassType()) {
	case EBSDF:
		if (m_bsdf)
			throw NoriException(
				"Shape: tried to register multiple BSDF instances!");
		m_bsdf = static_cast<BSDF *>(obj);
		break;

	case EEmitter: {
		Emitter *emitter = static_cast<Emitter *>(obj);
		if (m_emitter)
			throw NoriException(
				"Shape: tried to register multiple Emitter instances!");
		m_emitter = emitter;
		m_emitter->setShape(this);
	}
		break;

	default:
		throw NoriException("Shape::addChild(<%s>) is not supported!",
			classTypeName(obj->getClassType()));
	}
}

std::string Shape::toString() const {
	return tfm::format(
		"Shape[\n"
		"  bsdf = %s,\n"
		"  emitter = %s\n"
		"]",
		m_bsdf ? indent(m_bsdf->toString()) : std::string("null"),
		m_emitter ? indent(m_emitter->toString()) : std::string("null")
	);
}

std::string Intersection::toString() const {
	if (!shape)
		return "Intersection[invalid]";

	return tfm::format(
		"Intersection[\n"
		"  p = %s,\n"
		"  t = %f,\n"
		"  uv = %s,\n"
		"  shFrame = %s,\n"
		"  geoFrame = %s,\n"
		"  shape = %s\n"
		"]",
		p.toString(),
		t,
		uv.toString(),
		indent(shFrame.toString()),
		indent(geoFrame.toString()),
		shape ? shape->toString() : std::string("null")
	);
}

NORI_NAMESPACE_END