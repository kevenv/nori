/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/accel.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

void Accel::addShape(Shape *shape) {
    if (m_shape)
        throw NoriException("Accel: only a single shape is supported!");
	m_bvh.addShape(shape);
    m_shape = m_bvh.getShape(0);
    m_bbox = m_bvh.getBoundingBox();
}

void Accel::build() {
	m_bvh.build();
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
	return m_bvh.rayIntersect(ray_, its, shadowRay);
}

NORI_NAMESPACE_END

