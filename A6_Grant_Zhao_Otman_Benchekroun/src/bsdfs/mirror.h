/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include "core/core.h"

TR_NAMESPACE_BEGIN

/**
 * Perfectly specular reflectance model
 */
struct MirrorBSDF : BSDF {

    const float DeltaEpsilon = 0.00001f;

    MirrorBSDF(const WorldData& scene, const Config& config, const size_t& matID) : BSDF(scene, config, matID) {
        components.push_back(EDeltaReflection);

        combinedType = 0;
        for (size_t i = 0; i < components.size(); ++i)
            combinedType |= components[i];
    }

    v3f eval(const SurfaceInteraction& i) const override {
        if (Frame::cosTheta(i.wi) <= 0 || Frame::cosTheta(i.wo) <= 0)
            return v3f(0.0, 0.0, 0.0);

        if (Frame::cosTheta(i.wi) * Frame::cosTheta(i.wo) >= 0) {   //are wi wo on same side


            if (std::abs(dot(reflect(i.wo), i.wi) - 1) > DeltaEpsilon)
                return v3f(0.0);
        }

        return v3f(1.0, 1.0, 1.0);
    }

    float pdf(const SurfaceInteraction& i) const override {
        if (Frame::cosTheta(i.wi) <= 0 || Frame::cosTheta(i.wo) <= 0)
            return 0.0f;

        return 1.0;
    }

    v3f sample(SurfaceInteraction& i, const v2f& sample, float* pdf) const override {
        if (Frame::cosTheta(i.wo) <= 0)
            return v3f(0.0f);

        i.sampledComponent = 0;
        i.sampledType = EDeltaReflection;

        i.wi = reflect(i.wo);
        if (pdf) *pdf = Warp::squareToCosineHemispherePdf(i.wi);

        return v3f(1.0, 1.0, 1.0);
    }

    inline v3f reflect(const v3f& d) const {
        return v3f(-d.x, -d.y, d.z);
    }

    std::string toString() const override { return "Mirror"; }
};

TR_NAMESPACE_END