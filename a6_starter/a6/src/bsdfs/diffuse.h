/*
    This file is part of TinyRender, an educative rendering system.
    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include "core/core.h"

TR_NAMESPACE_BEGIN

/**
 * Perfectly diffuse, Lambertian reflectance model
 */
struct DiffuseBSDF : BSDF {
    std::unique_ptr<Texture < v3f>> albedo;

    DiffuseBSDF(const WorldData& scene, const Config& config, const size_t& matID) : BSDF(scene, config, matID) {
        const tinyobj::material_t& mat = scene.materials[matID];

        if (mat.diffuse_texname.empty())
            albedo = std::unique_ptr<Texture<v3f>>(new ConstantTexture3f(glm::make_vec3(mat.diffuse)));
        else
            albedo = std::unique_ptr<Texture<v3f>>(new BitmapTexture3f(config, mat.diffuse_texname));

        components.push_back(EDiffuseReflection);

        combinedType = 0;
        for (size_t i = 0; i < components.size(); ++i)
            combinedType |= components[i];
    }

    v3f eval(const SurfaceInteraction& i) const override {
		v3f val(0.f);
		// TODO: Add previous assignment code (if needed)
		v3f rho = albedo->eval(worldData, i);
		v3f wo = i.wo;
		v3f wi = i.wi;
		float z_in = Frame::cosTheta(wi);
		float z_out = Frame::cosTheta(wo);

		if (z_in > 0 && z_out > 0) //check correct direction				
			val = rho / M_PI;
		return val;
    }

    float pdf(const SurfaceInteraction& i) const override {
		// TODO: Implement this
		float pdf = Warp::squareToCosineHemispherePdf(i.wi);
		return pdf;
    }

    v3f sample(SurfaceInteraction& i, const v2f& sample, float* pdf_p) const override {
		v3f val(0.f);
		i.wi = Warp::squareToCosineHemisphere(sample);

		float pdf_val = pdf(i);
		*pdf_p = pdf_val;

		v3f brdf_factor = eval(i);

		if (pdf_val > 0.f)
			val = brdf_factor / pdf_val;
		else
			val = v3f(0);
		return val;
    }

    std::string toString() const override { return "Diffuse"; }
};

TR_NAMESPACE_END