/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include "core/core.h"

TR_NAMESPACE_BEGIN

/**
 * Modified Phong reflectance model
 */
struct PhongBSDF : BSDF {

    std::unique_ptr<Texture < v3f>> specularReflectance;
    std::unique_ptr<Texture < v3f>> diffuseReflectance;
    std::unique_ptr<Texture < float>> exponent;
    float specularSamplingWeight;
    float scale;

    PhongBSDF(const WorldData& scene, const Config& config, const size_t& matID) : BSDF(scene, config, matID) {
        const tinyobj::material_t& mat = scene.materials[matID];

        if (mat.specular_texname.empty())
            specularReflectance = std::unique_ptr<Texture<v3f>>(new ConstantTexture3f(glm::make_vec3(mat.specular)));
        else
            specularReflectance = std::unique_ptr<Texture<v3f>>(new BitmapTexture3f(config, mat.specular_texname));

        if (mat.diffuse_texname.empty())
            diffuseReflectance = std::unique_ptr<Texture<v3f>>(new ConstantTexture3f(glm::make_vec3(mat.diffuse)));
        else
            diffuseReflectance = std::unique_ptr<Texture<v3f>>(new BitmapTexture3f(config, mat.diffuse_texname));

        exponent = std::unique_ptr<Texture<float>>(new ConstantTexture1f(mat.shininess));

        //get scale value to ensure energy conservation
        v3f maxValue = specularReflectance->getMax() + diffuseReflectance->getMax();
        float actualMax = max(max(maxValue.x, maxValue.y), maxValue.z);
        scale = actualMax > 1.0f ? 0.99f * (1.0f / actualMax) : 1.0f;

        float dAvg = getLuminance(diffuseReflectance->getAverage() * scale);
        float sAvg = getLuminance(specularReflectance->getAverage() * scale);
        specularSamplingWeight = sAvg / (dAvg + sAvg);

        components.push_back(EGlossyReflection);
        components.push_back(EDiffuseReflection);

        combinedType = 0;
        for (unsigned int component : components)
            combinedType |= component;
    }

    inline v3f reflect(const v3f& d) const {
        return v3f(-d.x, -d.y, d.z);
    }

	v3f eval(const SurfaceInteraction& i) const override {
		v3f val(0.f);
		// TODO: Implement 
		v3f wo = i.wo;
		v3f wi = i.wi;
		if (Frame::cosTheta(wo) < 0.f || Frame::cosTheta(wi) < 0.f) {
			return val;
		}
		float cosWi = Frame::cosTheta(wi);
		float cosWo = Frame::cosTheta(wo);
		v3f diffuseR = diffuseReflectance->eval(worldData, i);
		v3f specularF = specularReflectance->eval(worldData, i);
		float n = exponent->eval(worldData, i);
		//a dot b = |a||b|cos(theta)
		float cosa = fmax(glm::dot(glm::normalize(reflect(wi)), glm::normalize(wo)), 0.f);

		val = diffuseR * INV_PI + specularF * (n + 2)*INV_TWOPI*pow(cosa, n);

		//val = specularF * (n + 2)*INV_TWOPI*pow(cosa, n);
		return val * cosWi*scale;
	}

	float pdf(const SurfaceInteraction& i) const override {
		float pdf = 0.f;
		v3f reflection = reflect(i.wo);
		v3f wi = glm::toMat4(glm::quat(reflection, v3f(0, 0, 1))) * v4f(i.wi, 1);
		pdf = Warp::squareToPhongLobePdf(wi, exponent->eval(worldData, i));
		return pdf;
	}

	v3f sample(SurfaceInteraction& i, const v2f& _sample, float* pdf1) const override {
		v3f val(0.f);
		// TODO: Implement this
		v3f wi = Warp::squareToPhongLobe(_sample, exponent->eval(worldData, i));
		v3f reflection = reflect(i.wo);
		wi = glm::toMat4(glm::quat(v3f(0, 0, 1), reflection)) * v4f(wi, 1);
		i.wi = wi;
		float Pdf = pdf(i);
		*pdf1 = Pdf;
		if (Pdf == 0.f) {
			return val;
		}
		v3f brdf = eval(i);
		val = brdf / Pdf;
		return val;
	}


    std::string toString() const override { return "Phong"; }
};

TR_NAMESPACE_END