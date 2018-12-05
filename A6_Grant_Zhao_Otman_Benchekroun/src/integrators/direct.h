/*
	This file is part of TinyRender, an educative rendering system.
	Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
	Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Direct illumination integrator with MIS
 */
	struct DirectIntegrator : Integrator {
	explicit DirectIntegrator(const Scene& scene) : Integrator(scene) {
		m_emitterSamples = scene.config.integratorSettings.di.emitterSamples;
		m_bsdfSamples = scene.config.integratorSettings.di.bsdfSamples;
		m_samplingStrategy = scene.config.integratorSettings.di.samplingStrategy;
	}

	static inline float balanceHeuristic(float nf, float fPdf, float ng, float gPdf) {
		float f = nf * fPdf, g = ng * gPdf;
		return f / (f + g);
	}

	void sampleSphereByCosineHemisphere(const p2f& sample,
		const v3f& n,
		const p3f& pShading,
		const v3f& emitterCenter,
		float emitterRadius,
		v3f& wiW,
		float& pdf) const {
		// TODO: Implement this

	}

	void sampleSphereByArea(const p2f& sample,
		const p3f& pShading,
		const v3f& emitterCenter,
		float emitterRadius,
		v3f& pos,
		v3f& ne,
		v3f& wiW,
		float& pdf) const {
		// TODO: Implement this
		ne = Warp::squareToUniformSphere(sample);//unit sphere
		v3f samplePoint = ne * emitterRadius;//sphere with the emitter radius
		pos = samplePoint + emitterCenter;//world
		wiW = glm::normalize(pos - pShading);//world
		pdf = INV_FOURPI / pow(emitterRadius, 2);


	}

	void sampleSphereBySolidAngle(const p2f& sample,
		const p3f& pShading,
		const v3f& emitterCenter,
		float emitterRadius,
		v3f& wiW,
		float& pdf) const {
		// TODO: Implement this

		float sinThetaMax2 = emitterRadius * emitterRadius / glm::distance2(emitterCenter, pShading);
		float cosThetaMax = glm::sqrt(max(0.f, 1 - sinThetaMax2));
		v3f ne = Warp::squareToUniformCone(sample, cosThetaMax);
		wiW = glm::mat4(glm::quat(v3f(0, 0, 1), glm::normalize(emitterCenter - pShading)))*v4f(ne, 1);
		wiW = glm::normalize(wiW);
		pdf = Warp::squareToUniformConePdf(cosThetaMax);
	}

	v3f renderArea(const Ray& ray, Sampler& sampler) const {
		v3f Lr(0.f);
		// TODO: Implement this
		int emitterSamples = m_emitterSamples;
		SurfaceInteraction i;
		SurfaceInteraction shadow;
		if (scene.bvh->intersect(ray, i)) {
			v3f emission = getEmission(i);
			if (emission != v3f(0.f)) {
				return emission;
			}
			for (int j = 0; j < emitterSamples; j++) {
				const p2f sample = sampler.next2D();
				const p3f pShading = i.p;
				float emPdf;
				size_t id = selectEmitter(sampler.next(), emPdf);
				const Emitter& em = getEmitterByID(id);
				const v3f emitterCenter = scene.getShapeCenter(em.shapeID);
				float emitterRadius = scene.getShapeRadius(em.shapeID);
				v3f ne, pos, wiW;
				float pdf;
				sampleSphereByArea(sample, pShading, emitterCenter, emitterRadius, pos, ne, wiW, pdf);
				i.wi = i.frameNs.toLocal(wiW);
				Ray shadowRay = Ray(i.p, wiW);
				float rayAngle = glm::dot(glm::normalize(emitterCenter - i.p), ne);
				if (rayAngle >= 0) {
					if (scene.bvh->intersect(shadowRay, shadow)) {

						v3f emission = getEmission(shadow);
						if (emission != v3f(0.f)) {

							float pdfOmega = glm::distance2(i.p, pos)*pdf / glm::abs(glm::dot(wiW, ne));
							Lr += getBSDF(i)->eval(i)*emission / pdfOmega / emPdf;
						}
					}

				}


			}


		}
		return Lr / emitterSamples;
	}

	v3f renderCosineHemisphere(const Ray& ray, Sampler& sampler) const {
		v3f Lr(0.f);
		// TODO: Implement this
		int emitterSamples = m_emitterSamples;
		SurfaceInteraction i;
		SurfaceInteraction shadow;
		if (scene.bvh->intersect(ray, i)) {
			v3f emission = getEmission(i);
			if (emission != v3f(0.f)) {
				return emission;
			}

			for (int j = 0; j < emitterSamples; j++) {

				v3f wi = Warp::squareToCosineHemisphere(sampler.next2D());
				float pdf = Warp::squareToCosineHemispherePdf(wi);
				i.wi = wi;
				wi = glm::normalize(i.frameNs.toWorld(wi));
				Ray shadowRay = Ray(i.p, wi);
				if (scene.bvh->intersect(shadowRay, shadow)) {//check if the ray from x was blocked
					v3f emission = getEmission(shadow);
					if (emission != v3f(0.f)) {
						v3f brdf = getBSDF(i)->eval(i);
						Lr += brdf * emission / pdf;
					}
				}

			}
		}
		return Lr / emitterSamples;
	}

	v3f renderBSDF(const Ray& ray, Sampler& sampler) const {
		v3f Lr(0.f);
		// TODO: Implement this
		int emitterSamples = m_emitterSamples;
		SurfaceInteraction i;
		SurfaceInteraction shadow;
		if (scene.bvh->intersect(ray, i)) {
			v3f emission = getEmission(i);
			if (emission != v3f(0.f)) {
				return emission;
			}

			for (int j = 0; j < emitterSamples; j++) {
				float pdf;
				v3f brdf = getBSDF(i)->sample(i, sampler.next2D(), &pdf);

				v3f wi = i.wi;
				wi = glm::normalize(i.frameNs.toWorld(wi));

				Ray shadowRay = Ray(i.p, wi);
				if (scene.bvh->intersect(shadowRay, shadow)) {//check if the ray from x was blocked
					v3f emission = getEmission(shadow);
					if (emission != v3f(0.f)) {
						//v3f brdf = getBSDF(i)->eval(i);
						Lr += brdf * emission;
					}
				}

			}
		}
		return Lr / emitterSamples;

	}

	v3f renderSolidAngle(const Ray& ray, Sampler& sampler) const {
		v3f Lr(0.f);
		// TODO: Implement this
		int emitterSamples = m_emitterSamples;
		SurfaceInteraction i;
		SurfaceInteraction shadow;
		if (scene.bvh->intersect(ray, i)) {
			v3f emission = getEmission(i);
			if (emission != v3f(0.f)) {
				return emission;
			}
			for (int j = 0; j < emitterSamples; j++) {
				const p2f sample = sampler.next2D();
				const p3f pShading = i.p;
				float emPdf;
				size_t id = selectEmitter(sampler.next(), emPdf);
				const Emitter& em = getEmitterByID(id);
				const v3f emitterCenter = scene.getShapeCenter(em.shapeID);
				float emitterRadius = scene.getShapeRadius(em.shapeID);
				float pdf;
				v3f wiW;
				sampleSphereBySolidAngle(sample, pShading, emitterCenter, emitterRadius, wiW, pdf);
				i.wi = i.frameNs.toLocal(wiW);
				Ray shadowRay = Ray(i.p, wiW);
				if (scene.bvh->intersect(shadowRay, shadow)) {
					v3f emission = getEmission(shadow);
					if (emission != v3f(0.f)) {
						Lr += emission * getBSDF(i)->eval(i) / emPdf / pdf;
					}
				}
			}
		}
		return Lr / emitterSamples;
	}

	v3f renderMIS(const Ray& ray, Sampler& sampler) const {
		v3f Lr(0.f);
		v3f LrEmitter(0.f);
		v3f LrBsdf(0.f);
		// TODO: Implement this
		int emitterSamples = m_emitterSamples;
		int bsdfSamples = m_bsdfSamples;
		SurfaceInteraction i;
		SurfaceInteraction shadow;
		if (scene.bvh->intersect(ray, i)) {
			v3f emission = getEmission(i);
			if (emission != v3f(0.f)) {
				return emission;
			}


			for (int j = 0; j < emitterSamples; j++) {
				const p2f sample = sampler.next2D();
				const p3f pShading = i.p;
				float emPdf;
				size_t id = selectEmitter(sampler.next(), emPdf);
				const Emitter& em = getEmitterByID(id);
				const v3f emitterCenter = scene.getShapeCenter(em.shapeID);
				float emitterRadius = scene.getShapeRadius(em.shapeID);
				float pdf;
				v3f wiW;
				sampleSphereBySolidAngle(sample, pShading, emitterCenter, emitterRadius, wiW, pdf);
				i.wi = i.frameNs.toLocal(wiW);
				Ray shadowRay = Ray(i.p, wiW);
				if (scene.bvh->intersect(shadowRay, shadow)) {
					v3f emission = getEmission(shadow);
					if (emission != v3f(0.f)) {
						float pdfBsdf = getBSDF(i)->pdf(i);
						float pdfSA = pdf * emPdf;
						float we = balanceHeuristic(emitterSamples, pdfSA, bsdfSamples, pdfBsdf);
						LrEmitter += emission * getBSDF(i)->eval(i)*we* Frame::cosTheta(i.wi) / pdfSA;
					}
				}
			}


			for (int j = 0; j < bsdfSamples; j++) {
				float pdf;
				const p2f sample = sampler.next2D();
				v3f brdf = getBSDF(i)->sample(i, sample, &pdf)*(Frame::cosTheta(i.wi));

				v3f wi = i.wi;
				wi = glm::normalize(i.frameNs.toWorld(wi));

				Ray shadowRay = Ray(i.p, wi);
				if (scene.bvh->intersect(shadowRay, shadow)) {
					v3f emission = getEmission(shadow);
					if (emission != v3f(0.f)) {
						const p3f pShading = i.p;
						float emPdf = 1.0 / scene.emitters.size();
						size_t id = selectEmitter(sampler.next(), emPdf);
						const Emitter& em = getEmitterByID(getEmitterIDByShapeID(shadow.shapeID));
						const v3f emitterCenter = scene.getShapeCenter(em.shapeID);
						float emitterRadius = scene.getShapeRadius(em.shapeID);
						float sinThetaMax2 = emitterRadius * emitterRadius / glm::distance2(emitterCenter, pShading);
						float cosThetaMax = glm::sqrt(max(0.f, 1 - sinThetaMax2));
						float pdfSA = Warp::squareToUniformConePdf(cosThetaMax)*emPdf;
						float wb = balanceHeuristic(bsdfSamples, pdf, emitterSamples, pdfSA);
						LrBsdf += brdf * emission*wb;
					}
				}

			}
			if (emitterSamples <= 0) {
				Lr = LrBsdf / bsdfSamples;
			}
			else if (bsdfSamples <= 0) {
				Lr = LrEmitter / emitterSamples;
			}
			else {
				Lr = LrEmitter / emitterSamples + LrBsdf / bsdfSamples;
			}

		}


		return Lr;
	}

	v3f render(const Ray& ray, Sampler& sampler) const override {
		if (m_samplingStrategy == "mis")
			return this->renderMIS(ray, sampler);
		else if (m_samplingStrategy == "area")
			return this->renderArea(ray, sampler);
		else if (m_samplingStrategy == "solidAngle")
			return this->renderSolidAngle(ray, sampler);
		else if (m_samplingStrategy == "cosineHemisphere")
			return this->renderCosineHemisphere(ray, sampler);
		else if (m_samplingStrategy == "bsdf")
			return this->renderBSDF(ray, sampler);
		std::cout << "Error: wrong strategy" << std::endl;
		exit(EXIT_FAILURE);
	}

	size_t m_emitterSamples;     // Number of emitter samples
	size_t m_bsdfSamples;        // Number of BSDF samples
	string m_samplingStrategy;   // Sampling strategy to use
};

TR_NAMESPACE_END