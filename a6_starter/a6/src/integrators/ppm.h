
/*
	This file is part of TinyRender, an educative rendering system.
	Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
	Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include <core/core.h>
#include <kdtree.h>
#include "direct.h"

TR_NAMESPACE_BEGIN

/* This is the photon
 * The power is not compressed so the
 * size is 28 bytes
*/
typedef struct Photon {
	v3f pos;                  // photon position
	v3f dir;      // incoming direction
	v3f power;   // photon power (uncompressed)
	v3f n;
} Photon;


/**
 * Photon mapping integrator
 */
struct PPMIntegrator : Integrator {

	std::vector<Photon> m_photonMap;

	typedef unsigned int PhotonMapIdx;
	typedef GenericKDTreeNode<v3f, PhotonMapIdx> PhotonKDTreeNode;
	PointKDTree<PhotonKDTreeNode> m_KDTree;

	//1st pass
	int m_photonCount;
	float m_photonRrProb;
	int m_photonRrDepth;

	//2nd pass
	int m_emittedPhotonCount;
	bool m_usePhotonsForDirect;
	float m_radiusSearch;
	int m_nbPhotonsSearch;
	bool m_useFinalGather;
	int m_nbFinalGather;

	std::unique_ptr<DirectIntegrator> m_directIntegrator;

	explicit PPMIntegrator(const Scene& scene) : Integrator(scene)
	{
		m_directIntegrator = std::unique_ptr<DirectIntegrator>(new DirectIntegrator(scene));

		m_directIntegrator->m_emitterSamples = scene.config.integratorSettings.pm.emitterSamplesCount;
		m_directIntegrator->m_bsdfSamples = scene.config.integratorSettings.pm.emitterSamplesCount;
		m_directIntegrator->m_samplingStrategy = "mis";

		//1st pass
		m_photonCount = scene.config.integratorSettings.pm.photonCount;
		m_photonRrDepth = scene.config.integratorSettings.pm.photonRrDepth;
		m_photonRrProb = scene.config.integratorSettings.pm.photonRrProb;

		//2nd pass
		m_radiusSearch = scene.config.integratorSettings.pm.searchRadius;
		m_nbPhotonsSearch = scene.config.integratorSettings.pm.photonsSearchCount;
		m_useFinalGather = scene.config.integratorSettings.pm.useFinalGather;
		m_nbFinalGather = scene.config.integratorSettings.pm.finalGatherSamplesCount;
		m_usePhotonsForDirect = scene.config.integratorSettings.pm.usePhotonsForDirect;
	}

	bool init() override {
		Integrator::init();

		std::cout << "Start emitting photons. " << std::endl;
		generatePhotonMap();
		printf("first bounce hit:%d", firstBH);
		return true;
	}

	void generatePhotonMap() {
		// TODO: Implement this
		float totalLightArea = 0.f;
		Sampler sampler = Sampler(260563769);

		for (unsigned int i = 0; i < scene.emitters.size(); i++) {
			const Emitter& em = scene.emitters[i];
			totalLightArea += em.area;
		}
		/* loop through all lights
		for (int i = 0; i < scene.emitters.size(); i++) {
			const Emitter& em = scene.emitters[i];
			float lightArea = em.area;
			int emitterPhotonNum = (int)ceil(m_photonCount*lightArea / totalLightArea);
			//initial energy
			v3f energy = lightArea / double(emitterPhotonNum)*em.getPower();
			for (int j = 0; j <= emitterPhotonNum; j++) {
			}
		}
		*/
		//randomly select light
		for (int i = 0; i < m_photonCount; i++) {
			float emPdf;
			size_t id = selectEmitter(sampler.next(), emPdf);
			const Emitter& em = getEmitterByID(id);
			float lightArea = em.area;
			//initial photon from emitter
			v3f n, pos, dir;
			float pdfPos, pdfDir;
			sampleEmitterPosition(sampler, em, n, pos, pdfPos);
			sampleEmitterDirection(sampler, em, n, dir, pdfDir);
			int emitterPhotonNum = (int)ceil(m_photonCount*lightArea / totalLightArea);
			v3f energy = em.getPower() / emitterPhotonNum;
			tracePhoton(sampler, pos, dir, energy, 0);
		}
		m_KDTree.build();
	}
	int firstBH = 0;
	void tracePhoton(Sampler& sampler, const v3f& pos, const v3f& dir, const v3f& energy, int bounces) {

		float rrProb = 1.f;
		if (bounces >= m_photonRrDepth) {
			if (sampler.next() >= m_photonRrProb) {
				return;
			}
			else {
				//return;
				rrProb = m_photonRrProb;
			}
		}
		//shot a ray in the direction
		Ray ray = Ray(pos, dir);
		SurfaceInteraction hit;
		if (scene.bvh->intersect(ray, hit)) {
			if (bounces == 0) {
				firstBH += 1;
			}
			v3f emission = getEmission(hit);
			//if not hitting emitter
			if (emission == v3f(0.f)) {//also determine whether surface is diffuse, but how?
				v3f wiW = glm::normalize(pos - hit.p);//light pos - hit pos
				Photon p;
				v3f normal = hit.frameNs.n;
				p.pos = hit.p;
				p.dir = dir;
				p.power = energy;
				p.n = normal;
				//add the photon to the back of the list
				//id = size of the list 
				m_photonMap.push_back(p);
				PhotonMapIdx curr = m_photonMap.size();
				const PhotonKDTreeNode currNode = PhotonKDTreeNode(pos, curr - 1);
				m_KDTree.push_back(currNode);
				//recurse photon mapping
				hit.wi = hit.frameNs.toLocal(wiW);
				float pdf;
				getBSDF(hit)->sample(hit, sampler.next2D(), &pdf);
				v3f wiW2 = glm::normalize(hit.frameNs.toWorld(hit.wi));
				hit.wo = hit.wi;
				hit.wi = hit.frameNs.toLocal(wiW);
				v3f bsdf = getBSDF(hit)->eval(hit);
				//TODO how to get the normal of the hit
				v3f power = energy * glm::abs(glm::dot(normal, wiW2)) *bsdf / rrProb / pdf;
				tracePhoton(sampler, hit.p, wiW2, power, bounces + 1);
				//select direction
			}


			//create a new photon
		}
	}
	v3f render(const Ray& ray, Sampler& sampler) const override {
		v3f throughput(0.f);
		//get Kd
		SurfaceInteraction hit;
		//shoot ray into scene.

		if (scene.bvh->intersect(ray, hit)) {
			v3f emission = getEmission(hit);
			if (emission != v3f(0.f)) {
				return emission;
			}
			
						
			const double num = m_nbPhotonsSearch;
			PointKDTree<PhotonKDTreeNode>::SearchResult results[501];

			float searchr =  m_radiusSearch* m_radiusSearch;
			m_KDTree.nnSearch(hit.p, searchr, m_nbPhotonsSearch, results);

			for (int i = 0; i < m_nbPhotonsSearch; i++) {//for all the nearest neighbors
				int index = results[i].index;
				if (index < m_photonMap.size()) {
					Photon p = m_photonMap[index]; //TODO  the photon normal check
					if (glm::abs(glm::dot(normalize(hit.frameNs.n), normalize(p.n))) -1 < 0.1) {
						hit.wo = -hit.wi;
						hit.wi = hit.frameNs.toLocal(p.dir);
						//hit.frameNs.n = p.n;
						v3f eval_bsdf = getBSDF(hit)->eval(hit);
						//eval_bsdf = v3f(1.f);
						//float radius2 = glm::distance2(p.pos, hit.p);
						float radius2 = m_radiusSearch * m_radiusSearch;
		
						throughput += eval_bsdf * p.power*INV_PI / (radius2);
					}

				}

				//throughput = v3f(1.f);
			}

		}


		// TODO: Implement this

		return throughput;
	}
};

TR_NAMESPACE_END

