
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
	float m_radiusSqr;
	int randomSeed = 265;
	std::unique_ptr<DirectIntegrator> m_directIntegrator;

	explicit PPMIntegrator(const Scene& scene) : Integrator(scene)
	{
		m_directIntegrator = std::unique_ptr<DirectIntegrator>(new DirectIntegrator(scene));

		m_directIntegrator->m_emitterSamples = scene.config.integratorSettings.pm.emitterSamplesCount;
		m_directIntegrator->m_bsdfSamples = scene.config.integratorSettings.pm.emitterSamplesCount;
		m_directIntegrator->m_samplingStrategy = "mis";
		// = 0;
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
	float radiusCalculator() {
		static int iter = 0;
		static float radiusSqr = m_radiusSearch*m_radiusSearch;
		
		float factor;
		if (iter == 0) {
			factor = 1.f;
		}
		else {
			factor = (((2.f / 3.f) + (float)(iter-1.f)) / (1.f + (float)(iter-1.f)));
		}
		
		radiusSqr = radiusSqr * factor;
		iter++;
		return radiusSqr;
	}
	float randomSeedCalculator() {
		static int iter = 1;
		static int raseed = 51;
		raseed+=iter;
		iter++;
		return raseed;
	}
	bool init() override {
		Integrator::init();
		//m_radiusSearch = iter / 10;
		//iter++;
		m_radiusSqr = radiusCalculator();
		randomSeed = randomSeedCalculator();
		std::vector<Photon>().swap(m_photonMap);
		m_photonMap.clear();
		m_photonMap.shrink_to_fit();
		m_KDTree.clear();
		

		printf("\nradius square: %f , new random seed: %d\n",m_radiusSqr,randomSeed);
		std::cout << "Start emitting photons. " << std::endl;
		

		generatePhotonMap();
		return true;
	}
	bool haveEnoughPhotons = false;

	void generatePhotonMap() {
		// TODO: Implement this
		float totalLightArea = 0.f;
		Sampler sampler = Sampler(randomSeed);
		for (unsigned int i = 0; i < scene.emitters.size(); i++) {
			const Emitter& em = scene.emitters[i];
			totalLightArea += em.area;
		}
		//randomly select light
		int emittedPhotons = 0;
		try {
			m_photonMap.reserve(m_photonCount);
			m_KDTree.reserve(m_photonCount);
		}
		catch (const std::bad_alloc& e) {
			std::cout << "Allocation failed: " << e.what() << '\n';
		}
		m_photonMap.reserve(m_photonCount);
		m_KDTree.reserve(m_photonCount);
		while (!haveEnoughPhotons) {
			float emPdf;
			size_t id = selectEmitter(sampler.next(), emPdf);
			const Emitter& em = getEmitterByID(id);
			float lightArea = em.area;
			//initial photon from emitter
			v3f n, pos, dir;
			float pdfPos, pdfDir;
			sampleEmitterPosition(sampler, em, n, pos, pdfPos);
			sampleEmitterDirection(sampler, em, n, dir, pdfDir);
			emittedPhotons++;
			int emitterPhotonNum = (int)ceil(m_emittedPhotonCount*lightArea / totalLightArea);
			v3f energy = em.getPower();

			tracePhoton(sampler, pos, dir, energy, 0);

		}
		cout << "total number of emitted photons: " << emittedPhotons << endl;
		for (int i = 0; i < m_photonMap.size(); i++) {
			m_photonMap[i].power = m_photonMap[i].power / emittedPhotons;
		}
	
		m_KDTree.build();

	}

	bool survivedRussianRoullette(v3f newEnergy, v3f oldEnergy) {
		float p = fmin(1.f, newEnergy.length() / oldEnergy.length());
		Sampler rand = Sampler(260663493);
		double num = rand.next();
		if (num > p) {
			return false;
		}
		else {
			newEnergy = newEnergy / p;
			return true;
		}

	}

	void tracePhoton(Sampler& sampler, const v3f& pos, const v3f& dir, const v3f& energy, int bounces) {

		float rrProb = 1.f;
		if (bounces >= m_photonRrDepth) {
			if (sampler.next() >= m_photonRrProb) {
				return;
			}
			else {
				rrProb = m_photonRrProb;
			}
		}
		/*if (bounces > 0) {
			return;
		}*/
		if (m_photonMap.size() >= m_photonCount) {
			cout << "Photon Map is full" << endl;
			haveEnoughPhotons = true;
			return;
		}

		//shot a ray in the direction
		Ray ray = Ray(pos, dir);
		SurfaceInteraction hit;
		if (scene.bvh->intersect(ray, hit)) {
			v3f emission = getEmission(hit);
			//if not hitting emitter
			if (emission == v3f(0.f)) {//also determine whether surface is diffuse, but how?
				//v3f wiW = glm::normalize(pos - hit.p);//light pos - hit pos
				Photon p;
				v3f normal = (hit.frameNs.n);
				float epsilon = 0;
				p.pos = hit.p; //+epsilon * hit.frameNs.n;
				p.n = normal;
				p.dir = dir;
				p.power = energy;
				float pdf;
				const BSDF* bsdf = getBSDF(hit);
				v3f brdf = bsdf->sample(hit, sampler.next2D(), &pdf);
				v3f wiW2 = glm::normalize(hit.frameNs.toWorld(hit.wi));
				v3f power = energy * brdf* glm::abs(glm::dot(normal, normalize(wiW2))) / rrProb;
				PhotonMapIdx curr = m_photonMap.size();
				const PhotonKDTreeNode currNode = PhotonKDTreeNode(p.pos, curr);
				if (m_usePhotonsForDirect||m_useFinalGather) {
					m_photonMap.push_back(p);
					m_KDTree.push_back(currNode);
				}
				else if(bounces > 0) {
					m_photonMap.push_back(p);
					m_KDTree.push_back(currNode);
				}
				
				tracePhoton(sampler, hit.p, wiW2, power, bounces + 1);
			}

		}
		//create a new photon

	}




	v3f render(const Ray& ray, Sampler& sampler) const override {
		v3f throughput(0.f);
		v3f direct = v3f(0.f);
		SurfaceInteraction hit;

		if (scene.bvh->intersect(ray, hit)) {

			v3f emission = getEmission(hit);
			if (emission != v3f(0.f)) {
				return emission;
			}
			//see if ray intersects geometry
			//if it does, go a closest neighbor search
			
			float radius = m_radiusSqr;
			if (m_useFinalGather) {
				direct = m_directIntegrator->render(ray, sampler);
				
				for (int j = 0; j < m_nbFinalGather; j++) {
					float bsdfSamplePdf;
					SurfaceInteraction indirectHit;
					v3f brdfDirect =  getBSDF(hit)->sample(hit, sampler.next2D(), &bsdfSamplePdf);
					v3f wiW = glm::normalize(hit.frameNs.toWorld(hit.wi));
					
					Ray indRay = Ray(hit.p, wiW);
					if (scene.bvh->intersect(indRay, indirectHit)) {
						const double num_photons = m_nbPhotonsSearch;
						PointKDTree<PhotonKDTreeNode>::SearchResult results[501];
						size_t n = m_KDTree.nnSearch(indirectHit.p, radius, m_nbPhotonsSearch, results);
						
						for (int i = 0; i < n; i++) {

							PointKDTree<PhotonKDTreeNode>::SearchResult result = results[i];
							//find the photons in the photon map;
							Photon photon = m_photonMap[m_KDTree[result.index].data];
							//add the contribution of the photon to the throughput
							//get the direction from hit point to photon
							indirectHit.wi = glm::normalize(indirectHit.frameNs.toLocal(-photon.dir));
							indirectHit.wo = glm::normalize(indirectHit.frameNs.toLocal(-wiW));
							
							const BSDF* brdf = getBSDF(indirectHit);

							// brdf_factor
							v3f brdf_factor = brdf->eval(indirectHit);

							double dot = glm::dot(photon.dir, (indirectHit.frameNs.n));
							if (dot < 0 && n > 2) { //approximately in the same direction. Photon is on the same surface
								throughput += brdfDirect * brdf_factor * photon.power*INV_PI / radius / m_nbFinalGather * Frame::cosTheta(hit.wi);

							}
						}
					}
				}
			}
			else {
				const double num_photons = m_nbPhotonsSearch;
				PointKDTree<PhotonKDTreeNode>::SearchResult results[501];
				size_t n = m_KDTree.nnSearch(hit.p, radius, m_nbPhotonsSearch, results);
				if (!m_usePhotonsForDirect) {
					direct = m_directIntegrator->render(ray, sampler);
				}
				for (int i = 0; i < n; i++) {

					PointKDTree<PhotonKDTreeNode>::SearchResult result = results[i];
					//find the photons in the photon map;

					Photon photon = m_photonMap[m_KDTree[result.index].data];
					//add the contribution of the photon to the throughput
					//get the direction from hit point to photon
					hit.wi = glm::normalize(hit.frameNs.toLocal(-photon.dir));
					hit.wo = glm::normalize(hit.frameNs.toLocal(-ray.d));
					//hit.p = photon.pos;
					const BSDF* brdf = getBSDF(hit);

					// brdf_factor
					v3f brdf_factor = brdf->eval(hit);

					double dot = glm::dot(photon.dir, (hit.frameNs.n));
					if (dot < 0 && n > 2) { //approximately in the same direction. Photon is on the same surface
						throughput += brdf_factor * photon.power*INV_PI / radius;

					}
				}
			}
			

		}	
		
		
		
		
		return throughput +direct;
	}


	TR_NAMESPACE_END;
}

