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

        return true;
    }

    void generatePhotonMap() {
        // TODO: Implement this
		float totalLightArea = 0.f;
		Sampler sampler = Sampler(260563769);
		
		for (int i = 0; i < scene.emitters.size();i++) {
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
			v3f energy = em.getPower() /emitterPhotonNum/emPdf/pdfPos/pdfDir;
			tracePhoton(sampler, pos, dir, n, energy, 0);
		}
    }
	void tracePhoton(Sampler& sampler, const v3f& pos, const v3f& dir, const v3f& n, const v3f& energy, int bounces){
		float rrProb = 1.f;
		if (bounces >= m_photonRrDepth) {
			if (sampler.next() >= m_photonRrProb) {
				return;
			}
			else {
				rrProb = m_photonRrProb;
			}
		}
		//shot a ray in the direction
		Ray ray = Ray(pos, dir);
		SurfaceInteraction hit;
		if (scene.bvh->intersect(ray, hit)) {
			v3f emission = getEmission(hit);
			//if not hitting emitter
			if (emission == v3f(0.f)) {//also determine whether surface is diffuse, but how?
				v3f wiW = glm::normalize(pos - hit.p);
				Photon p;
				p.pos = hit.p;
				p.dir = wiW;
				p.power = energy;
				p.n = n;
				//add the photon to the back of the list
				//id = size of the list 
				m_photonMap.push_back(p);
				PhotonMapIdx curr = m_photonMap.size();
				const PhotonKDTreeNode currNode = PhotonKDTreeNode(pos,curr-1);
				m_KDTree.push_back(currNode);
				//recurse photon mapping
				float pdf;
				v3f BSDF = getBSDF(hit)->sample(hit, sampler.next2D(), &pdf);
				v3f normal = v3f(0.f);
				//TODO how to get the normal of the hit
				v3f power = energy * glm::abs(glm::dot(normal, wiW))*BSDF / pdf;
				tracePhoton(sampler, hit.p, wiW, power, normal,bounces + 1);
				
				//select direction

			}
			
			
			//create a new photon
		}
	}
    v3f render(const Ray& ray, Sampler& sampler) const override {
        v3f throughput(0.f);
        // TODO: Implement this

        return throughput;
    }
};

TR_NAMESPACE_END