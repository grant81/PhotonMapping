/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#include <core/core.h>
#include <core/accel.h>
#include <core/renderer.h>
#include <GL/glew.h>

#ifdef __APPLE__
#include "SDL.h"
#include <OpenGL/gl.h>
#else
#ifdef _WIN32
#include <GL/gl.h>
#include "SDL.h"
#else
#include <GL/gl.h>
#include "SDL2/SDL.h"
#endif
#endif


#include <bsdfs/diffuse.h>

#include <integrators/normal.h>
#include <renderpasses/normal.h>

#include <integrators/simple.h>
#include <renderpasses/simple.h>
#include <bsdfs/phong.h>

#include <integrators/ao.h>
#include <integrators/ro.h>
#include <renderpasses/ssao.h>

#include <integrators/direct.h>

#include <integrators/path.h>
#include <renderpasses/gi.h>
#include <bsdfs/mixture.h>

#include <integrators/ppm.h>
#include <bsdfs/mirror.h>

TR_NAMESPACE_BEGIN

Renderer::Renderer(const Config& config) : scene(config) { }
bool Renderer::init(const bool isRealTime, bool nogui) {
    realTime = isRealTime;
    this->nogui = nogui;
    realTimeCameraFree = false;

    if (!scene.load(isRealTime)) return false;

    if (realTime) {
        if (scene.config.renderpass == ENormalRenderPass) {
            renderpass = std::unique_ptr<NormalPass>(new NormalPass(scene));
        }
        else if (scene.config.renderpass == EDirectRenderPass) {
            renderpass = std::unique_ptr<SimplePass>(new SimplePass(scene));
        }
        else if (scene.config.renderpass == ESSAORenderPass) {
            renderpass = std::unique_ptr<SSAOPass>(new SSAOPass(scene));
        }
        else if (scene.config.renderpass == EGIRenderPass) {
            renderpass = std::unique_ptr<GIPass>(new GIPass(scene));
        }
        else {
            throw std::runtime_error("Invalid renderpass type");
        }

        bool succ = renderpass.get()->initOpenGL(scene.config.width, scene.config.height);
        if (!succ) return false;

        return renderpass->init(scene.config);
    } else {
        if (scene.config.integrator == ENormalIntegrator) {
            integrator = std::unique_ptr<NormalIntegrator>(new NormalIntegrator(scene));
        }
        else if (scene.config.integrator == EAOIntegrator) {
            integrator = std::unique_ptr<AOIntegrator>(new AOIntegrator(scene));
        } else if (scene.config.integrator == EROIntegrator) {
            integrator = std::unique_ptr<ROIntegrator>(new ROIntegrator(scene));
        }
        else if (scene.config.integrator == ESimpleIntegrator) {
            integrator = std::unique_ptr<SimpleIntegrator>(new SimpleIntegrator(scene));
        }
        else if (scene.config.integrator == EDirectIntegrator) {
            integrator = std::unique_ptr<DirectIntegrator>(new DirectIntegrator(scene));
        }
        else if (scene.config.integrator == EPathTracerIntegrator) {
            integrator = std::unique_ptr<PathTracerIntegrator>(new PathTracerIntegrator(scene));
        }
        else if (scene.config.integrator == EPhotonMapperIntegrator) {
            integrator = std::unique_ptr<PPMIntegrator>(new PPMIntegrator(scene));
        }
        else {
            throw std::runtime_error("Invalid integrator type");
        }
		
		return integrator->init();
    }
}

void Renderer::render() {
	int progressive = 0;
	int progressivePass = 20;
	if (realTime) {
		/**
		 * Your real-time rendering loop solution from A1 here.
		 */
		SDL_Event e;
		
		for (;;) {
			SDL_PollEvent(&e);
			if (e.type == SDL_QUIT) {
				SDL_Log("Program quit after %i ticks", e.quit.timestamp);
				break;
			}
			renderpass->updateCamera(e);
			renderpass->render();
			SDL_GL_SwapWindow(renderpass->window);
		}
	}
	else {
		
		if (progressive == 1) {
			int width = scene.config.width;
			int height = scene.config.height;
			int sampleNum = scene.config.spp;
			integrator->rgb->clear();
			Sampler sampler = Sampler(260563769);
			v3f o = scene.config.camera.o;
			//1
			float cameraPerspective = scene.config.camera.fov;
			glm::mat4 inverseView = glm::lookAt(o, scene.config.camera.at, scene.config.camera.up);
			float aspectRatio = (float)width / (float)height;
			float scaling = tan((scene.config.camera.fov*deg2rad) / 2.f);
			std::unique_ptr<v3f[]> dataBuffer = std::unique_ptr<v3f[]>(new v3f[width * height]);
			
			for (int i = 0; i < progressivePass;i++){
				//2
				
				if (i != 0) {
					integrator = std::unique_ptr<PPMIntegrator>(new PPMIntegrator(scene));
					integrator->init();
				}
				
				for (int x = 0; x < width; ++x) {
					for (int y = 0; y < height; ++y) {
						v3f colorSum = v3f(0.f, 0.f, 0.f);
						if (sampleNum == 1) {
							float xNDC = (x + 0.5) / width;
							float yNDC = (y + 0.5) / height;
							float xScreen = 2 * xNDC - 1;
							float yScreen = 1 - 2 * yNDC;
							float xCamera = xScreen * aspectRatio*scaling;
							float yCamera = yScreen * scaling;
							v4f P = v4f(xCamera, yCamera, -1.f, 1.f);
							v4f pWorld = P * inverseView;
							//v4f eyeTransformed = v4f(o[0],o[1],o[2],1.f)*inverseView;
							Ray ray = Ray(o, v3f(pWorld));
							
							v3f color = integrator->render(ray, sampler);
							dataBuffer[y*width + x] += color/progressivePass;
						}
						else {
							for (int z = 0; z < sampleNum; z++) {
								float xNDC = (x + sampler.next()) / width;
								float yNDC = (y + sampler.next()) / height;
								float xScreen = 2 * xNDC - 1;
								float yScreen = 1 - 2 * yNDC;
								float xCamera = xScreen * aspectRatio*scaling;
								float yCamera = yScreen * scaling;
								v4f P = v4f(xCamera, yCamera, -1.f, 1.f);
								v4f pWorld = P * inverseView;
								Ray ray = Ray(o, v3f(pWorld));
								colorSum += integrator->render(ray, sampler);
							}
							dataBuffer[y*width + x] += colorSum / sampleNum / progressivePass;
						}

					}
				}
				
			}
			for (int x = 0; x < width; ++x) {
				for (int y = 0; y < height; ++y) {
					integrator->rgb->data[y*width + x] = dataBuffer[y*width + x];
				}
			}
			
		}
		else {
			int width = scene.config.width;
			int height = scene.config.height;
			int sampleNum = scene.config.spp;

			Sampler sampler = Sampler(260563769);
			v3f o = scene.config.camera.o;
			//1
			float cameraPerspective = scene.config.camera.fov;
			glm::mat4 inverseView = glm::lookAt(o, scene.config.camera.at, scene.config.camera.up);
			float aspectRatio = (float)width / (float)height;
			float scaling = tan((scene.config.camera.fov*deg2rad) / 2.f);
			//2
			integrator->rgb->clear();

			for (int x = 0; x < width; ++x) {
				for (int y = 0; y < height; ++y) {
					v3f colorSum = v3f(0.f, 0.f, 0.f);
					if (sampleNum == 1) {
						float xNDC = (x + 0.5) / width;
						float yNDC = (y + 0.5) / height;
						float xScreen = 2 * xNDC - 1;
						float yScreen = 1 - 2 * yNDC;
						float xCamera = xScreen * aspectRatio*scaling;
						float yCamera = yScreen * scaling;
						v4f P = v4f(xCamera, yCamera, -1.f, 1.f);
						v4f pWorld = P * inverseView;
						//v4f eyeTransformed = v4f(o[0],o[1],o[2],1.f)*inverseView;
						Ray ray = Ray(o, v3f(pWorld));
						v3f color = integrator->render(ray, sampler);
					
						integrator->rgb->data[y*width + x] = color;
					}
					else {
						for (int z = 0; z < sampleNum; z++) {
							float xNDC = (x + sampler.next()) / width;
							float yNDC = (y + sampler.next()) / height;
							float xScreen = 2 * xNDC - 1;
							float yScreen = 1 - 2 * yNDC;
							float xCamera = xScreen * aspectRatio*scaling;
							float yCamera = yScreen * scaling;
							v4f P = v4f(xCamera, yCamera, -1.f, 1.f);
							v4f pWorld = P * inverseView;
							Ray ray = Ray(o, v3f(pWorld));
							colorSum += integrator->render(ray, sampler);
						}
						integrator->rgb->data[y*width + x] = colorSum / sampleNum;
					}

				}
			}
		}



	}
}

/**
 * Post-rendering step.
 */
void Renderer::cleanUp() {
    if (realTime) {
        renderpass->cleanUp();
    } else {
        integrator->cleanUp();
    }
}

BSDF::BSDF(const WorldData& d, const Config& c, const size_t matID) : worldData(d), config(c) {
    emission = glm::make_vec3(worldData.materials[matID].emission);
}

Scene::Scene(const Config& config) : config(config) { }

bool Scene::load(bool isRealTime) {
    fs::path file(config.objFile);
    bool ret = false;
    std::string err;

    if (!file.is_absolute())
        file = (config.tomlFile.parent_path() / file).make_preferred();

    tinyobj::attrib_t* attrib_ = &worldData.attrib;
    std::vector<tinyobj::shape_t>* shapes_ = &worldData.shapes;
    std::vector<tinyobj::material_t>* materials_ = &worldData.materials;
    std::string* err_ = &err;
    const string filename_ = file.string();
    const string mtl_basedir_ = file.make_preferred().parent_path().string();
    ret = tinyobj::LoadObj(attrib_, shapes_, materials_, err_, filename_.c_str(), mtl_basedir_.c_str(), true);

    if (!err.empty()) { std::cout << "Error: " << err.c_str() << std::endl; }
    if (!ret) {
        std::cout << "Failed to load scene " << config.objFile << " " << std::endl;
        return false;
    }

    // Build list of BSDFs
    bsdfs = std::vector<std::unique_ptr<BSDF>>(worldData.materials.size());
    for (size_t i = 0; i < worldData.materials.size(); i++) {
        if (worldData.materials[i].illum == 5)
            bsdfs[i] = std::unique_ptr<BSDF>(new MirrorBSDF(worldData, config, i));
        if (worldData.materials[i].illum == 7)
            bsdfs[i] = std::unique_ptr<BSDF>(new DiffuseBSDF(worldData, config, i));
        if (worldData.materials[i].illum != 5 && worldData.materials[i].illum != 7 && worldData.materials[i].illum != 8)
            bsdfs[i] = std::unique_ptr<BSDF>(new PhongBSDF(worldData, config, i));
        if (worldData.materials[i].illum == 8)
            bsdfs[i] = std::unique_ptr<BSDF>(new MixtureBSDF(worldData, config, i));
    }

    // Build list of emitters (and print what has been loaded)
    std::string nbShapes = worldData.shapes.size() > 1 ? " shapes" : " shape";
    std::cout << "Found " << worldData.shapes.size() << nbShapes << std::endl;
    worldData.shapesCenter.resize(worldData.shapes.size());
    worldData.shapesAABOX.resize(worldData.shapes.size());

    for (size_t i = 0; i < worldData.shapes.size(); i++) {
        const tinyobj::shape_t& shape = worldData.shapes[i];
        const BSDF* bsdf = bsdfs[shape.mesh.material_ids[0]].get();
        std::cout << "Mesh " << i << ": " << shape.name << " ["
                  << shape.mesh.indices.size() / 3 << " primitives | ";

        if (bsdf->isEmissive()) {
            Distribution1D faceAreaDistribution;
            float shapeArea = getShapeArea(i, faceAreaDistribution);
            emitters.emplace_back(Emitter{i, shapeArea, bsdf->emission, faceAreaDistribution});
            std::cout << "Emitter]" << std::endl;
        } else {
            std::cout << bsdf->toString() << "]" << std::endl;
        }

        // Build world AABB and shape centers
        worldData.shapesCenter[i] = v3f(0.0);
        for (auto idx: shape.mesh.indices) {
            v3f p = {worldData.attrib.vertices[3 * idx.vertex_index + 0],
                     worldData.attrib.vertices[3 * idx.vertex_index + 1],
                     worldData.attrib.vertices[3 * idx.vertex_index + 2]};
            worldData.shapesCenter[i] += p;
            worldData.shapesAABOX[i].expandBy(p);
            aabb.expandBy(p);
        }
        worldData.shapesCenter[i] /= float(shape.mesh.indices.size());
    }

    // Build BVH
    bvh = std::unique_ptr<TinyRender::AcceleratorBVH>(new TinyRender::AcceleratorBVH(this->worldData));

    const clock_t beginBVH = clock();
    bvh->build();
    std::cout << "BVH built in " << float(clock() - beginBVH) / CLOCKS_PER_SEC << "s" << std::endl;

    return true;
}

float Scene::getShapeArea(const size_t shapeID, Distribution1D& faceAreaDistribution) {
    const tinyobj::shape_t& s = worldData.shapes[shapeID];

    for (size_t i = 0; i < s.mesh.indices.size(); i += 3) {
        const int i0 = s.mesh.indices[i + 0].vertex_index;
        const int i1 = s.mesh.indices[i + 1].vertex_index;
        const int i2 = s.mesh.indices[i + 2].vertex_index;
        const v3f v0{worldData.attrib.vertices[3 * i0 + 0], worldData.attrib.vertices[3 * i0 + 1],
                     worldData.attrib.vertices[3 * i0 + 2]};
        const v3f v1{worldData.attrib.vertices[3 * i1 + 0], worldData.attrib.vertices[3 * i1 + 1],
                     worldData.attrib.vertices[3 * i1 + 2]};
        const v3f v2{worldData.attrib.vertices[3 * i2 + 0], worldData.attrib.vertices[3 * i2 + 1],
                     worldData.attrib.vertices[3 * i2 + 2]};

        const v3f e1{v1 - v0};
        const v3f e2{v2 - v0};
        const v3f e3{glm::cross(e1, e2)};
        faceAreaDistribution.add(0.5f * std::sqrt(e3.x * e3.x + e3.y * e3.y + e3.z * e3.z));
    }
    const float area = faceAreaDistribution.cdf.back();
    faceAreaDistribution.normalize();
    return area;
}

v3f Scene::getFirstLightPosition() const {
    return worldData.shapesCenter[emitters[0].shapeID];
}

v3f Scene::getFirstLightIntensity() const {
    return emitters[0].getRadiance(); // point lights are defined by intensity not radiance
}

float Scene::getShapeRadius(const size_t shapeID) const {
    assert(shapeID < worldData.shapes.size());
    v3f emitterCenter = worldData.shapesCenter[shapeID];
    return worldData.shapesAABOX[shapeID].max.x - emitterCenter.x;
}

v3f Scene::getShapeCenter(const size_t shapeID) const {
    assert(shapeID < worldData.shapes.size());
    return worldData.shapesCenter[shapeID];
}

size_t Scene::getFirstLight() const {
    if (emitters.size() <= 0) return -1;
    return emitters[0].shapeID;
}

v3f Scene::getObjectVertexPosition(size_t objectIdx, size_t vertexIdx) const {
    const tinyobj::attrib_t& sa = worldData.attrib;
    const tinyobj::shape_t& s = worldData.shapes[objectIdx];

    int idx = s.mesh.indices[vertexIdx].vertex_index;
    float x = sa.vertices[3 * idx + 0];
    float y = sa.vertices[3 * idx + 1];
    float z = sa.vertices[3 * idx + 2];
    return v3f(x,y,z);
}

v3f Scene::getObjectVertexNormal(size_t objectIdx, size_t vertexIdx) const {
    const tinyobj::attrib_t& sa = worldData.attrib;
    const tinyobj::shape_t& s = worldData.shapes[objectIdx];

    int idx_n = s.mesh.indices[vertexIdx].normal_index;
    float nx = sa.normals[3 * idx_n + 0];
    float ny = sa.normals[3 * idx_n + 1];
    float nz = sa.normals[3 * idx_n + 2];
    return glm::normalize(v3f(nx,ny,nz));
}

size_t Scene::getObjectNbVertices(size_t objectIdx) const {
    return worldData.shapes[objectIdx].mesh.indices.size();
}

int Scene::getPrimitiveID(size_t vertexIdx) const {
    return vertexIdx / 3;
}

int Scene::getMaterialID(size_t objectIdx, int primID) const {
    return worldData.shapes[objectIdx].mesh.material_ids[primID];
}

TR_NAMESPACE_END
