/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

inline float safeSqrt(float v){
    return std::sqrt(std::max(float(0), v));
}

/**
 * Computes barycentric coordinates.
 */
template<class T>
inline T barycentric(const T& a, const T& b, const T& c, const float u, const float v) {
    return a * (1 - u - v) + b * u + c * v;
}

/**
 * Restricts a value to a given interval.
 */
template<class T>
inline T clamp(T v, T min, T max) {
    return std::min(std::max(v, min), max);
}

/**
 * Checks if vector is zero.
 */
inline bool isZero(const v3f v) {
    return glm::dot(v, v) < Epsilon;
}

/**
 * Generates coordinate system.
 */
inline void coordinateSystem(const v3f& a, v3f& b, v3f& c) {
    if (std::abs(a.x) > std::abs(a.y)) {
        float invLen = 1.f / std::sqrt(a.x * a.x + a.z * a.z);
        c = v3f(a.z * invLen, 0.f, -a.x * invLen);
    } else {
        float invLen = 1.f / std::sqrt(a.y * a.y + a.z * a.z);
        c = v3f(0.f, a.z * invLen, -a.y * invLen);
    }
    b = glm::cross(c, a);
}

/**
 * Converts RGB value to luminance.
 */
inline float getLuminance(const v3f& rgb) {
    return glm::dot(rgb, v3f(0.212671f, 0.715160f, 0.072169f));
}

/**
 * Pseudo-random sampler (Mersenne Twister 19937) structure.
 */
struct Sampler {
    std::mt19937 g;
    std::uniform_real_distribution<float> d;
    explicit Sampler(int seed) {
        g = std::mt19937(seed);
        d = std::uniform_real_distribution<float>(0.f, 1.f);
    }
    float next() { return d(g); }
    p2f next2D() { return {d(g), d(g)}; }
    void setSeed(int seed) {
        g.seed(seed);
        d.reset();
    }
};

/**
 * 1D discrete distribution.
 */
struct Distribution1D {
    std::vector<float> cdf{0};
    bool isNormalized = false;

    inline void add(float pdfVal) {
        cdf.push_back(cdf.back() + pdfVal);
    }

    size_t size() {
        return cdf.size() - 1;
    }

    float normalize() {
        float sum = cdf.back();
        for (float& v : cdf) {
            v /= sum;
        }
        isNormalized = true;
        return sum;
    }

    inline float pdf(size_t i) const {
        assert(isNormalized);
        return cdf[i + 1] - cdf[i];
    }

    int sample(float sample) const {
        assert(isNormalized);
        const auto it = std::upper_bound(cdf.begin(), cdf.end(), sample);
        return clamp(int(distance(cdf.begin(), it)) - 1, 0, int(cdf.size()) - 2);
    }
};


/**
 * Warping functions.
 */
namespace Warp {

	inline v3f squareToUniformSphere(const p2f& sample) {
		v3f v(0.f);
		// TODO: Implement this
		float wz = 1.f - 2.f * sample[0];
		float r = sqrt(1.f - wz * wz);
		float phi = 2 * M_PI *sample[1];
		float wx = r * cos(phi);
		float wy = r * sin(phi);
		v = v3f(wx, wy, wz);
		return v;
	}

	inline float squareToUniformSpherePdf() {
		float pdf = 0.f;
		// TODO: Implement this
		pdf = INV_FOURPI;
		return pdf;
	}

	inline v3f squareToUniformHemisphere(const p2f& sample) {
		v3f v(0.f);
		// TODO: Implement this
		float wz = sample.x;
		float r = sqrt(1.f - wz * wz);
		float phi = 2 * M_PI *sample.y;
		float wx = r * cos(phi);
		float wy = r * sin(phi);
		v = v3f(wx, wy, wz);
		return v;
	}

	inline float squareToUniformHemispherePdf(const v3f& v) {
		float pdf = 0.f;
		// TODO: Implement this
		pdf = INV_TWOPI;
		return pdf;
	}

	inline v2f squareToUniformDiskConcentric(const p2f& sample) {
		v2f v(0.f);
		// TODO: Implement this (optional)
		float phi, r, x, y;
		float a = 2 * sample.x - 1.f;
		float b = 2 * sample.y - 1.f;
		if (a > -b) {
			if (a > b) {
				r = a;
				phi = (M_PI / 4)*(b / a);
			}
			else {
				r = b;
				phi = (M_PI / 4)*(2 - (a / b));
			}
		}
		else {
			if (a < b) {
				r = -a;
				phi = (M_PI / 4)*(4 + b / a);

			}
			else {
				r = -b;
				if (b != 0) {
					phi = (M_PI / 4)*(6 - (a / b));
				}
				else {
					phi = 0;
				}
			}
		}
		x = r * cos(phi);
		y = r * sin(phi);
		v = v2f(x, y);
		return v;
	}

	inline v3f squareToCosineHemisphere(const p2f& sample) {
		v3f v(0.f);
		// TODO: Implement this
		v2f disk = squareToUniformDiskConcentric(sample);
		float z = sqrt(1.f - disk[0] * disk[0] - disk[1] * disk[1]);
		v = v3f(disk[0], disk[1], z);
		return v;
	}

	inline float squareToCosineHemispherePdf(const v3f& v) {
		float pdf = 0.f;
		// TODO: Implement this
		float cosTheta = v.z;
		pdf = cosTheta * INV_PI;
		return pdf;
		return pdf;
	}

	inline v3f squareToPhongLobe(const p2f& sample, const float exponent) {
		v3f v(0.f);

		// TODO: Implement this
		float theta = acos(pow((1 - sample.x), 1 / (exponent + 2)));
		float phi = 2 * M_PI*sample.y;
		float x = sin(theta)*cos(phi);
		float y = sin(theta)*sin(phi);
		float z = cos(theta);
		v = v3f(x, y, z);

		return v;
	}

	inline float squareToPhongLobePdf(const v3f& v, const float exponent) {
		float pdf = 0.f;
		// TODO: Implement this
		pdf = (exponent + 2)*INV_TWOPI*pow(v.z, exponent);
		return pdf;
	}



	inline v3f squareToUniformCone(const p2f& sample, float cosThetaMax) {
		v3f v(0.f);
		// TODO: Implement this
		float cosTheta = (1.f - sample.x) + sample.x*cosThetaMax;
		float sinTheta = glm::sqrt(1.f - cosTheta * cosTheta);
		float phi = sample.y * 2 * M_PI;
		v = v3f(glm::cos(phi)*sinTheta, glm::sin(phi)*sinTheta, cosTheta);
		return v;
	}

	inline float squareToUniformConePdf(float cosThetaMax) {
		float pdf = 0.f;
		// TODO: Implement this
		pdf = 1 * INV_TWOPI / (1 - cosThetaMax);
		return pdf;
	}

	inline v2f squareToUniformTriangle(const p2f& sample) {
		v2f v(0.f);
		float u = std::sqrt(1.f - sample.x);
		v = { 1 - u, u * sample.y };
		return v;
	}

}

TR_NAMESPACE_END