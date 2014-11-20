#ifndef LENS_MODEL_H
#define LENS_MODEL_H

#include <vector>
#include "geometry.h"
#include "transform.h"

#define X_AXIS_INDEX 0
#define Y_AXIS_INDEX 1

#define EYE

class LensModel {
public:
	LensModel();
	LensModel(float radius, float axisPosisiton, float asph, float refr, float aperture);
	void Init(float radius, float axisPosisiton, float asph, float refr, float aperture);
	// Whether a point is on lens
	//bool OnLens(Point &p) const;
	// Refract the ray through this lens
	bool Refract(const Ray &ray, Ray &newRay, const float inRefr, float &outRefr) const;
	// Draw a single lens interfaces
	void Draw(void);
	// Draw a ray shot through lens
	bool DrawRay(Ray ray, Ray &refractedRay, const float inRefr, float &outRefr, float zpos);
	// Refract a ray
	bool RefractRay(const Vector &inDir, const Vector &Normal, float mu, Vector &outDir) const;
	bool SphericalRefract(const Ray &inDir, Ray &outDir) const;
	bool isAir(void) const { return (refr == 1.0f)? true: false; }
	void ChangeUnit();

	float rad;		// Radius for the lens interface
	float asph;		// Aspherity
    float zpos;		// Lens interface position (1D);
	float refr;		// Refraction ratio
	float aper;		// Aperture (diameter) of the interface
	bool isPlane;	// Lens is a plane
	bool isActive;	// Lens is active
private:
	void DrawLens();
	float yToZ(float y);
};

#endif