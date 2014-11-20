#include "lensmodel.h"
#include "lensview.h"
#include "transform.h"

LensModel::LensModel() : rad(0), zpos(0), aper(0), isActive(false)
{
}

LensModel::LensModel(float radius, float axisPosition, float asph, float refr, float aperture) :
	rad(radius), zpos(axisPosition), asph(asph), refr(refr), aper(aperture),
	isActive(true)
{
	this->isPlane = radius == 0? true: false;
}

void LensModel::Init(float radius, float axisPosition, float asph, float refr, float aperture)
{
	this->rad = radius;
	this->zpos = axisPosition;
	this->asph = asph;
	this->refr = refr;
	this->aper = aperture;
	this->isActive = true;
	this->isPlane = radius == 0? true: false;
}

bool OnLens(Point &p, float aper, float rad)
{
    if ((p.y*p.y + p.x*p.x) >= (aper*aper / 4.f + 0.01))
		return false;

	if (abs(p.z) > (abs(rad)))
		return false;

    return true;
}

bool LensModel::RefractRay(const Vector &inDir, const Vector &Normal, float mu2, Vector &outDir) const
{
	Vector N = Normalize(Normal);
	// Change normal direction
	if ((rad < 0) != (inDir.z < 0)) N = -N;
	// mu = n1 / n2
	if (inDir.z < 0)
		mu2 = 1 / mu2;
	float cosI = -Dot(N, inDir);
	float sinT2 = mu2 * mu2 *(1 - cosI * cosI);
	if (sinT2 > 1.0)
		return false;
	float cosT = sqrt(1 - sinT2);
	outDir = mu2 * inDir + (mu2 * cosT - cosT) * N;
	return true;
}

void TransformToLens(const Ray ray, Ray &tRay, float zpos)
{
	tRay.o = ray.o - Vector(0, 0, zpos);
	tRay.o.z = -tRay.o.z;
	tRay.d = ray.d;
	tRay.d.z = -tRay.d.z;
}

void TransformToWorld(Ray &ray, float zpos)
{
	ray.o.z = zpos - ray.o.z;
	ray.d.z = -ray.d.z;
}

bool IntersectPoint(const Ray &tRay, float asph, float rad, float aper, Point &phit)
{
	// Compute equation x^2 + y^2 + z^2 - 2Rz - Qz^2 = 0
	Vector d = tRay.d;
    float A = d.x * d.x + d.y * d.y
              + (1 - asph) * d.z * d.z;
	float B = 2 * (tRay.o.x * d.x + tRay.o.y * d.y + (1 - asph) * tRay.o.z * d.z - rad * d.z);
    float C = tRay.o.x * tRay.o.x + tRay.o.y * tRay.o.y + (1 - asph) * tRay.o.z * tRay.o.z
              - 2 * rad * tRay.o.z;
    // Solve quadratic equation
    float t0, t1;
    if (!Quadratic(A, B, C, &t0, &t1))
    {
        //printf("Ray miss lens.\n");
        return false;
    }
    if (t0 > tRay.maxt)
    {
        //printf("Larger than max time\n");
        return false;
    }
    Point p0 = tRay(t0);
    Point p1 = tRay(t1);
    if (t0 > tRay.mint && OnLens(p0, aper, rad))
        phit = p0;
    else if (t1 < tRay.maxt && OnLens(p1, aper, rad))
        phit = p1;
    else
    {
        //printf("Intersection not on lens.\n");
        return false;
    }
}

bool IntersectPlane(Ray ray, Point &phit, float zpos, float aper)
{
	float thit = (zpos - ray.o.z) / ray.d.z;
	if (thit > ray.maxt)
		return false;
	phit = ray(thit);
	if (fabs(phit.y) > aper / 2)
	{
		//printf("Ray is blocked by aperture\n");
		return false;
	}
	return true;
}

bool LensModel::Refract(const Ray &ray, Ray &newRay, const float inRefr, float &outRefr) const
{
	outRefr = this->refr;
	// hit point
	Point phit;
	// Transform ray to len's Euclidean space
	Ray tRay = ray;
	TransformToLens(ray, tRay, zpos);

	// deal with the plane lens
	if (this->isPlane)
	{
		if (IntersectPlane(tRay, phit, 0, aper))
			return true;
		Vector N = Vector(0, 0, 1);
		newRay.o = phit;
		RefractRay(Normalize(ray.d), N, inRefr / outRefr, newRay.d);
		// transform ray back to world space
		TransformToWorld(newRay, zpos);
		return true;
	}
	// deal with normal lens
    if (!IntersectPoint(tRay, asph, rad, aper, phit))
		return false;
	if (rad > 0)
	{

	}

    //Vector N = Vector(phit.x, phit.y, phit.z - zpos + rad + asph);
	Vector N = Vector(phit.x, phit.y, (1 - asph) * phit.z - rad);
    N = Normalize(N);

	newRay.o = phit;
	//printf("mu: %f calculate: %f\n", mu, inRefr / outRefr);
	if (!RefractRay(Normalize(tRay.d), N, inRefr / outRefr, newRay.d))
		return false;
	// transform ray back to world space
	TransformToWorld(newRay, zpos);
	//drawLine(newRay.o.z + N.z, newRay.o.y - N.y, newRay.o.z, newRay.o.y);
    return true;
}

void LensModel::Draw()
{
	glColor3f(0.75f, 0.75f, 0.75f); glLineWidth(1);
	if (this->isPlane)
	{
		drawLine(zpos, -this->aper / 2, zpos, this->aper / 2);
		return;
	}
	float nextX, nextY, thisX, thisY;
	thisY = -this->aper / 2.0f;
	float rs = this->rad / fabs(this->rad);
	int i;
  
	float step = this->aper / LVSD;	
	thisX = yToZ(thisY);
	for (i = 0; i < LVSD; i++) {
		nextY = thisY + step;
		nextX = yToZ(nextY);
		drawLine(thisX, thisY, nextX, nextY);
		thisX = nextX; thisY = nextY;
	}
}

float LensModel::yToZ(float y)
{
	float posz1, posz2;
	if (!Quadratic(1 - asph, -2 * rad, y * y, &posz1, &posz2))
		return 1.f;
	return zpos - ((rad > 0)? posz1: posz2);
}

bool LensModel::DrawRay(Ray ray, Ray& refractedRay, const float inRefr, float &outRefr, float zpos)
{
	if (!Refract(ray, refractedRay, inRefr, outRefr))
		return false;
	if (refractedRay.o.z < zpos - 3)
		return false;
	drawLine(ray.o.z, ray.o.y, refractedRay.o.z, refractedRay.o.y);
	// in ray
	
	return true;
}

void LensModel::ChangeUnit()
{
	rad *= 0.001;
    zpos *= 0.001;
	aper *= 0.001;
}