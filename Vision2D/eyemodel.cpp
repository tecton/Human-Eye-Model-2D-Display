#include <fstream>
#include <iostream>
#include "lensview.h"
#include "eyemodel.h"

using namespace std;

EyeModel::EyeModel()
{
	testYPos = 30;
	testZPos = 800;
	pupilZpos = -3.66;
	exitPupilRad = 1.13;
	imagePlane = -38;
	dp = 0;
}

EyeModel::~EyeModel()
{
	lenses.clear();
}

bool EyeModel::Load(const char* file)
{
	ifstream eyeData(file);
	if (!eyeData)
		return false;

	this->lenses.clear();
	float lastRefractive = 1.f;
	float axisPos = 0.f;
	float radius;
	float asphericity;
	float thickness;
	float refractive;
	float aperture;
	while (!eyeData.eof())
	{
		LensModel lensModel;
		string data;
		eyeData >> data;
		// Remove comments and empty lines
		if (data[0] == '#' || !data.compare(""))
		{
			getline(eyeData, data);
			continue;
		}
		radius = (float)atof(data.c_str());
		eyeData >> data;
		asphericity = (float)atof(data.c_str());
		eyeData >> data;
		thickness = (float)atof(data.c_str());
		eyeData >> data;
		if (data == "P")
			refractive = -1;
		else if (data == "A")
			refractive = -2;
		else
			refractive = (float)atof(data.c_str());
		eyeData >> data;
		aperture = (float)atof(data.c_str());
		lensModel.Init(radius, axisPos, asphericity,
					   refractive, aperture);
		this->lenses.push_back(lensModel);
		lastRefractive = refractive;
		axisPos -= thickness;
	}
	//this->imagePlane = this->lenses[this->lenses.size() - 1].zpos
	//				   - thickness;
	//CalculateImagePlane();
	//CalculateExitPupil();
	//testBalance();
	SetDP(dp);
	return true;
}

void EyeModel::SetDP(double A)
{
	lenses[2].rad = 10.2 - 1.75 * log(A + 1);
	lenses[3].rad = -6 + 0.2294 * log(A + 1);
	lenses[2].zpos = lenses[1].zpos - (3.05 - 0.05 * log(A + 1));
	lenses[3].zpos = lenses[2].zpos - (4 + 0.1 * log(A + 1));
	lenses[4].zpos = lenses[3].zpos - 16.4;
	lenses[2].refr = 1.42 + 9e-5 * (10 * A + A * A);
	lenses[2].asph = -3.1316 - 0.34 * log(A + 1);
	lenses[3].asph = -1 - 0.125 * log(A + 1);
}

float EyeModel::MaxAper()
{
	float aper = 0.0f;
	vector<LensModel>::iterator lensIter;
	for (lensIter = lenses.begin(); lensIter != lenses.end(); lensIter++)
		if (lensIter->aper > aper)
			aper = lensIter->aper;
	return aper;
}

void EyeModel::Draw()
{
	// Draw lenses
	vector<LensModel>::iterator lensIter;
	for (lensIter = lenses.begin(); lensIter != lenses.end(); lensIter++)
		lensIter->Draw();
#ifdef TEST_RAY_RIGHT
	DrawLeftTestRays();
#else
	DrawLeftTestRays();
#endif // TEST_RAY
}

void EyeModel::DrawRay(Ray &ray)
{
	drawLine(ray.o.z, ray.o.y, ray.o.z + ray.d.z * 500, ray.o.y + ray.d.y * 500);
}

void EyeModel::DrawLeftTestRays()
{
	//draw rays
	//float y = exitPupilRad / 4;
	float y = 4;
	float step = -2 * y / LVSD;
	glColor3f(0, 0, 1);
	Ray inRay, outRay;
	float inRefr = 1, outRefr = 1;
	float g, b;
	// from left to right
	for (int i = 0; i < LVSD; ++i)
	{
		int j = 0;
		float g, b;
		//inRay = Ray(Point(0, y + step * i, 10), Vector(0, 0, -1));
		inRefr = 1;
		inRay.o = Point(0, testYPos, testZPos);
		inRay.d = Point(0, y + step * i, 0) - inRay.o;
		//printf("ray %f, %f\n", y + step * i, 100);
		//outfile << "ray " << y + step * i << " " << "100\n";
		//printf("new ray\n");
		inRay.d /= inRay.d.Length();
		if (inRay.o.z > 10)
		{
			double t = (Point(0, 0, 10) - inRay.o).z / inRay.d.z;
			inRay.o = inRay(t);
		}
		for (; j < lenses.size(); ++j)
		//for (; j < 1; ++j)
		{
			g = (float)(j % 2);
			b = (float)((j + 1) % 2);
			glColor3f(0, g, b);
			if (lenses[j].DrawRay(inRay, outRay, inRefr, outRefr, lenses[j].zpos))
			{
				if (j == 2 && outRay.o.z < lenses[2].zpos - 1)
					break;
				printf("z: %f\n", outRay.o.z);
				if (lenses[j].rad > 0)
				{
					if (outRay.o.z < lenses[j].zpos -3)
						break;
				}
				inRay = outRay;
				inRefr = outRefr;
			}
			else
				break;
		}
		printf("\n");
		if (j == lenses.size())
		{
			glColor3f(0, b, g);
			DrawRay(inRay);
			//printf("focus point %f\n", -inRay.o.y / inRay.d.y * inRay.d.z + inRay.o.z);
		}
	}
}

void EyeModel::DrawTestRays()
{
	////draw rays
	//float step = -exitPupilRad / LVSD;
	//float y = exitPupilRad;
	//glColor3f(0, 0, 1);
	//Ray inRay, outRay;
	//float g, b;

	////y = this->lenses[this->lenses.size() - 1].aper / 2;
	//y = this->exitPupilRad;
	//step = -(float)(2 * y / LVSD);
	//Point imagePoint(0, testYPos, this->imagePlane + this->testZPos);
	//for (int i = 0; i < LVSD; ++i)
	//{
	//	int j = this->lenses.size() - 1;
	//	//Point lensPoint(0, y + step * i, this->lenses[this->lenses.size() - 1].zpos);
	//	Point lensPoint(0, y + step * i, pupilZpos);
	//	Vector direction = lensPoint - imagePoint;
	//	inRay = Ray(imagePoint, direction);
	//	for (;j >= 0; --j)
	//	{
	//		g = (float)(j % 2);
	//		b = (float)((j + 1) % 2);
	//		glColor3f(0, g, b);
	//		if (lenses[j].DrawRay(inRay, outRay))
	//			inRay = outRay;
	//		else
	//			break;
	//	}
	//	if (j == -1)
	//	{
	//		glColor3f(0, b, g);
	//		DrawRay(inRay);
	//		//printf("%f %f\n", inRay.d.y, inRay.d.z);
	//		float time = (1500 - inRay.o.z) /inRay.d.z;
	//		Point p = inRay(time);
	//		printf("%f %f %f\n", p.x, p.y, p.z);
	//	}
	//	else
	//		printf("miss\n");
	//}
}

void EyeModel::DrawGPURays()
{
	////draw rays
	//float step = 10.0 / LVSD;
	//float y = 0;
	//glColor3f(0, 0, 1);
	//Ray inRay, outRay;
	//float g, b;

	////Point imagePoint(0, testYPos, this->imagePlane + this->testZPos);
	//Point lensPoint(0, 0, this->pupilZpos);
	//for (int i = 0; i < LVSD; ++i)
	//{
	//	int j = this->lenses.size() - 1;
	//	//Point lensPoint(0, y + step * i, this->lenses[this->lenses.size() - 1].zpos);
	//	Point imagePoint(0, y + step * i, this->imagePlane);
	//	Vector direction = lensPoint - imagePoint;
	//	inRay = Ray(imagePoint, direction);
	//	for (;j >= 0; --j)
	//	{
	//		g = (float)(j % 2);
	//		b = (float)((j + 1) % 2);
	//		glColor3f(0, g, b);
	//		if (lenses[j].DrawRay(inRay, outRay))
	//			inRay = outRay;
	//		else
	//			break;
	//	}
	//	if (j == -1)
	//	{
	//		glColor3f(0, b, g);
	//		DrawRay(inRay);
	//	}
	//}
}

bool EyeModel::TraceRay(Ray inRay, Ray &outRay) const
{
//	Ray ray = inRay;
//	for (int i = this->lenses.size() - 1; i >= 0; --i)
//	{
//#ifndef EYE
//		if (this->lenses[i].SphericalRefract(ray, outRay))
//#else
//		if (this->lenses[i].Refract(ray, outRay))
//#endif
//			ray = outRay;
//		else
//			return false;
//	}
//	outRay = ray;
	return true;
}

bool EyeModel::TraceLeftRay(Ray inRay, Ray &outRay) const
{
//	Ray ray = inRay;
//	for (int i = 0; i < this->lenses.size(); ++i)
//	{
//#ifndef EYE
//		if (this->lenses[i].SphericalRefract(ray, outRay))
//#else
//		if (this->lenses[i].Refract(ray, outRay))
//#endif
//			ray = outRay;
//		else
//			return false;
//	}
//	outRay = ray;
	return true;
}

void EyeModel::ChangeUnit()
{
	testYPos *= 0.001;
	testZPos *= 0.001;
	imagePlane *= 0.001;
	pupilZpos *= 0.001;
	exitPupilRad *= 0.001;
	for (int i = 0; i < lenses.size(); ++i)
	{
		lenses[i].ChangeUnit();
	}
}