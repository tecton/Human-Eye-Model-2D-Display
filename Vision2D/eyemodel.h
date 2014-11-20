#ifndef EYE_MODEL_H
#define EYE_MODEL_H

#include <vector>
#include "lensmodel.h"
#define TEST_RAY_RIGHT

class EyeModel
{
public:
	EyeModel();
	~EyeModel();
	// load model data file
	bool Load(const char *file);
	void SetDP(double A);
	float MaxAper();
	bool TraceLeftRay(Ray inRay, Ray &outRay) const;
	bool TraceRay(Ray inRay, Ray &outRay) const;
	void Draw();
	void MoveUp() {testYPos += 10;}//printf("y: %f\n", testYPos);}
	void MoveDown() {testYPos -= 10;}//printf(": %f\n", testYPos);}
	void MoveLeft() {testZPos += 50;}//printf("zpos: %f\n", testZPos + imagePlane);}
	void MoveRight() {testZPos -= 50;}//printf("zpos: %f\n", testZPos + imagePlane);}
	void MoreDP() {dp += 0.2; if (dp > 12) dp = 12; SetDP(dp); }
	void LessDP() {dp -= 0.2; if (dp < 0) dp = 0; SetDP(dp); }
	vector<LensModel> lenses;
	void ChangeUnit();
private:
	void testBalance();
	void DrawRay(Ray &ray);
	void DrawTestRays();
	void DrawLeftTestRays();
	void DrawGPURays();
	float testYPos, testZPos;
	float imagePlane;
	float pupilZpos;
	float exitPupilRad;
	float dp;
};

#endif