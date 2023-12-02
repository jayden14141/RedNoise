#include <algorithm>
#include <CanvasTriangle.h>
#include <CanvasPoint.h>
#include <Colour.h>
#include <DrawingWindow.h>
#include <map>
#include <ModelTriangle.h>
#include <TextureMap.h>
#include <RayTriangleIntersection.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <cmath>

#define WIDTH 320
#define HEIGHT 240
#define PI 3.14159265358979323846

std::map<std::string, Colour> colourMap;
TextureMap tM;
std::vector<ModelTriangle> modelTriangles;
std::vector<ModelTriangle> mirrors;
std::vector<ModelTriangle> glasses;
glm::vec3 cameraPosition = glm::vec3(0.0, 0.0, 6.0);
float focalLength = 2.0;
glm::mat3 cameraOrientation = glm::mat3(1.0);
bool orbits = false;
bool lightMove = false;
int draw_mode = 2;
int light_mode = 2;
glm::vec3 light(0, 0.2, 0.7);
std::vector<glm::vec3> lights;
bool shadow_soft = true;
bool mirror_mode = true;
bool glass_mode = true;

void sort(bool yAxis, CanvasTriangle &t) {
	if (yAxis) {
		if (t[0].y > t[1].y) std::swap(t[0], t[1]);
		if (t[0].y > t[2].y) std::swap(t[0], t[2]);
		if (t[1].y > t[2].y) std::swap(t[1], t[2]); 
	} else {
		if (t[0].x > t[1].x) std::swap(t[0], t[1]);
		if (t[0].x > t[2].x) std::swap(t[0], t[2]);
		if (t[1].x > t[2].x) std::swap(t[1], t[2]); 
	}
}

bool compareVertices(glm::vec3 &v1, glm::vec3 &v2) {
	return ((v1.x == v2.x) && (v1.y == v2.y) && (v1.z == v2.z));
}

void initialize_lights() {
	float radius = 1.0f;
		for (int j = 1; j <= 2; j++) {
			for (int i = 0; i < 8; i++) {
			float angle = PI / 4;
			float x = light.x + j * radius * glm::cos(i * angle);
			float y = light.y;
			float z = light.z + j * radius * glm::sin(i * angle);
			lights.push_back(glm::vec3(x,y,z));
			}
		}
		for (int j = 1; j <= 2; j++) {
			for (int i = 0; i < 8; i++) {
			float angle = PI / 4;
			float x = light.x + j * radius * glm::cos(i * angle);
			float y = light.y + j * radius * glm::sin(i * angle);
			float z = light.z;
			lights.push_back(glm::vec3(x,y,z));
			}
		}
}

// Adds/deletes 8 lightsource having the 'light' as the core
void change_lightSize(bool expand) {
	if (expand) {
		int steps = (lights.size() - 1) / 8;
		float radius = (steps + 1) * 1.0f;
		for (int i = 0; i < 8; i++) {
			float angle = PI / 4;
			float x = light.x + radius * glm::cos(i * angle);
			float y = light.y;
			float z = light.z + radius * glm::sin(i * angle);
			lights.push_back(glm::vec3(x,y,z));
		}
		for (int i = 0; i < 8; i++) {
			float angle = PI / 4;
			float x = light.x + radius * glm::cos(i * angle);
			float y = light.y + radius * glm::sin(i * angle);
			float z = light.z;
			lights.push_back(glm::vec3(x,y,z));
		}
		
	} else {
		if (lights.size() > 1) {
			for (int i = 0; i < 16; i++) {
				lights.pop_back();
			}
		}
	}
}

void lookAt() {
	glm::vec3 forward = glm::normalize(cameraPosition);
	glm::vec3 right = glm::normalize(glm::cross(glm::vec3(0,1,0), forward));
	glm::vec3 up = glm::normalize(glm::cross(forward, right));

	cameraOrientation[0] = right;
	cameraOrientation[1] = up;
	cameraOrientation[2] = forward;
}

void change_orientation(bool x_axis, float degree) {
	glm::mat3 m;
	if (x_axis) {
		m = glm::mat3(
   			1, 0, 0,
   			0, cos(degree), sin(degree),
   			0, -sin(degree), cos(degree));
	} else {
		m = glm::mat3(
   			cos(degree), 0, -sin(degree),
  			0, 1, 0,
   			sin(degree), 0, cos(degree));
	}
	cameraOrientation = cameraOrientation * m;
}

void rotate_camera(bool x_axis, float degree) {
	glm::mat3 m;
	if (x_axis) {
		m = glm::mat3(
   			1, 0, 0,
   			0, cos(degree), sin(degree),
   			0, -sin(degree), cos(degree));
	} else {
		m = glm::mat3(
   			cos(degree), 0, -sin(degree),
  			0, 1, 0,
   			sin(degree), 0, cos(degree));
	}
	cameraPosition = cameraPosition * m;
}

void translate_light(int where, bool positive) {

	//translate by x
	if (where == 0) {
		if (positive) {
			light += glm::vec3(0.1, 0, 0);
			for (glm::vec3 l : lights) l += glm::vec3(0.1, 0, 0);
		} else {
			light -= glm::vec3(0.1, 0, 0);
			for (glm::vec3 l : lights) l -= glm::vec3(0.1, 0, 0);
		}
	//translate by y
	} else if (where == 1) {
		if (positive) {
			light += glm::vec3(0, 0.1, 0);
			for (glm::vec3 l : lights) l += glm::vec3(0, 0.1, 0);
		} else {
			light -= glm::vec3(0, 0.1, 0);
			for (glm::vec3 l : lights) l -= glm::vec3(0, 0.1, 0);
		}
	//translate by z
	} else {
		if (positive) {
			light += glm::vec3(0, 0, 0.1);
			for (glm::vec3 l : lights) l += glm::vec3(0, 0, 0.1);
		} else {
			light -= glm::vec3(0, 0, 0.1);
			for (glm::vec3 l : lights) l += glm::vec3(0, 0, 0.1);
		}
	}
	std::cout << light.x << ", " << light.y << ", " << light.z << std::endl;
}

void translate_camera(int where, bool positive) {

	//translate by x
	if (where == 0) {
		if (positive) {
			cameraPosition += glm::vec3(0.2, 0, 0);
		} else {
			cameraPosition -= glm::vec3(0.2, 0, 0);
		}
	//translate by y
	} else if (where == 1) {
		if (positive) {
			cameraPosition += glm::vec3(0, 0.2, 0);
		} else {
			cameraPosition -= glm::vec3(0, 0.2, 0);
		}
	//translate by z
	} else {
		if (positive) {
			cameraPosition += glm::vec3(0, 0, 0.2);
		} else {
			cameraPosition -= glm::vec3(0, 0, 0.2);
		}
	}
}

void orbit() {
	if (orbits) {
		rotate_camera(false, -PI / 200);
		lookAt();
	}
}

// Interpolation for colour (rgb) => vec3
std::vector<glm::vec3> interpolateThreeElementValues (glm::vec3 from, glm::vec3 to, int numberOfValues) {
	std::vector<glm::vec3> v;
	glm::vec3 dif = to - from;
	glm::vec3 unit = dif / (float)(numberOfValues - 1);
	// Scalar * vec3 is not accepted.
	for(int i = 0; i < numberOfValues; i++) {
		glm::vec3 a((float)i, (float)i, (float)i);
		v.push_back(from + unit * a);
	}
	return v;
}

std::vector<CanvasPoint> interpolateCanvasPoint (CanvasPoint from, CanvasPoint to, int numberOfValues) {
	float xUnit = (to.x - from.x) / (float)(numberOfValues - 1);
	float yUnit = (to.y - from.y) / (float)(numberOfValues - 1);
	std::vector<CanvasPoint> v;
	if (from.depth == 0 & to.depth == 0) {
		for(int i = 0; i < numberOfValues; i++) v.push_back(CanvasPoint(from.x+xUnit*i, from.y+yUnit*i));
		return v;
	} else {
		float zUnit = (1 / to.depth - 1 / from.depth) / (float)(numberOfValues - 1);
		for(int i = 0; i < numberOfValues; i++) v.push_back(CanvasPoint(from.x+xUnit*i, from.y+yUnit*i, 1 / (1 / from.depth+zUnit*i)));
		return v;
	}
}

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
	float unit = (to - from) / (float)(numberOfValues - 1);
	std::vector<float> v;
	for(int i = 0; i < numberOfValues; i++) v.push_back(from + unit * i);
	return v;
}

// Using struct CanvasPoint
void line(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour c) {
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float numberOfSteps = std::max(abs(xDiff), abs(yDiff));
	float xStepSize = xDiff / numberOfSteps;
	float yStepSize = yDiff / numberOfSteps;
	uint32_t colour = (255 << 24) + (c.red << 16) + (c.green << 8) + c.blue;

	for (float i = 0.0; i <= numberOfSteps; i++) {
		float x = from.x + (xStepSize * i);
		float y = from.y + (yStepSize * i);
		window.setPixelColour(round(x), round(y), colour);
	}
}

void line(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour c, std::vector<std::vector<float>> &depthBuffer) {
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float zDiff = to.depth - from.depth;
	float numberOfSteps = std::max(abs(xDiff), abs(yDiff));
	numberOfSteps = std::max(numberOfSteps, zDiff);
	float xStepSize = xDiff / numberOfSteps;
	float yStepSize = yDiff / numberOfSteps;
	float zStepSize = zDiff / numberOfSteps;


	uint32_t colour = (255 << 24) + (c.red << 16) + (c.green << 8) + c.blue;

	for (float i = 0.0; i < numberOfSteps; i++) {
		float x = from.x + (xStepSize * i);
		float y = from.y + (yStepSize * i);
		float z = from.depth + (zStepSize * i);

		int x_int = round(x);
		int y_int = round(y);

		if (x_int >= 0 && x_int < WIDTH && y_int >= 0 && y_int < HEIGHT) {
			if (depthBuffer[round(x)][round(y)] <= 1.0 / z) {
			window.setPixelColour(round(x), round(y), colour);
			depthBuffer[round(x)][round(y)] = 1.0 / z;
			}
		}
	}
}

void fill(DrawingWindow &window, bool bottomFlat, CanvasPoint left, CanvasPoint right, CanvasPoint p, Colour c) {
	int numberOfVal = abs(left.y - p.y) + 1;
	if(bottomFlat) {
		std::vector<float> v1 = interpolateSingleFloats(p.x, left.x, numberOfVal);
		std::vector<float> v2 = interpolateSingleFloats(p.x, right.x, numberOfVal);
		int i = 0;
		for(float y = p.y; y < left.y; y++) {
			line(window, CanvasPoint(v1[i], y), CanvasPoint(v2[i], y), c);
			i++;
		}
	} else {
		std::vector<float> v1 = interpolateSingleFloats(left.x, p.x, numberOfVal);
		std::vector<float> v2 = interpolateSingleFloats(right.x, p.x, numberOfVal);
		int i = 0;
		for(float y = left.y; y < p.y; y++) {
			line(window, CanvasPoint(v1[i], y), CanvasPoint(v2[i], y), c);
			i++;
		}
	}
}


// Override function having a depth buffer
void fill(DrawingWindow &window, bool bottomFlat, CanvasPoint left, CanvasPoint right, CanvasPoint p, 
		Colour c, std::vector<std::vector<float>> &depthBuffer) {
	int numberOfVal = abs(left.y - p.y) + 1;
	if(bottomFlat) {
		std::vector<float> v1 = interpolateSingleFloats(p.x, left.x, numberOfVal);
		std::vector<float> v2 = interpolateSingleFloats(p.x, right.x, numberOfVal);
		std::vector<float> v1_z = interpolateSingleFloats(p.depth, left.depth, numberOfVal);
		std::vector<float> v2_z = interpolateSingleFloats(p.depth, right.depth, numberOfVal);
		int i = 0;
		for(float y = p.y; y < left.y; y++) {
			line(window, CanvasPoint(v1[i], y, v1_z[i]), CanvasPoint(v2[i], y, v2_z[i]), c, depthBuffer);
			i++;
		}
	} else {
		std::vector<float> v1 = interpolateSingleFloats(left.x, p.x, numberOfVal);
		std::vector<float> v2 = interpolateSingleFloats(right.x, p.x, numberOfVal);
		std::vector<float> v1_z = interpolateSingleFloats(left.depth, p.depth, numberOfVal);
		std::vector<float> v2_z = interpolateSingleFloats(right.depth, p.depth, numberOfVal);
		int i = 0;
		for(float y = left.y; y < p.y; y++) {
			line(window, CanvasPoint(v1[i], y, v1_z[i]), CanvasPoint(v2[i], y, v2_z[i]), c, depthBuffer);
			i++;
		}
	}
}

std::vector<uint32_t> textureColour (TextureMap tm, CanvasPoint from, CanvasPoint to, int numberOfValues) {
	std::vector<CanvasPoint> texturePoint = interpolateCanvasPoint(from, to, numberOfValues);
	std::vector<uint32_t> colour;
	for (int i = 0; i < numberOfValues; i++) {
		colour.push_back(tm.pixels[round(texturePoint[i].x) + round(texturePoint[i].y) * tm.width]);
	}
	return colour;
}


// All the rake rows are passed as v1 and v2
void fillTexture(DrawingWindow &window, TextureMap tm, CanvasTriangle texture, CanvasTriangle canvas) {
	int numberOfRow = abs(canvas[0].y - canvas[1].y) + 1;

	std::vector<CanvasPoint> c_left = interpolateCanvasPoint(canvas[0], canvas[1], numberOfRow);
	std::vector<CanvasPoint> c_right = interpolateCanvasPoint(canvas[0], canvas[2], numberOfRow);

	std::vector<CanvasPoint> t_left = interpolateCanvasPoint(texture[0], texture[1], numberOfRow);
	std::vector<CanvasPoint> t_right = interpolateCanvasPoint(texture[0], texture[2], numberOfRow);
	
	for (int i = 0; i < numberOfRow; i++) {
			int numberOfCol = round(c_right[i].x - c_left[i].x);
			std::vector<uint32_t> colour = textureColour(tm, t_left[i], t_right[i], numberOfCol);
			int start = round(c_left[i].x);
			for (int j = 0; j < numberOfCol; j++) {
				window.setPixelColour(start++, (c_left[i].y), colour[j]);
			}
	}
}

// Helper function that specifies left and right line within the pointer
void leftAndRight(CanvasTriangle &t, CanvasPoint &left, CanvasPoint &right) {

	sort(true, t);
	// yDiff : xDiff = y1 - y0: alpha - x0
	// alpha = (xDiff * (y1-y0) / yDiff) + x0
	// alpha = xDiff * proportion + x0         (proportion = (y1-y0) / yDiff)
	// zeta = zDiff * proportion + z0
	float xDiff = t[2].x - t[0].x;
	float yDiff = t[2].y - t[0].y;
	float zDiff = t[2].depth - t[0].depth;
	float proportion = (t[1].y-t[0].y) / yDiff;
	float alpha = (xDiff * proportion) + t[0].x;
	float zeta = (zDiff * proportion) + t[0].depth;

	if (t[1].x < alpha) {
		left = t[1];
		right = CanvasPoint(alpha, t[1].y, zeta);
	} else {
		left = CanvasPoint(alpha, t[1].y, zeta);
		right = t[1];
	}
}

void unfilledTriangle(DrawingWindow &window) {
	CanvasPoint v0 = CanvasPoint(rand() % window.width - 1, rand() % window.height - 1);
	CanvasPoint v1 = CanvasPoint(rand() % window.width - 1, rand() % window.height - 1);
	CanvasPoint v2 = CanvasPoint(rand() % window.width - 1, rand() % window.height - 1);
	CanvasTriangle t = CanvasTriangle(v0, v1, v2);

	Colour c = Colour(rand() % 256, rand() % 256, rand() % 256);
	line(window, t.v0(), t.v1(), c);
	line(window, t.v1(), t.v2(), c);
	line(window, t.v2(), t.v0(), c);
}

void unfilledTriangle(DrawingWindow &window, CanvasTriangle &t, Colour colour) {
	line(window, t.v0(), t.v1(), colour);
	line(window, t.v1(), t.v2(), colour);
	line(window, t.v2(), t.v0(), colour);
}


void filledTriangle(DrawingWindow &window) {
	CanvasPoint v0 = CanvasPoint(rand() % window.width - 1, rand() % window.height - 1);
	CanvasPoint v1 = CanvasPoint(rand() % window.width - 1, rand() % window.height - 1);
	CanvasPoint v2 = CanvasPoint(rand() % window.width - 1, rand() % window.height - 1);
	CanvasTriangle t = CanvasTriangle(v0, v1, v2);
	
	CanvasPoint left, right;
	leftAndRight(t, left, right);

	Colour c = Colour(rand() % 256, rand() % 256, rand() % 256);
	fill(window, true, left, right, t[0], c);
	fill(window, false, left, right, t[2], c);
	// line(window, left, right, Colour(255, 255, 255));

	unfilledTriangle(window, t, Colour(255, 255, 255));
}

// Overriding the function that has designated colour & coordinates
void filledTriangle(DrawingWindow &window, CanvasTriangle &t, Colour c, std::vector<std::vector<float>> &depthBuffer) {
	CanvasPoint left, right;
	leftAndRight(t, left, right);

	fill(window, true, left, right, t[0], c, depthBuffer);
	fill(window, false, left, right, t[2], c, depthBuffer);
}

void textureMapping(DrawingWindow &window, CanvasTriangle texture, CanvasTriangle canvas) {
	CanvasPoint tLeft, tRight;

	CanvasPoint cLeft, cRight;
	leftAndRight(canvas, cLeft, cRight);
	leftAndRight(texture, tLeft, tRight);
	
	float txDiff = texture[2].x - texture[0].x;
	float tyDiff = texture[2].y - texture[0].y;
	// Proportion = length from top to left(right) / length of whole side
	float proportion = (cRight.x - canvas[0].x) / (canvas[2].x - canvas[0].x);
	// When rake row starts(ends) in left vertex
	if(cLeft.x == canvas[1].x) {
		tLeft = texture[1];
		tRight = CanvasPoint(texture[0].x + txDiff * proportion, texture[0].y + tyDiff * proportion);
	} else {
		tRight = texture[1];
		tLeft = CanvasPoint(texture[0].x + txDiff * proportion, texture[1].y + tyDiff * proportion);
	}

	fillTexture(window, tM, CanvasTriangle(texture[0], tLeft, tRight), CanvasTriangle(canvas[0], cLeft, cRight));
	fillTexture(window, tM, CanvasTriangle(texture[2], tLeft, tRight), CanvasTriangle(canvas[2], cLeft, cRight));
}

CanvasPoint getCanvasIntersectionPoint(glm::vec3 &cameraPosition, glm::vec3 &vertexPosition, float focalLength) {

	glm::vec3 direction = vertexPosition - cameraPosition;
	direction = direction * cameraOrientation;

	float canvasX = HEIGHT * (focalLength * -direction[0] / direction[2]) + WIDTH / 2;
	float canvasY = HEIGHT * (focalLength * direction[1] / direction[2]) + HEIGHT / 2;
	float depth =  -direction[2];
	return CanvasPoint(canvasX, canvasY, depth);
}

// Returns vectors of vertexNormals of the given RayTriangleIntersection
std::vector<glm::vec3> vertexNormal(RayTriangleIntersection &rti) {
	std::vector<glm::vec3> vertexNormals;
	for (int i = 0; i < 3; i++) {
		glm::vec3 v(0,0,0);
		// Makes it count itself as well
		int count = 0;
		for (ModelTriangle compare : modelTriangles) {
			for (int j = 0; j < 3; j++) {
				if (compareVertices(rti.intersectedTriangle.vertices[i], compare.vertices[j])) {
					count++;
					v += compare.normal;
				}
			}
		}
		v /= (float)count;
		v = glm::normalize(v);
		vertexNormals.push_back(v);
	}
	return vertexNormals;
}

// Returns barycentric coordinate of point P which is intersection point of rti
glm::vec3 barycentric(RayTriangleIntersection &rti) {
	glm::vec3 v0 = rti.intersectedTriangle.vertices[1] - rti.intersectedTriangle.vertices[0];
	glm::vec3 v1 = rti.intersectedTriangle.vertices[2] - rti.intersectedTriangle.vertices[0];
	glm::vec3 v2 = rti.intersectionPoint - rti.intersectedTriangle.vertices[0];

    float d00 = glm::dot(v0, v0);
    float d01 = glm::dot(v0, v1);
    float d11 = glm::dot(v1, v1);
    float d20 = glm::dot(v2, v0);
    float d21 = glm::dot(v2, v1);

    float size = d00 * d11 - d01 * d01;
    float u = (d11 * d20 - d01 * d21) / size;
    float v = (d00 * d21 - d01 * d20) / size;
    float w = (float)1 - u - v;
	return glm::vec3(w, u, v);	
}

// Returns a directional vector from the 2D canvaspoint to the camera Position
glm::vec3 getDirectionVector(CanvasPoint canvasPosition) {
	// glm::vec3 direction;
	// Bug : - Camera on 1/4 th quadrant : direction = glm::vec3(canvasPosition.x, canvasPosition.y, 0) + cameraPosition;
	// 		 - Camera on 2/3 th quadrant : direction = glm::vec3(canvasPosition.x, canvasPosition.y, 0) - cameraPosition;
	//		 - Having a large z value solves the issue (Why..?)
	glm::vec3 direction =  glm::vec3(canvasPosition.x, canvasPosition.y, 100) + cameraPosition;
	float dZ = -1 * direction[2];
	float dX = -dZ * (direction[0] - WIDTH / 2) / (HEIGHT * focalLength);
	float dY = dZ * (direction[1] - HEIGHT / 2) / (HEIGHT * focalLength);
	return glm::normalize(glm::vec3(dX, dY, dZ));
}

RayTriangleIntersection getClosestIntersection (glm::vec3 &source, glm::vec3 &rayDirection, bool rotate) {
	RayTriangleIntersection rti;
	rti.distanceFromCamera = std::numeric_limits<float>::infinity();
	if (rotate) rayDirection = cameraOrientation * rayDirection;
	int index = 0;
	for (ModelTriangle mT : modelTriangles) {
		glm::vec3 e0 = mT.vertices[1] - mT.vertices[0];
		glm::vec3 e1 = mT.vertices[2] - mT.vertices[0];
		glm::vec3 spVector = source - mT.vertices[0];
		glm::mat3 deMatrix(-rayDirection, e0, e1);
		glm::vec3 possibleSol = glm::inverse(deMatrix) * spVector;
		float t = possibleSol[0];
		float u = possibleSol[1];
		float v = possibleSol[2];
		glm::vec3 intersection = mT.vertices[0] + u*(mT.vertices[1] - mT.vertices[0]) + v*(mT.vertices[2]- mT.vertices[0]);

		if ((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0) {
			if (rti.distanceFromCamera > t && t > 0.000001f) {
				rti.distanceFromCamera = t;
				rti.intersectedTriangle = mT;
				rti.intersectionPoint = intersection;
				rti.triangleIndex = index;
				}
		}
		index++;
	}
	return rti;
}

// Checking whether the rti.intersectionPoint is a shadow
bool is_shadow (glm::vec3 lightSource, RayTriangleIntersection rti) {

	glm::vec3 lightDirection = lightSource - rti.intersectionPoint;
	float length = glm::length(lightDirection);

	int index = 0;
	for (ModelTriangle mT : modelTriangles) {
		glm::vec3 e0 = mT.vertices[1] - mT.vertices[0];
		glm::vec3 e1 = mT.vertices[2] - mT.vertices[0];
		glm::vec3 spVector = rti.intersectionPoint - mT.vertices[0];
		glm::mat3 deMatrix(-lightDirection, e0, e1);
		glm::vec3 possibleSol = glm::inverse(deMatrix) * spVector;
		float t = possibleSol[0];
		float u = possibleSol[1];
		float v = possibleSol[2];

		// t > 0.01 to deal with shadow acne problem
		if ((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0) {
			if (t > 0.01 && index != rti.triangleIndex && t < length) {
				return true;
			}
		}
		index++;
	}
	return false;
}

void drawWireframe(DrawingWindow &window) {
	std::vector<std::vector<float>> depthBuffer(WIDTH, std::vector<float>(HEIGHT, 0));
	for (ModelTriangle &t : modelTriangles) {
		CanvasPoint v0 = getCanvasIntersectionPoint(cameraPosition, t.vertices[0], focalLength);
		CanvasPoint v1 = getCanvasIntersectionPoint(cameraPosition, t.vertices[1], focalLength);
		CanvasPoint v2 = getCanvasIntersectionPoint(cameraPosition, t.vertices[2], focalLength);
		CanvasTriangle ct = CanvasTriangle(v0, v1, v2);
		
		unfilledTriangle(window, ct, Colour(255, 255, 255));
	}
	orbit();
}

void renderRasterised(DrawingWindow &window, float focalLength) {

	std::vector<std::vector<float>> depthBuffer(WIDTH, std::vector<float>(HEIGHT, 0));
	for (ModelTriangle &t : modelTriangles) {
		CanvasPoint v0 = getCanvasIntersectionPoint(cameraPosition, t.vertices[0], focalLength);
		CanvasPoint v1 = getCanvasIntersectionPoint(cameraPosition, t.vertices[1], focalLength);
		CanvasPoint v2 = getCanvasIntersectionPoint(cameraPosition, t.vertices[2], focalLength);
		CanvasTriangle ct = CanvasTriangle(v0, v1, v2);
		
		filledTriangle(window, ct, t.colour, depthBuffer);
		// unfilledTriangle(window, ct, Colour(255, 255, 255));
	}
}

void drawRasterised(DrawingWindow &window) {
	renderRasterised(window, focalLength);
	orbit();
}

glm::vec2 flatLighting(glm::vec3 lightSource, RayTriangleIntersection &rti) {
	glm::vec3 lightVector = lightSource - rti.intersectionPoint;

	// Diffuse Lighting
	float distance = glm::length(lightVector);
	float illuminity = 25.0 / (4 * PI * distance * distance);
	float incidence = glm::dot(rti.intersectedTriangle.normal, glm::normalize(lightVector));
	incidence = std::fmax(incidence, 0);
	float diffuse = illuminity * incidence;

	// Specular Lighting
	glm::vec3 viewVector = glm::normalize(cameraPosition - rti.intersectionPoint);
	glm::vec3 normalVector = glm::normalize(rti.intersectedTriangle.normal);
	// r = d - 2(d*n)n where r -> reflectionVector, d -> lightVector, n -> normalVector (normalized)
	glm::vec3 reflectionVector = glm::normalize(glm::normalize(lightVector) - (2 * glm::dot(lightVector, normalVector)) * normalVector);
	float viewAndReflection = glm::dot(-reflectionVector, viewVector);
	viewAndReflection = std::fmax(viewAndReflection, 0);
	float specular = pow(viewAndReflection, 256);

	return glm::vec2(diffuse, specular);
}

glm::vec2 gouraudLighting(glm::vec3 lightSource, RayTriangleIntersection &rti) {
	std::vector<glm::vec3> vN = vertexNormal(rti);
	glm::vec3 weight = barycentric(rti);
	glm::vec3 lightVector = lightSource - rti.intersectionPoint;

	// Diffuse Lighting
	float distance = glm::length(lightVector);
	float illuminity = 25.0 / (4 * PI * distance * distance);
	float incidence = 0;
	for (int i = 0; i < 3; i++) {
		incidence += weight[i] * glm::dot(vN[i], glm::normalize(lightVector));
	}
	incidence = std::fmax(incidence, 0);
	float diffuse = illuminity * incidence;

	// Specular Lighting
	glm::vec3 viewVector = glm::normalize(cameraPosition - rti.intersectionPoint);
	// r = d - 2(d*n)n where r -> reflectionVector, d -> lightVector, n -> normalVector (normalized)
	float specular = 0;
	for (int i = 0; i < 3; i++) {
		glm::vec3 reflectionVector = glm::normalize(glm::normalize(lightVector) - (2 * glm::dot(lightVector, vN[i])) * vN[i]);
		float viewAndReflection = glm::dot(-reflectionVector, viewVector);
		// viewAndReflection = std::fmax(viewAndReflection, 0);
		viewAndReflection = glm::clamp(viewAndReflection, 0.0f, 1.0f);
		specular += weight[i] * pow(viewAndReflection, 256);
	}

	return glm::vec2(diffuse, specular);
}

glm::vec2 phongLighting(glm::vec3 lightSource, RayTriangleIntersection &rti) {
	std::vector<glm::vec3> vN = vertexNormal(rti);
	glm::vec3 weight = barycentric(rti);
	glm::vec3 lightVector = lightSource - rti.intersectionPoint;
	glm::vec3 interpolatedNormal = glm::normalize((weight.x * vN[0]) + (weight.y * vN[1]) + (weight.z * vN[2]));

	// Diffuse Lighting
	float distance = glm::length(lightVector);
	float illuminity = 25.0 / (4 * PI * distance * distance);
	float incidence = glm::dot(interpolatedNormal, glm::normalize(lightVector));
	incidence = glm::clamp(incidence,0.0f,1.0f);
	float diffuse = illuminity * incidence;

	// Specular Lighting
	glm::vec3 viewVector = glm::normalize(cameraPosition - rti.intersectionPoint);
	// r = d - 2(d*n)n where r -> reflectionVector, d -> lightVector, n -> normalVector (normalized)
	glm::vec3 reflectionVector;
	reflectionVector =  glm::normalize(glm::normalize(lightVector) - (2 * glm::dot(lightVector, interpolatedNormal)) * interpolatedNormal);
	float viewAndReflection = glm::dot(-reflectionVector, viewVector);
	// viewAndReflection = std::fmax(viewAndReflection, 0);
	viewAndReflection = glm::clamp(viewAndReflection,0.0f,1.0f);
	float specular = pow(viewAndReflection, 256);

	return glm::vec2(diffuse, specular);
}

glm::vec2 lightingMode(glm::vec3 lightSource, RayTriangleIntersection &rti) {
	if (light_mode == 0) return flatLighting(lightSource, rti);
	else if (light_mode == 1) return gouraudLighting(lightSource, rti);
	else if (light_mode == 2) return phongLighting(lightSource, rti);
	return glm::vec2(0,0);
}

glm::vec2 soft_shadow(RayTriangleIntersection &rti) {
	int count = 0;
	glm::vec2 scale(0,0);
	for (glm::vec3 &l : lights) {
		if (!is_shadow(l, rti)) {
			glm::vec2 singlePoint = lightingMode(l, rti);
			singlePoint[0] = std::fmax(singlePoint[0], 0.0);
			singlePoint[0] = std::fmin(singlePoint[0], 1);
			singlePoint[1] = std::fmax(singlePoint[1], 0.0);
			singlePoint[1] = std::fmin(singlePoint[1], 1);
			scale += singlePoint;
			count++;
		}
	}
	return scale / (float)count;
}

// Returns a vector that has the certian incidence with the normal vector, which lies in the plane created by incident and normal vector
glm::vec3 get_vector(glm::vec3 incident, glm::vec3 normal, float sin_theta) {
	// Ensures that normal vector is same side with the refraction vector
	if (glm::dot(normal, incident) > 0) normal = -normal; 
	float theta = PI / 2 - asin(sin_theta);
	glm::vec3 axis = glm::normalize(glm::cross(incident, normal));
	glm::vec3 projection = glm::normalize((glm::dot(incident, normal) / glm::dot(normal, normal)) * normal);
	if (glm::dot(incident, projection) < 0) projection = -projection;
	// Rodrigues' rotation formula : vrot = v*cos(theta) + (k * v)sin(theta) + k(k.v)(1-cos(theta))
	// Where v = normal, k = axis
	glm::vec3 rotated = normal*cos(theta) + glm::cross(axis, normal)*sin(theta) + axis*glm::dot(axis, normal)*(1 - cos(theta));
	glm::vec3 rotatedOtherSide = normal*cos(-theta) + glm::cross(axis, normal)*sin(-theta) + axis*glm::dot(axis, normal)*(1 - cos(-theta));
	return (glm::dot(rotated, normal) > 0) ? glm::normalize(rotated) : glm::normalize(rotatedOtherSide);
}

// Returns sin value of the given two vectors (positive only)
float get_sin(glm::vec3 &v1, glm::vec3 &v2) {
	glm::vec3 crossed = glm::cross(v1, v2);
	float sin = glm::length(crossed) / (glm::length(v1) * glm::length(v2));
	return std::abs(sin);
}

// void glass(DrawingWindow &window, std::vector<CanvasPoint> &glassPoints, std::vector<RayTriangleIntersection> &glassRti, float air, float glass) {
// 	// sin(air) {angle of incidence} / sin (glass) {angle of refraction} = glass_index / air_index = refractive_ratio 
// 	// => 1.0f * sin(air) = 1.4f * sin(glass)
// 	float refractive_ratio = glass/air;
// 	int index = 0;
// 	for (CanvasPoint &c: glassPoints) {
// 		RayTriangleIntersection glass = glassRti[index];
// 		glm::vec3 incidentVector = glm::normalize(glass.intersectionPoint - cameraPosition);
// 		glm::vec3 glassNormal = glass.intersectedTriangle.normal;
// 		float incidence_sin = get_sin(incidentVector, glassNormal);
// 		float refraction_sin = incidence_sin / refractive_ratio;
// 		glm::vec3 refractedVector = get_vector(incidentVector, glassNormal, refraction_sin);
// 		RayTriangleIntersection rti_glass = getClosestIntersection(glass.intersectionPoint, refractedVector, false);
// 		// if (rti_glass.intersectedTriangle.colour.name == glass.intersectedTriangle.colour.name) {
// 		RayTriangleIntersection	rti_new = getClosestIntersection(rti_glass.intersectionPoint, incidentVector, false);
// 		if (rti_new.distanceFromCamera < std::numeric_limits<float>::infinity()) {
// 			CanvasPoint cp = getCanvasIntersectionPoint(cameraPosition, rti_new.intersectionPoint, focalLength);
// 			uint32_t r = window.getPixelColour(cp.x, cp.y);
// 			window.setPixelColour(c.x, c.y, r);
// 		}
// 		index++;
// 	}
// }

void glass(DrawingWindow &window, std::vector<CanvasPoint> &glassPoints, std::vector<RayTriangleIntersection> &glassRti, float air, float glass) {
	// sin(air) {angle of incidence} / sin (glass) {angle of refraction} = glass_index / air_index = eta
	// eta = n1/n2, c = cos(t1) = -N.V(in)
	// V(ref) = eta * V(in) + (eta*c - sqrt(1-eta^2(1-c^2)))*N
	float eta = glass / air;
	int index = 0;
	for (CanvasPoint &c : glassPoints) {
		RayTriangleIntersection glass = glassRti[index];
		glm::vec3 incidentVector = glm::normalize(glass.intersectionPoint - cameraPosition);
		glm::vec3 glassNormal = glass.intersectedTriangle.normal;
		float cos = abs(glm::dot(glassNormal, incidentVector));
		glm::vec3 refractedVector = eta*incidentVector + (eta*cos - sqrt(1 - eta*eta*(1.0f - cos*cos))) * glassNormal;
		refractedVector = refractedVector;
		RayTriangleIntersection rti_glass = getClosestIntersection(glass.intersectionPoint, refractedVector, false);

		glassNormal = -glassNormal;
		float newCos = abs(glm::dot(glassNormal, refractedVector));
		float newEta = 1 / eta;
		glm::vec3 finalVector = newEta*refractedVector + (newEta*newCos - sqrt(1 - newEta*newEta*(1.0f - newCos*newCos))) * glassNormal;
		finalVector = -finalVector;
		RayTriangleIntersection	rti_new = getClosestIntersection(rti_glass.intersectionPoint, finalVector, false);
		if (rti_new.distanceFromCamera < std::numeric_limits<float>::infinity()) {
			CanvasPoint cp = getCanvasIntersectionPoint(cameraPosition, rti_new.intersectionPoint, focalLength);
			uint32_t r = window.getPixelColour(cp.x, cp.y);
			window.setPixelColour(c.x, c.y, r);
		}
		index++;
	}
}

// MirrorPoints are points that are selected as a mirror and visible from the camera (not hidden)
void mirror(DrawingWindow &window, std::vector<CanvasPoint> &mirrorPoints, std::vector<RayTriangleIntersection> &mirrorRti) {
	int index = 0;
	for (CanvasPoint &c : mirrorPoints) {
		RayTriangleIntersection mirror = mirrorRti[index];
		glm::vec3 normalVector = mirror.intersectedTriangle.normal;
		glm::vec3 viewVector = glm::normalize(cameraPosition - mirror.intersectionPoint);
		viewVector = -viewVector;
		glm::vec3 reflectionVector = glm::normalize(glm::normalize(viewVector) - (2 * glm::dot(viewVector, normalVector)) * normalVector);
		RayTriangleIntersection rti_mirror = getClosestIntersection(mirror.intersectionPoint, reflectionVector, false);
		// uint32_t white = (255 << 24) + (255 << 16) + (255 << 8) + 255;
		// uint32_t black = (255 << 24);
		if (rti_mirror.distanceFromCamera < std::numeric_limits<float>::infinity()) {
			CanvasPoint cp = getCanvasIntersectionPoint(cameraPosition, rti_mirror.intersectionPoint, focalLength);
			uint32_t r = window.getPixelColour(cp.x, cp.y);
			window.setPixelColour(c.x, c.y, r);
		}
		index++;
	}
}

// Finds the closest intersection from the cameraPoisition to every pixel
void drawRayTrace(DrawingWindow &window) {
	std::vector<CanvasPoint> mirrorPoints;
	std::vector<RayTriangleIntersection> mirrorRti;
	std::vector<CanvasPoint> glassPoints;
	std::vector<RayTriangleIntersection> glassRti;
	for (int y = 0; y < window.height; y++) {
		for (int x = 0; x < window.width; x++) {
			glm::vec3 rayDirection = getDirectionVector(CanvasPoint(x,y,0));
			// if(x==0 && y==0) std::cout << rayDirection.x << ", " << rayDirection.y << ", " << rayDirection.z << std::endl;
			RayTriangleIntersection rti = getClosestIntersection(cameraPosition, rayDirection, true);
			if(rti.distanceFromCamera < std::numeric_limits<float>::infinity()) {
				if (mirrors.size() >= 1) {
					bool is_mirror = false;
					for (ModelTriangle t : mirrors) {
						is_mirror = is_mirror || compareVertices(rti.intersectedTriangle.vertices[0], t.vertices[0]);
					}
					if (is_mirror) {
						mirrorPoints.push_back(CanvasPoint(x,y,0));
						mirrorRti.push_back(rti);
					}
				}
				if (glasses.size() >= 1) {
					bool is_glass = false;
					for (ModelTriangle t : glasses) {
						is_glass = is_glass || compareVertices(rti.intersectedTriangle.vertices[0], t.vertices[0]);
					}
					if (is_glass) {
						glassPoints.push_back(CanvasPoint(x,y,0));
						glassRti.push_back(rti);
					}
				}
				Colour c = rti.intersectedTriangle.colour;

				glm::vec2 diffuseAndSpecualar = lightingMode(light, rti);
				float diffuse = diffuseAndSpecualar[0];
				float specular = diffuseAndSpecualar[1];
				// if(x==180 && y==120) std::cout << diffuse << ", " << specular << std::endl;
				float ambiant = 0.2;

				// Diffuse and ambiant contribute to the original object colour
				float illuminity = diffuse + ambiant;
				illuminity = std::fmin(illuminity, 1.0f);
				glm::vec3 objColour(c.red*illuminity, c.green*illuminity, c.blue*illuminity);

				// Specular lighting contribute to the light colour (white)
				specular *= 0.2f;
				glm::vec3 lightColour(255.0f*specular, 255.0f*specular , 255.0f*specular);

				glm::vec3 finalColour = objColour + lightColour;
				// glm::vec3 finalColour = lightColour;

				for (int i = 0; i < 3; i++) finalColour[i] = std::fmin(finalColour[i], 255.0f);

				uint32_t colour = (255 << 24) + (int(finalColour[0]) << 16) + (int(finalColour[1]) << 8) + (int(finalColour[2]));

				if (is_shadow(light, rti)) {
					float scale = 0;
					uint32_t shadow;
					if (shadow_soft) {
						glm::vec2 factor = soft_shadow(rti);
						scale += 0.4f * factor[0] + ambiant;
						factor[1] *= 0.2f;
						uint32_t shadowLight = (255 << 24) + (int(c.red * scale) << 16) + (int(c.green * scale) << 8) + (int(c.blue * scale));
						uint32_t shadowSpecular = (255 << 24) + (int(255 * factor[1]) << 16) + (int(255 * factor[1]) << 8) + (int(255 * factor[1]));
						shadow = shadowLight + shadowSpecular;
					} else {
						scale = ambiant;
						shadow = (255 << 24) + (int(c.red * scale) << 16) + (int(c.green * scale) << 8) + (int(c.blue * scale));
					}
					window.setPixelColour(x, y, shadow);
				}
				else {
					window.setPixelColour(x, y, colour);
				}
			}
		}
	}

	float air_refrective = 1.0f;
	float glass_refrective = 1.0f;
	if (glass_mode) glass(window, glassPoints, glassRti, air_refrective, glass_refrective);

	if (mirror_mode) mirror(window, mirrorPoints, mirrorRti);
	orbit();
}

void draw(DrawingWindow &window) {
	window.clearPixels();

	if (draw_mode == 0) drawWireframe(window);
	else if (draw_mode == 1) drawRasterised(window);
	else if (draw_mode == 2) drawRayTrace(window);
}

std::map<std::string, Colour> parseMtl (const std::string &filename) {
	std::map<std::string, Colour> colourMap;
	std::ifstream inputStream(filename);
	std::string nextLine;

	while(std::getline(inputStream, nextLine)) {
		if (!nextLine.empty()) {
			if (nextLine.at(0) == 'n') {
				std::string c;
				std::vector<std::string> s = split(nextLine, ' ');
				c = s[1];
				if (c == "Cobbles") {
					std::getline(inputStream, nextLine);
					std::vector<std::string> kd = split(nextLine, ' ');
					Colour colour = Colour(c, 255 * std::stof(kd[1]), 255 * std::stof(kd[2]), 255 * std::stof(kd[3]));
					colourMap.insert(std::pair<std::string, Colour>(c, colour));

					std::getline(inputStream, nextLine);
					std::vector<std::string> map_kd = split(nextLine, ' ');
					std::string textureFile = map_kd[1];
					tM = TextureMap(textureFile);
				} else {
					std::getline(inputStream, nextLine);
					std::vector<std::string> kd = split(nextLine, ' ');
					Colour colour = Colour(c, 255 * std::stof(kd[1]), 255 * std::stof(kd[2]), 255 * std::stof(kd[3]));
					colourMap.insert(std::pair<std::string, Colour>(c, colour));
				}
			}
		}
	}
	inputStream.close();
	return colourMap;
}

std::vector<ModelTriangle> parseObj (const std::string &filename, std::map<std::string, Colour> cMap, float sF) {
	std::ifstream inputStream(filename);
	std::string nextLine;
	std::vector<glm::vec3> vertices;
	std::vector<TexturePoint> textures;
	std::vector<ModelTriangle> mT;
	Colour currentColour;

	// Match the index of obj file
	vertices.push_back(glm::vec3(0,0,0));
	textures.push_back(TexturePoint());

	while(std::getline(inputStream, nextLine)) {
		if (!nextLine.empty()) {
			if (nextLine.at(0) == 'u') {
				std::vector<std::string> s = split(nextLine, ' ');
				currentColour = cMap.at(s[1]);
			}
			if (nextLine.at(0) == 'v' && nextLine.at(1) != 't') {
				std::vector<std::string> s = split(nextLine, ' ');
				vertices.push_back(glm::vec3(sF * std::stof(s[1]), sF * std::stof(s[2]), sF * std::stof(s[3])));
			} else if (nextLine.at(0) == 'v' && nextLine.at(1) != 't') {
				std::vector<std::string> vt = split(nextLine, ' ');
				TexturePoint tp = TexturePoint(std::stof(vt[1]), std::stof(vt[2]));
				textures.push_back(tp);
			} else if (nextLine.at(0) == 'f') {
				if (currentColour.name == "Cobbles") {
					std::vector<std::string> s = split(nextLine, ' ');
					std::vector<int> facet;
					std::vector<int> tFacet;
					for (int i = 0; i < 3; i++) {
						facet.push_back(std::stoi(split(s[i + 1], '/').at(0)));
						tFacet.push_back(std::stoi(split(s[i + 1], '/').at(1)));
					}
					ModelTriangle m = ModelTriangle(vertices[facet[0]], vertices[facet[1]], vertices[facet[2]], currentColour);
					m.texturePoints = {textures[tFacet[0]], textures[tFacet[1]], textures[tFacet[2]]};
					m.normal = glm::normalize(glm::cross(m.vertices[1] - m.vertices[0], m.vertices[2] - m.vertices[0]));
					mT.push_back(m);
					std::cout << m << std::endl;
				} else {
					std::vector<std::string> s = split(nextLine, ' ');
					std::vector<int> facet;
					for (int i = 0; i < 3; i++) facet.push_back(std::stoi(split(s[i + 1], '/').at(0)));
					ModelTriangle m = ModelTriangle(vertices[facet[0]], vertices[facet[1]], vertices[facet[2]], currentColour);
					m.normal = glm::normalize(glm::cross(m.vertices[1] - m.vertices[0], m.vertices[2] - m.vertices[0]));
					mT.push_back(m);
					// if (facet[0] == 62) mirrors.push_back(m);
					if (currentColour.name == "Magenta") mirrors.push_back(m);
					if (currentColour.name == "Red") glasses.push_back(m);
				}
			} else {
				continue;
			}
		}
	}
	inputStream.close();
	return mT;
}

void parseFiles(const std::string &objFile, const std::string &mtlFile, float scalingFactor) {
	colourMap = parseMtl(mtlFile);
	modelTriangles = parseObj(objFile, colourMap, scalingFactor);
}


void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_z) rotate_camera(false, -PI / 16);
		else if (event.key.keysym.sym == SDLK_x) rotate_camera(false, PI / 16);
		else if (event.key.keysym.sym == SDLK_c) rotate_camera(true, -PI / 16);
		else if (event.key.keysym.sym == SDLK_v) rotate_camera(true, PI / 16);
		else if (event.key.keysym.sym == SDLK_w) lightMove? translate_light(1, true): translate_camera(1, true);
		else if (event.key.keysym.sym == SDLK_a) lightMove? translate_light(0, false): translate_camera(0, false);
		else if (event.key.keysym.sym == SDLK_s) lightMove? translate_light(1, false): translate_camera(1, false);
		else if (event.key.keysym.sym == SDLK_d) lightMove? translate_light(0, true): translate_camera(0, true);
		else if (event.key.keysym.sym == SDLK_q) lightMove? translate_light(2, true): translate_camera(2, true);
		else if (event.key.keysym.sym == SDLK_e) lightMove? translate_light(2, false): translate_camera(2, false);
		else if (event.key.keysym.sym == SDLK_j) {
			if (lightMove) change_lightSize(false);
			std::cout << "Number of lights: " << lights.size() << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_k) {
			if (lightMove) change_lightSize(true);
			std::cout << "Number of lights: " << lights.size() << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_LEFT) change_orientation(false, PI / 16);
		else if (event.key.keysym.sym == SDLK_RIGHT) change_orientation(false, -PI / 16);
		else if (event.key.keysym.sym == SDLK_UP) change_orientation(true, PI / 16);
		else if (event.key.keysym.sym == SDLK_DOWN) change_orientation(true, -PI / 16);
		else if (event.key.keysym.sym == SDLK_1) {
			draw_mode = 0;
			std::cout << "Changed to Wireframe drawing" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_2) {
			draw_mode = 1;
			std::cout << "Changed to Rasterised drawing" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_3) {
			draw_mode = 2;
			std::cout << "Changed to RayTrace drawing" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_4) {
			light_mode = 0;
			std::cout << "Flat lighting applied" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_5) {
			light_mode = 1;
			std::cout << "Gouraud lighting applied" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_6) {
			light_mode = 2;
			std::cout << "Phong lighting applied" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_7) {
			if (shadow_soft) std::cout << "Soft shadow off" << std::endl;
			else std::cout << "Soft shadow on" << std::endl;
			shadow_soft = !shadow_soft;
		}
		else if (event.key.keysym.sym == SDLK_9) {
			if (glass_mode) std::cout << "Glass off" << std::endl;
			else std::cout << "Glass on" << std::endl;
			glass_mode = !glass_mode;
		}
		else if (event.key.keysym.sym == SDLK_0) {
			if (mirror_mode) std::cout << "Mirror off" << std::endl;
			else std::cout << "Mirror on" << std::endl;
			mirror_mode = !mirror_mode;
		}

		else if (event.key.keysym.sym == SDLK_o) orbits = !orbits;
		else if (event.key.keysym.sym == SDLK_p) {
			lightMove = !lightMove;
			if (lightMove) {
				std::cout << "Translation now changes Light Position" << std::endl;
			} else {
				std::cout << "Translation now changes Camera Position" << std::endl;
			}
		}
		else if (event.key.keysym.sym == SDLK_l) lookAt();
		else if (event.key.keysym.sym == SDLK_u) unfilledTriangle(window);
		else if (event.key.keysym.sym == SDLK_f) filledTriangle(window);
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	// std::string objFile = "logo.obj";
	// std::string objFile = "sphere.obj";
	std::string objFile = "cornell-box.obj";
	std::string mtlFile = "cornell-box.mtl";
	std::string textureFile = "texture.ppm";
	parseFiles(objFile, mtlFile, 0.35);

	lights.push_back(light);
	initialize_lights();
	// draw(window);
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
