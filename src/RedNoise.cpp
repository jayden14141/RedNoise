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
std::vector<ModelTriangle> modelTriangles;
glm::vec3 cameraPosition = glm::vec3(0.0, 0.0, 6.0);
float focalLength = 2.0;
glm::mat3 cameraOrientation = glm::mat3(1.0);
bool orbits = false;
bool lightMove = false;
int draw_mode = 2;
int light_mode = 2;
glm::vec3 light(0, 0, 0.7);

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

// Translates either camera poisition or light position depending on 'lightMove' variable
void translate(glm::vec3 &what, int where, bool positive) {

	//translate by x
	if (where == 0) {
		if (positive) {
			what += glm::vec3(0.2, 0, 0);
		} else {
			what -= glm::vec3(0.2, 0, 0);
		}
	//translate by y
	} else if (where == 1) {
		if (positive) {
			what += glm::vec3(0, 0.2, 0);
		} else {
			what -= glm::vec3(0, 0.2, 0);
		}
	//translate by z
	} else {
		if (positive) {
			what += glm::vec3(0, 0, 0.2);
		} else {
			what -= glm::vec3(0, 0, 0.2);
		}
	}
	std::cout << light.x << ", " << light.y << ", " << light.z << std::endl;
}

void lookAt() {
	glm::vec3 forward = glm::normalize(cameraPosition);
	glm::vec3 right = glm::normalize(glm::cross(glm::vec3(0,1,0), forward));
	glm::vec3 up = glm::normalize(glm::cross(forward, right));

	cameraOrientation[0] = right;
	cameraOrientation[1] = up;
	cameraOrientation[2] = forward;
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

void textureMapping(const std::string &filename, DrawingWindow &window, CanvasTriangle texture, CanvasTriangle canvas) {
	TextureMap tM = TextureMap(filename);
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

	// unfilledTriangle(window, canvas, Colour(255, 255, 255));
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
		glm::normalize(v);
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

RayTriangleIntersection getClosestIntersection (glm::vec3 &rayDirection) {
	RayTriangleIntersection rti;
	rti.distanceFromCamera = std::numeric_limits<float>::infinity();
	rayDirection = cameraOrientation * rayDirection;
	
	int index = 0;
	for (ModelTriangle mT : modelTriangles) {
		glm::vec3 e0 = mT.vertices[1] - mT.vertices[0];
		glm::vec3 e1 = mT.vertices[2] - mT.vertices[0];
		glm::vec3 spVector = cameraPosition - mT.vertices[0];
		glm::mat3 deMatrix(-rayDirection, e0, e1);
		glm::vec3 possibleSol = glm::inverse(deMatrix) * spVector;
		float t = possibleSol[0];
		float u = possibleSol[1];
		float v = possibleSol[2];
		glm::vec3 intersection = mT.vertices[0] + u*(mT.vertices[1] - mT.vertices[0]) + v*(mT.vertices[2]- mT.vertices[0]);


		if ((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0) {
			if (rti.distanceFromCamera > t && t >= 0) {
				// std::string a("White");
				// if (a.compare(mT.colour.name) == 0) {
				// 	for (int i = 0; i < 3; i++) {
				// 			std::cout << rayDirection.x << ", " << rayDirection.y << ", " << rayDirection.z << std::endl;
				// 	}
				// 	// std::cout <<  << std::endl;
				// }
				rti.distanceFromCamera = t;
				rti.intersectedTriangle = mT;
				rti.intersectionPoint = intersection;
				rti.triangleIndex = index;
				// std::cout << rti << std::endl;
			}
		}
		index++;
	}
	return rti;
}

// Checking whether the rti.intersectionPoint is a shadow
bool is_shadow (RayTriangleIntersection rti) {

	glm::vec3 lightDirection = light - rti.intersectionPoint;
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
	for (ModelTriangle t : modelTriangles) {
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
	for (ModelTriangle t : modelTriangles) {
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

float flatLighting(RayTriangleIntersection &rti) {
	glm::vec3 lightVector = light - rti.intersectionPoint;

	// Diffuse Lighting
	float distance = glm::length(lightVector);
	float illuminity = 40.0 / (4 * PI * distance * distance);
	float incidence = glm::dot(rti.intersectedTriangle.normal, glm::normalize(lightVector));
	incidence = std::fmax(incidence, 0);
	float diffuse = illuminity * incidence;

	// Specular Lighting
	glm::vec3 viewVector = glm::normalize(cameraPosition - rti.intersectionPoint);
	glm::vec3 normalVector = glm::normalize(rti.intersectedTriangle.normal);
	// r = d - 2(d*n)n where r -> reflectionVector, d -> lightVector, n -> normalVector (normalized)
	glm::vec3 reflectionVector = glm::normalize(glm::normalize(lightVector) - (2 * glm::dot(lightVector, normalVector)) * normalVector);
	float viewAndReflection = dot(reflectionVector, viewVector);
	viewAndReflection = std::fmax(viewAndReflection, 0);
	float specular = pow(viewAndReflection, 256);

	float brightness = diffuse + specular;

	return brightness;
}

float gouraudLighting(RayTriangleIntersection &rti) {
	std::vector<glm::vec3> vN = vertexNormal(rti);
	glm::vec3 weight = barycentric(rti);
	glm::vec3 lightVector = light - rti.intersectionPoint;

	// Diffuse Lighting
	float distance = glm::length(lightVector);
	float illuminity = 40.0 / (4 * PI * distance * distance);
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
		float viewAndReflection = glm::dot(reflectionVector, viewVector);
		viewAndReflection = std::fmax(viewAndReflection, 0);
		specular += weight[i] * pow(viewAndReflection, 256);
	}

	float brightness = diffuse + specular;
	return brightness;
}

float phongLighting(RayTriangleIntersection &rti) {
	std::vector<glm::vec3> vN = vertexNormal(rti);
	glm::vec3 weight = barycentric(rti);
	glm::vec3 lightVector = light - rti.intersectionPoint;
	glm::vec3 interpolatedNormal = glm::normalize((weight.x * vN[0]) + (weight.y * vN[1]) + (weight.z * vN[2]));

	// Diffuse Lighting
	float distance = glm::length(lightVector);
	float illuminity = 40.0 / (4 * PI * distance * distance);
	float incidence = glm::dot(interpolatedNormal, glm::normalize(lightVector));
	incidence = std::fmax(incidence, 0);
	float diffuse = illuminity * incidence;

	// Specular Lighting
	glm::vec3 viewVector = glm::normalize(cameraPosition - rti.intersectionPoint);
	// r = d - 2(d*n)n where r -> reflectionVector, d -> lightVector, n -> normalVector (normalized)
	glm::vec3 reflectionVector;
	reflectionVector =  glm::normalize(glm::normalize(lightVector) - (2 * glm::dot(lightVector, interpolatedNormal)) * interpolatedNormal);
	float viewAndReflection = glm::dot(reflectionVector, viewVector);
	viewAndReflection = std::fmax(viewAndReflection, 0);
	float specular = pow(viewAndReflection, 256);

	float brightness = diffuse + specular;
	return brightness;
}

float lightingMode(RayTriangleIntersection &rti) {
	if (light_mode == 0) return flatLighting(rti);
	else if (light_mode == 1) return gouraudLighting(rti);
	else if (light_mode == 2) return phongLighting(rti);
	return 0.0f;
}

// Finds the closest intersection from the cameraPoisition to every pixel
void drawRayTrace(DrawingWindow &window) {
	for (int y = 0; y < window.height; y++) {
		for (int x = 0; x < window.width; x++) {
			glm::vec3 rayDirection = getDirectionVector(CanvasPoint(x,y,0));
			RayTriangleIntersection rti = getClosestIntersection(rayDirection);
			if(rti.distanceFromCamera < std::numeric_limits<float>::infinity()) {
				Colour c = rti.intersectedTriangle.colour;

				float brightness = lightingMode(rti);

				float ambiant = 0.3;
				brightness = std::fmax(brightness, ambiant);
				brightness = std::fmin(brightness, 1);

				uint32_t colour = (255 << 24) + (int(c.red * brightness) << 16) + (int(c.green * brightness) << 8) + (int(c.blue * brightness));
				if (is_shadow(rti)) {
					uint32_t shadow = (255 << 24) + (int(c.red * ambiant) << 16) + (int(c.green * ambiant) << 8) + (int(c.blue * ambiant));
					window.setPixelColour(x, y, shadow);
				}
				else window.setPixelColour(x, y, colour);
			}
		}
	}
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

				std::getline(inputStream, nextLine);
				std::vector<std::string> kd = split(nextLine, ' ');
				Colour colour = Colour(c, 255 * std::stof(kd[1]), 255 * std::stof(kd[2]), 255 * std::stof(kd[3]));
				colourMap.insert(std::pair<std::string, Colour>(c, colour));
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
	std::vector<std::vector<int>> facets;
	std::vector<ModelTriangle> mT;
	Colour currentColour;

	// Match the index of obj file
	vertices.push_back(glm::vec3(0,0,0));

	while(std::getline(inputStream, nextLine)) {
		if (!nextLine.empty()) {
			if (nextLine.at(0) == 'u') {
				std::vector<std::string> s = split(nextLine, ' ');
				currentColour = cMap.at(s[1]);
			}
			if (nextLine.at(0) == 'v') {
				std::vector<std::string> s = split(nextLine, ' ');
				vertices.push_back(glm::vec3(sF * std::stof(s[1]), sF * std::stof(s[2]), sF * std::stof(s[3])));
			} else if (nextLine.at(0) == 'f') {
				std::vector<std::string> s = split(nextLine, ' ');
				std::vector<int> facet;
				facet.push_back(std::stoi(split(s[1], '/').at(0)));
				facet.push_back(std::stoi(split(s[2], '/').at(0)));
				facet.push_back(std::stoi(split(s[3], '/').at(0)));
				// facets.push_back(facet);
				ModelTriangle m = ModelTriangle(vertices[facet[0]], vertices[facet[1]], vertices[facet[2]], currentColour);
				m.normal = glm::normalize(glm::cross(m.vertices[1] - m.vertices[0], m.vertices[2] - m.vertices[0]));
				mT.push_back(m);
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
		else if (event.key.keysym.sym == SDLK_w) lightMove? translate(light, 1, true): translate(cameraPosition, 1, true);
		else if (event.key.keysym.sym == SDLK_a) lightMove? translate(light, 0, false): translate(cameraPosition, 0, false);
		else if (event.key.keysym.sym == SDLK_s) lightMove? translate(light, 1, false): translate(cameraPosition, 1, false);
		else if (event.key.keysym.sym == SDLK_d) lightMove? translate(light, 0, true): translate(cameraPosition, 0, true);
		else if (event.key.keysym.sym == SDLK_q) lightMove? translate(light, 2, true): translate(cameraPosition, 2, true);
		else if (event.key.keysym.sym == SDLK_e) lightMove? translate(light, 2, false): translate(cameraPosition, 2, false);
		else if (event.key.keysym.sym == SDLK_LEFT) change_orientation(false, -PI / 16);
		else if (event.key.keysym.sym == SDLK_RIGHT) change_orientation(false, PI / 16);
		else if (event.key.keysym.sym == SDLK_UP) change_orientation(true, -PI / 16);
		else if (event.key.keysym.sym == SDLK_DOWN) change_orientation(true, PI / 16);
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
	std::string objFile = "cornell-box.obj";
	// std::string objFile = "sphere.obj";
	std::string mtlFile = "cornell-box.mtl";
	parseFiles(objFile, mtlFile, 0.35);
	// draw(window);
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
