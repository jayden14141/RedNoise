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
int draw_mode = 2;

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

// Using plain float
// void line(DrawingWindow &window, float fromX, float fromY, float toX, float toY) {
// 	float xDiff = toX - fromX;
// 	float yDiff = toY - fromY;
// 	float numberOfSteps = std::max(abs(xDiff), abs (yDiff));
// 	float xStepSize = xDiff / numberOfSteps;
// 	float yStepSize = yDiff / numberOfSteps;

// 	for (float i = 0.0; i <= numberOfSteps; i++) {
// 		float x = fromX + (xStepSize * i);
// 		float y = fromY + (yStepSize * i);
// 		window.setPixelColour(round(x), round(y), Colour(255, 255, 255));
// 	}
// }

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

std::vector<uint32_t> textureColour (TextureMap tm, CanvasPoint from, CanvasPoint to, int numberOfValues) {
	std::vector<CanvasPoint> texturePoint = interpolateCanvasPoint(from, to, numberOfValues);
	std::vector<uint32_t> colour;
	for (int i = 0; i < numberOfValues; i++) {
		colour.push_back(tm.pixels[round(texturePoint[i].x) + round(texturePoint[i].y) * tm.width]);
	}
	return colour;
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
				mT.push_back(ModelTriangle(vertices[facet[0]], vertices[facet[1]], vertices[facet[2]], currentColour));
			} else {
				continue;
			}
		}
	}
	inputStream.close();

	// for (std::vector<int> f : facets) {
	// 	mT.push_back(ModelTriangle(vertices[f[0]], vertices[f[1]], vertices[f[2]], Colour()));
	// }
	return mT;
}

CanvasPoint getCanvasIntersectionPoint(glm::vec3 &cameraPosition, glm::vec3 &vertexPosition, float focalLength) {

	glm::vec3 direction = vertexPosition - cameraPosition;
	direction = direction * cameraOrientation;

	float canvasX = HEIGHT * (focalLength * -direction[0] / direction[2]) + WIDTH / 2;
	float canvasY = HEIGHT * (focalLength * direction[1] / direction[2]) + HEIGHT / 2;
	float depth =  -direction[2];
	// std::cout << depth << std::endl;
	return CanvasPoint(canvasX, canvasY, depth);
}

// Returns a directional vector from the 2D canvaspoint to the camera Position
glm::vec3 getDirectionVector(CanvasPoint canvasPosition) {
	glm::vec3 direction = glm::vec3(canvasPosition.x, canvasPosition.y, 0) + cameraPosition;
	float dZ = -1 * direction[2];
	float dX = -dZ * (direction[0] - WIDTH / 2) / (HEIGHT * focalLength);
	float dY = dZ * (direction[1] - HEIGHT / 2) / (HEIGHT * focalLength);
	return glm::normalize(glm::vec3(dX, dY, dZ));
}

RayTriangleIntersection getClosestIntersection (glm::vec3 rayDirection) {
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
		glm::vec3 intersection = mT.vertices[0] + u * (mT.vertices[1] - mT.vertices[0]) + v*(mT.vertices[2]- mT.vertices[0]);


		if ((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0) {
			if (rti.distanceFromCamera > t && t >= 0) {
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

bool is_shadow (RayTriangleIntersection rti, glm::vec3 rayDirection) {

	glm::vec3 direction = rayDirection - rti.intersectionPoint;
	int index = 0;
	for (ModelTriangle mT : modelTriangles) {
		glm::vec3 e0 = mT.vertices[1] - mT.vertices[0];
		glm::vec3 e1 = mT.vertices[2] - mT.vertices[0];
		glm::vec3 spVector = rti.intersectionPoint - mT.vertices[0];
		glm::mat3 deMatrix(-direction, e0, e1);
		glm::vec3 possibleSol = glm::inverse(deMatrix) * spVector;
		float t = possibleSol[0];
		float u = possibleSol[1];
		float v = possibleSol[2];

		if ((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0) {
			if (t >= 0 && index != rti.triangleIndex) {
				return true;
			}
		}
		index++;
	}
	return false;
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

void parseFiles(const std::string &objFile, const std::string &mtlFile, float scalingFactor) {
	colourMap = parseMtl(mtlFile);
	modelTriangles = parseObj(objFile, colourMap, scalingFactor);
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

void drawRayTrace(DrawingWindow &window) {
	for (int y = 0; y < window.height; y++) {
		for (int x = 0; x < window.width; x++) {
			RayTriangleIntersection rti = getClosestIntersection(getDirectionVector(CanvasPoint(x,y,0)));
			if(!isinf(rti.distanceFromCamera)) {
				Colour c = rti.intersectedTriangle.colour;
				uint32_t colour = (255 << 24) + (c.red << 16) + (c.green << 8) + c.blue;
				if (is_shadow(rti, glm::vec3(0,0,1))) {
					uint32_t shadow = (255 << 24);
					window.setPixelColour(x, y, shadow);
				}
				else window.setPixelColour(x, y, colour);
				// window.setPixelColour(x, y, colour);
				// Checking the light location (150, 45)
				// if (colour == (255 << 24) + (255 << 16) + (255 << 8) + 255) std::cout << x << ", " << y << std::endl;
			}
		}
	}
	orbit();
}

void drawRasterised(DrawingWindow &window) {
	// for (size_t y = 0; y < window.height; y++) {
	// 	for (size_t x = 0; x < window.width; x++) {
	// 		float red = rand() % 256;
	// 		float green = 0.0;
	// 		float blue = 0.0;
	// 		uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
	// 		window.setPixelColour(x, y, colour);
	// 	}
	// }

	// Wk2 Task03
	// RGB (0,0,0) => Black || RGB (255, 255, 255) => White
	// (0,0) => Top Left / (255, 255) => Bottom Right

	// std::vector<float> v = interpolateSingleFloats(255, 0, window.width);
	// for (size_t y = 0; y < window.height; y++) {
	// 	float black = 0.0;
	// 	for (size_t x = 0; x < window.width; x++) {
	// 		black = v.at(x);
	// 		uint32_t colour = (255 << 24) + (int(black) << 16) + (int(black) << 8) + int(black);
	// 		window.setPixelColour(x, y, colour);
	// 	}
	// }

	// Wk2 Task05
	// glm::vec3 topLeft(255, 0, 0);        // red 
	// glm::vec3 topRight(0, 0, 255);       // blue 
	// glm::vec3 bottomRight(0, 255, 0);    // green 
	// glm::vec3 bottomLeft(255, 255, 0);   // yellow

	// std::vector<glm::vec3> left = interpolateThreeElementValues(topLeft, bottomLeft, window.width);
	// std::vector<glm::vec3> right = interpolateThreeElementValues(topRight, bottomRight, window.width);
	// for (size_t y = 0; y < window.height; y++) {
	// 	std::vector<glm::vec3> col = interpolateThreeElementValues(left[y], right[y], window.width);
	// 	for (size_t x = 0; x < window.width; x++) {
	// 		uint32_t colour = (255 << 24) + (int(col[x][0]) << 16) + (int(col[x][1]) << 8) + int(col[x][2]);
	// 		window.setPixelColour(x, y, colour);
	// 	}
	// }


	// Wk3 Task02
	// Colour c = Colour(255, 255, 255);
	// line(window, CanvasPoint(0,0), CanvasPoint(window.width/2, window.height/2), c);
	// line(window, CanvasPoint(window.width - 1, 0), CanvasPoint(window.width/2, window.height/2), c);
	// line(window, CanvasPoint(window.width/2, 0), CanvasPoint(window.width/2, window.height - 1), c);
	// line(window, CanvasPoint(window.width/3, window.height/2), CanvasPoint(window.width*2/3, window.height/2), c);


	// Wk3 Task06
	// CanvasPoint t0 = CanvasPoint(195, 5);
	// CanvasPoint t1 = CanvasPoint(395, 380);
	// CanvasPoint t2 = CanvasPoint(65, 330);
	// CanvasTriangle t = CanvasTriangle(t0, t1, t2);

	// CanvasPoint c0 = CanvasPoint(160, 10);
	// CanvasPoint c1 = CanvasPoint(300, 230);
	// CanvasPoint c2 = CanvasPoint(10, 150);
	// CanvasTriangle c = CanvasTriangle(c0, c1, c2);

	// std::string filename = "texture.ppm";
	// textureMapping(filename, window, t, c);


	// Wk4 Task02, 03
	renderRasterised(window, focalLength);
	orbit();
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

void draw(DrawingWindow &window) {
	window.clearPixels();
	if (draw_mode == 0) drawWireframe(window);
	else if (draw_mode == 1) drawRasterised(window);
	else if (draw_mode == 2) drawRayTrace(window);
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_z) rotate_camera(false, -PI / 16);
		else if (event.key.keysym.sym == SDLK_x) rotate_camera(false, PI / 16);
		else if (event.key.keysym.sym == SDLK_c) rotate_camera(true, -PI / 16);
		else if (event.key.keysym.sym == SDLK_v) rotate_camera(true, PI / 16);
		else if (event.key.keysym.sym == SDLK_w) translate_camera(1, true);
		else if (event.key.keysym.sym == SDLK_a) translate_camera(0, false);
		else if (event.key.keysym.sym == SDLK_s) translate_camera(1, false);
		else if (event.key.keysym.sym == SDLK_d) translate_camera(0, true);
		else if (event.key.keysym.sym == SDLK_q) translate_camera(2, true);
		else if (event.key.keysym.sym == SDLK_e) translate_camera(2, false);
		else if (event.key.keysym.sym == SDLK_LEFT) change_orientation(false, -PI / 16);
		else if (event.key.keysym.sym == SDLK_RIGHT) change_orientation(false, PI / 16);
		else if (event.key.keysym.sym == SDLK_UP) change_orientation(true, -PI / 16);
		else if (event.key.keysym.sym == SDLK_DOWN) change_orientation(true, PI / 16);
		else if (event.key.keysym.sym == SDLK_1) draw_mode = 0;
		else if (event.key.keysym.sym == SDLK_2) draw_mode = 1;
		else if (event.key.keysym.sym == SDLK_3) draw_mode = 2;
		else if (event.key.keysym.sym == SDLK_o) orbits = !orbits;
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
	std::string mtlFile = "cornell-box.mtl";
	parseFiles(objFile, mtlFile, 0.35);
	
	// draw(window);
	// std::vector<float> result;
	// result = interpolateSingleFloats(2.2, 8.5, 7);
	// for(size_t i=0; i<result.size(); i++) std::cout << result[i] << " ";

	// std::vector<glm::vec3> result;
	// glm::vec3 from(1.0, 4.0, 9.2);
	// glm::vec3 to(4.0, 1.0, 9.8);
	// result = interpolateThreeElementValues(from, to, 4);
	// for(size_t i=0; i<result.size(); i++) {
	// 	std::cout << "(" << result[i][0] << "," << result[i][1] << "," << result[i][2] << ")" << " ";
	// 	std::cout << std::endl;
	// }
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		draw(window);
		// drawRasterised(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
