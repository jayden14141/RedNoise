#include <algorithm>
#include <CanvasTriangle.h>
#include <CanvasPoint.h>
#include <Colour.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>

#define WIDTH 320
#define HEIGHT 240

bool compare(CanvasPoint p1, CanvasPoint p2) {
	return p1.y < p2.y;
}

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
	float numberOfSteps = std::max(abs(xDiff), abs (yDiff));
	float xStepSize = xDiff / numberOfSteps;
	float yStepSize = yDiff / numberOfSteps;
	uint32_t colour = (255 << 24) + (c.red << 16) + (c.green << 8) + c.blue;

	for (float i = 0.0; i <= numberOfSteps; i++) {
		float x = from.x + (xStepSize * i);
		float y = from.y + (yStepSize * i);
		window.setPixelColour(round(x), round(y), colour);
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

void unfilledTriangle(DrawingWindow &window, CanvasPoint v0, CanvasPoint v1, CanvasPoint v2, Colour colour) {
	CanvasTriangle t = CanvasTriangle(v0, v1, v2);

	line(window, t.v0(), t.v1(), colour);
	line(window, t.v1(), t.v2(), colour);
	line(window, t.v2(), t.v0(), colour);
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

void filledTriangle(DrawingWindow &window) {
	CanvasPoint v0 = CanvasPoint(rand() % window.width - 1, rand() % window.height - 1);
	CanvasPoint v1 = CanvasPoint(rand() % window.width - 1, rand() % window.height - 1);
	CanvasPoint v2 = CanvasPoint(rand() % window.width - 1, rand() % window.height - 1);
	// CanvasTriangle t = CanvasTriangle(v0, v1, v2);

	Colour c = Colour(rand() % 256, rand() % 256, rand() % 256);
	std::vector<CanvasPoint> v;
	v.push_back(v0);
	v.push_back(v1);
	v.push_back(v2);
	sort(v.begin(), v.end(), compare);

	// yDiff : xDiff = y1 - y0: alpha - x0
	// alpha = (xDiff * (y1-y0) / yDiff) + x0
	float xDiff = v[2].x - v[0].x;
	float yDiff = v[2].y - v[0].y;
	float alpha = (xDiff * (v[1].y-v[0].y) / yDiff) + v[0].x;
	CanvasPoint left, right;

	if (v[1].x < alpha) {
		left = v[1];
		right = CanvasPoint(alpha, v[1].y);
	} else {
		left = CanvasPoint(alpha, v[1].y);
		right = v[1];
	}

	// line(window, v1_a, v1_b, Colour(255, 255, 255));
	fill(window, true, left, right, v[0], c);
	fill(window, false, left, right, v[2], c);

	unfilledTriangle(window, v0, v1, v2, Colour(255, 255, 255));
}

void draw(DrawingWindow &window) {
	window.clearPixels();
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
	// line(window, CanvasPoint(0,0), CanvasPoint(window.width/2, window.height/2));
	// line(window, CanvasPoint(window.width - 1, 0), CanvasPoint(window.width/2, window.height/2));
	// line(window, CanvasPoint(window.width/2, 0), CanvasPoint(window.width/2, window.height - 1));
	// line(window, CanvasPoint(window.width/3, window.height/2), CanvasPoint(window.width*2/3, window.height/2));

}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
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
		// draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
