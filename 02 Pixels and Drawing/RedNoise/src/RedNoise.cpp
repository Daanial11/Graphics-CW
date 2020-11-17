#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>

#define WIDTH 320
#define HEIGHT 240

std::vector<float> interpolateSingleFloats(float firstNo, float secondNo, int listSize);
std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues);

void draw(DrawingWindow &window) {
	window.clearPixels();

	glm::vec3 topLeft(255, 0, 0);        // red 
	glm::vec3 topRight(0, 0, 255);       // blue 
	glm::vec3 bottomRight(0, 255, 0);    // green 
	glm::vec3 bottomLeft(255, 255, 0);   // yellow

	std::vector<glm::vec3> leftSide = interpolateThreeElementValues(topLeft, bottomLeft, 240);
	std::vector<glm::vec3> rightSide = interpolateThreeElementValues(topRight, bottomRight, 240);

	
	for (size_t y = 0; y < window.height; y++) {
		std::vector<glm::vec3> currentRow = interpolateThreeElementValues(leftSide[y], rightSide[y], WIDTH);
		for (size_t x = 0; x < window.width; x++) {
			float red = (currentRow[x])[0];
			float green = (currentRow[x])[1];
			float blue = (currentRow[x])[2];
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}

void update(DrawingWindow &window) {
	// Function for performing animation (shifting artifacts or moving the camera)
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
	} else if (event.type == SDL_MOUSEBUTTONDOWN) window.savePPM("output.ppm");
}

std::vector<float> interpolateSingleFloats(float firstNo, float secondNo, int listSize) {
	std::vector<float> result;
	float jumps = (secondNo - firstNo)/(listSize-1);
	
	for (int x = 0; x<listSize; x++) {
		float number = firstNo + (x*jumps);
		result.push_back(number);
	}

	return result;

}
std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues) {
	std::vector<glm::vec3> result;

	std::vector<float> firstRow = interpolateSingleFloats(from[0], to[0], numberOfValues);
	std::vector<float> secondRow = interpolateSingleFloats(from[1], to[1], numberOfValues);
	std::vector<float> thirdRow = interpolateSingleFloats(from[2], to[2], numberOfValues);

	

	for (int x = 0; x < numberOfValues; x++) {
		glm::vec3 vector(firstRow[x], secondRow[x], thirdRow[x]);
		result.push_back(vector);
	}

	return result;

}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

	
		

	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		update(window);
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}

