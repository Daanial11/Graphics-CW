#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <CanvasPoint.h>
#include <Colour.h>
#include <CanvasTriangle.h>
#include <TextureMap.h>

#define WIDTH 320
#define HEIGHT 240

std::vector<float> interpolateSingleFloats(float firstNo, float secondNo, int listSize);
std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues);

void drawline(CanvasPoint startPoint, CanvasPoint endPoint, uint32_t colour, DrawingWindow &window){
	 float xDiff = endPoint.x - startPoint.x;
	 float yDiff = endPoint.y - startPoint.y;

	 float steps = std :: max(abs(xDiff), abs(yDiff));

	 float xStep = xDiff/steps;
	 float yStep = yDiff/steps;

	 for (float i=0.0; i<steps; i++) {
		 float x = startPoint.x + (xStep*i);
		 float y = startPoint.y + (yStep*i);

		 window.setPixelColour(round(x), round(y), colour);
	 }
}

void drawTriangle(DrawingWindow &window){

	float red = rand()%256;
	float green = rand()%256;
	float blue = rand()%256;

	uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);

	std::cout << "got here" << std::endl;
	

	CanvasTriangle triangle(CanvasPoint(rand()%320 + 1, rand()%240 + 1), CanvasPoint(rand()%320 + 1, rand()%240 + 1), CanvasPoint(rand()%320 + 1, rand()%240 + 1));
	
	drawline(triangle.v0(), triangle.v1(), colour, window);
	drawline(triangle.v1(), triangle.v2(), colour, window);
	drawline(triangle.v2(), triangle.v0(), colour, window);

	//window.renderFrame();

}

std::vector<glm::vec2> interpolateLineCords(CanvasPoint startPoint, CanvasPoint endPoint){

	std::vector<glm::vec2> result;

	float xDiff = endPoint.x - startPoint.x;
	float yDiff = endPoint.y - startPoint.y;

	float steps = std :: max(abs(xDiff), abs(yDiff));

	float xStep = xDiff/steps;
	float yStep = yDiff/steps;

	for (float i=0.0; i<steps; i++) {
		float x = startPoint.x + (xStep*i);
		float y = startPoint.y + (yStep*i);
		glm::vec2 currentCord(round(x), round(y));

		result.push_back(currentCord);
	}

	return result;

}

bool insideTriangle(CanvasTriangle &triangle, const CanvasPoint &point) {

	float alpha = ((triangle.v1().y - triangle.v2().y)*(point.x - triangle.v2().x) + (triangle.v2().x - triangle.v1().x)*(point.y - triangle.v2().y)) / 
	((triangle.v1().y - triangle.v2().y)*(triangle.v0().x - triangle.v2().x) + (triangle.v2().x - triangle.v1().x)*(triangle.v0().y - triangle.v2().y));
									
	float beta = ((triangle.v2().y - triangle.v0().y)*(point.x - triangle.v2().x) + (triangle.v0().x - triangle.v2().x)*(point.y - triangle.v2().y)) / 
	((triangle.v1().y - triangle.v2().y)*(triangle.v0().x - triangle.v2().x) + (triangle.v2().x - triangle.v1().x)*(triangle.v0().y - triangle.v2().y));
		
	float gamma = 1.0f - alpha - beta;

	if(alpha >=0 && beta >=0 && gamma >=0){
		return true;
	} else return false;
} 

void drawFilledTriangle(DrawingWindow &window){
	float red = rand()%256;
	float green = rand()%256;
	float blue = rand()%256;

	uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);

	CanvasTriangle triangle(CanvasPoint(rand()%320 + 1, rand()%240 + 1), CanvasPoint(rand()%320 + 1, rand()%240 + 1), CanvasPoint(rand()%320 + 1, rand()%240 + 1));

	drawline(triangle.v0(), triangle.v1(), colour, window);
	drawline(triangle.v1(), triangle.v2(), colour, window);
	drawline(triangle.v2(), triangle.v0(), colour, window);

	std::vector<CanvasPoint> yRange;

	yRange.push_back(triangle.v0());
	yRange.push_back(triangle.v1());
	yRange.push_back(triangle.v2());

	//order from lowest to highest Y
	if(yRange[1].y < yRange[0].y){
		std::swap(yRange[1], yRange[0]);
	}

	if(yRange[2].y < yRange[0].y){
		std::swap(yRange[2], yRange[0]);
	}

	if(yRange[2].y < yRange[1].y){
		std::swap(yRange[2], yRange[1]);
	}

	std::vector<CanvasPoint> xRange = yRange;

	//sort xrange
	if(xRange[1].x < xRange[0].x){
		std::swap(xRange[1], xRange[0]);
	}

	if(xRange[2].x < xRange[0].x){
		std::swap(xRange[2], xRange[0]);
	}

	if(xRange[2].x < xRange[1].x){
		std::swap(xRange[2], xRange[1]);
	}

	TextureMap ttx(texture.ppm);

	
	

	for (int y = yRange[0].y; y<=yRange[2].y; y++){
		for (int x = xRange[0].x; x<=xRange[2].x; x++){
			bool inside = true;
			inside = insideTriangle(triangle, CanvasPoint(x, y));
			
			if(inside){
				//std::cout << inside << std::endl;
				window.setPixelColour(x, y, colour);
			} 

		}
	}


}

void draw(DrawingWindow &window) {
	//window.clearPixels();



	

	//drawTriangle(window);

	/*for (size_t y = 0; y < window.height; y++) {
		
		for (size_t x = 0; x < window.width; x++) {
			float red = 0;
			float green =0;
			float blue = 0;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);

			if(window.getPixelColour(x, y) != colour ){
				window.setPixelColour(x, y, colour);
			}
		}
	}*/
	
	
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
		else if (event.key.keysym.sym == SDLK_u) {
			std::cout << "u" << std::endl; 
			drawTriangle(window);
		}
		else if (event.key.keysym.sym == SDLK_f) {
			std::cout << "f" << std::endl; 
			drawFilledTriangle(window);
		}
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

