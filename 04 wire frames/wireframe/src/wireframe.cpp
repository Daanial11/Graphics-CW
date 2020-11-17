#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <ModelTriangle.h>
#include <sstream>

#define WIDTH 320
#define HEIGHT 240

std::vector<float> interpolateSingleFloats(float firstNo, float secondNo, int listSize);
std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues);
bool checkDepth(CanvasPoint point, std::vector<std::vector<float>> &depthBuffer, CanvasTriangle triangle, Colour colour);
std::vector<CanvasPoint> getLineCords(CanvasPoint startPoint, CanvasPoint endPoint);


//Draws line between two points, uses depthbuffer if depth values are given in depthbuffer or Canvaspoints.
void drawline(CanvasPoint startPoint, CanvasPoint endPoint, uint32_t colour, DrawingWindow &window, std::vector<std::vector<CanvasPoint>> &depthBuffer ){
	 std::vector<CanvasPoint> cords(getLineCords(startPoint, endPoint));
		
		for(auto&element:cords){
			if(element.depth >= depthBuffer[element.x][element.y].depth){
				 window.setPixelColour(round(element.x), round(element.y), colour);
			}
		}
		
}
//loads OBJ material data such as a colour palette
std::vector<Colour> loadOBJmaterial(const std::string &filename){
	std::ifstream in(filename, std::ios::in);

	if (!in) {
        std::cerr << "Cannot open " << filename << std::endl;
        exit(1);

    }

	std::string line;
	std::vector<Colour> palette;

	while(std::getline(in, line)){
		Colour currentColour;
		if(line.substr(0,7) == "newmtl "){
			std::istringstream c;
			c.str(line.substr(7));

			std::string colourName;
			c>>colourName;
			currentColour.name = colourName;

			std::getline(in, line);

			std::istringstream rgbLine;

			rgbLine.str(line.substr(2));
			float r, g, b;
			rgbLine >>r >> g >> b;
			currentColour.blue = b*255;
			currentColour.red = r*255;
			currentColour.green = g*255;

			//std::cout <<currentColour<<std::endl;

		}

		palette.push_back(currentColour);
	}

	return palette;


	
}

//loads OBJ model data
std::vector<ModelTriangle> loadOBJ(const std::string &filename){
	std::vector<Colour> palette = loadOBJmaterial("cornell-box.mtl");
	std::ifstream in(filename, std::ios::in);

	if (!in) {
        std::cerr << "Cannot open " << filename << std::endl;
        exit(1);

    }

	std::string line;
	std::vector<glm::vec3> vertexList;

	std::vector<ModelTriangle> result;

	std::string currentColour;

	while(std::getline(in, line)){

		
		if(line.substr(0,2) == "v "){
			
			std::istringstream v;
			v.str(line.substr(2));

			
			double x, y, z;
			v>>x;
			v>>y;
			v>>z;

			glm::vec3 vertex(x*0.17, y*0.17, z*0.17);

			vertexList.push_back(vertex);




		}
		else if(line.substr(0,7) == "usemtl "){
			std::istringstream c;
			c.str(line.substr(7));

			std::string colourName;
			c>>colourName;
			currentColour = colourName;
		}

		else if(line.substr(0, 2) == "f "){
			std::istringstream f;
			f.str(line.substr(2));

			int a, b, c;

			f>>a;
			f.ignore(1, '/');
			f>>b;
			f.ignore(1, '/');
			f>>c;
			

			ModelTriangle triangle;
			triangle.vertices = {vertexList[a -1], vertexList[b -1], vertexList[c -1]};

			for(auto&element:palette){
				if(currentColour.compare(element.name) == 0){
					triangle.colour=element;
				}

			}
			result.push_back(triangle);

			//std::cout << triangle.colour << std::endl;
		}

		
	}


	
	return result;

}

//checks if a given point lies within a CanvasTriangle.
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

//transforms point in 3D space to a point in 2D space which represents the original from a certain camera angle and focallenght
CanvasPoint getUV(glm::vec3 point){

	glm::vec3 cameraPosition(0.0, 0.0, 2);
	float imageScaler = 180;
	float focalLength = 2;

	glm::vec3 cameraSystemPoint = point - cameraPosition;

	float u = imageScaler*(focalLength * (cameraSystemPoint.x / cameraSystemPoint.z)) + (WIDTH/2);
    float v = imageScaler*(focalLength * (cameraSystemPoint.y / cameraSystemPoint.z)) + (HEIGHT/2);

	//HERE CALCULATE THE Z VALUE 
 
	//width-xvalue to flip everything into correct perspective, each canvas point also has a depth 1/z which will be used later
    return CanvasPoint(WIDTH - int(u), int(v), 1 / point.z);
}

//converts a Colour colour into uint32_t value used by setpixel function
uint32_t colourConvert(Colour colour){
	uint32_t thisColour = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);

	return thisColour;

}
//draws empty triangle/outline, uses depth buffer if values are populated in the buffer and canvas points
void drawTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour, std::vector<std::vector<CanvasPoint>> &depthBuffer) {


	uint32_t thisColour = colourConvert(colour);


	drawline(triangle.v0(), triangle.v1(), thisColour, window, depthBuffer);
	drawline(triangle.v1(), triangle.v2(), thisColour, window, depthBuffer);
	drawline(triangle.v2(), triangle.v0(), thisColour, window, depthBuffer);

	

}


//interpolates depth values 1/z between two points.
std::vector<CanvasPoint> getLineCords(CanvasPoint startPoint, CanvasPoint endPoint){

	std::vector<CanvasPoint> result;
	float xDiff = endPoint.x - startPoint.x;
	float yDiff = endPoint.y - startPoint.y;

	float steps = std :: max(abs(xDiff), abs(yDiff));

	float xStep = xDiff/steps;
	float yStep = yDiff/steps;
	std::vector<float> zValues(interpolateSingleFloats(startPoint.depth, endPoint.depth, steps));

	for (float i=0.0; i<steps; i++) {
		float x = startPoint.x + (xStep*i);
		float y = startPoint.y + (yStep*i);
		CanvasPoint cord(round(x), round(y));
		
		cord.depth = zValues[i];

		result.push_back(cord);
	}

	return result;
}

//produces a height*width depth buffer with just edge values of the triangle in question, all other depth values=0, used later to interpolate depth values for pixels inside triangle.
std::vector<std::vector<CanvasPoint>> getEdgeDepthBuffer(CanvasTriangle triangle){

	std::vector<std::vector<CanvasPoint>> depthBuffer(WIDTH);

	for ( int i = 0 ; i < WIDTH ; i++ ){
   		depthBuffer[i].resize(HEIGHT);
	}
	//getting depth data for points along edges
	std::vector<CanvasPoint> v0v1(getLineCords(triangle.v0(), triangle.v1()));
	std::vector<CanvasPoint> v1v2(getLineCords(triangle.v1(), triangle.v2()));
	std::vector<CanvasPoint> v2v0(getLineCords(triangle.v2(), triangle.v0()));
	
	//add interpolated edge depths to local depth buffer
	for(auto&element:v0v1){
		depthBuffer[element.x][element.y] = element;
	}

	for(auto&element:v1v2){
		depthBuffer[element.x][element.y] = element;
	}

	for(auto&element:v2v0){
		depthBuffer[element.x][element.y] = element;
	}

	return depthBuffer;
}

//draws a filled triangle with colour of choice, updates master depthbuffer with new values if needed.
void drawFilledTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour, std::vector<std::vector<CanvasPoint>> &depthBuffer){
	
	

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

	//order from lowest to highest X, left to right
	if(xRange[1].x < xRange[0].x){
		std::swap(xRange[1], xRange[0]);
	}

	if(xRange[2].x < xRange[0].x){
		std::swap(xRange[2], xRange[0]);
	}

	if(xRange[2].x < xRange[1].x){
		std::swap(xRange[2], xRange[1]);
	}


	std::vector<std::vector<CanvasPoint>> localDepthBuffer = getEdgeDepthBuffer(triangle);
	
	//go through imaginery square around xDiff and yDiff, interpolate the depth values and draw each pixel in a given row starting from the entry pixel, ending with the exit pixel out of the triangle.
	for (int y = yRange[0].y; y<=yRange[2].y; y++){
		std::vector<CanvasPoint> currentRowInside;
		

		for (int x = xRange[0].x; x<=xRange[2].x; x++){

			if(localDepthBuffer[x][y].depth != 0){
				currentRowInside.push_back(localDepthBuffer[x][y]);
			}



	

		}

		//getting interpolated depth values between entry and exit pixels into triangle for given row
		currentRowInside = getLineCords(currentRowInside[0], currentRowInside[currentRowInside.size()-1]);

		//compare local depth buffer with data of current triangle we'r trying to draw with master depth buffer values. if something is already ahead of it (higher 1/z/closer to camera) then it won't setpixel.
		for(auto&element:currentRowInside){
			
			if(depthBuffer[element.x][element.y].depth == 0){
				depthBuffer[element.x][element.y].depth = element.depth;
				window.setPixelColour(element.x, element.y, colourConvert(colour));
			}

			else if(depthBuffer[element.x][element.y].depth < element.depth){
				depthBuffer[element.x][element.y].depth = element.depth;
				window.setPixelColour(element.x, element.y, colourConvert(colour));
			}

			else if(depthBuffer[element.x][element.y].depth < 0){
				depthBuffer[element.x][element.y].depth = element.depth;
				window.setPixelColour(element.x, element.y, colourConvert(colour));
			}
		}
	}

	//draw outline of each triangle, done at end as the depthbuffer is fully populated
	drawTriangle(window, triangle, colour, depthBuffer);


}
/*
bool checkDepth(CanvasPoint point, std::vector<std::vector<float>> &depthBuffer, CanvasTriangle triangle, Colour colour){

	double d1 = sqrt(pow(triangle.v0().x - point.x, 2) + pow(triangle.v0().y - point.y, 2));
	double d2 = sqrt(pow(triangle.v1().x - point.x, 2) + pow(triangle.v1().y - point.y, 2));
	double d3 = sqrt(pow(triangle.v2().x - point.x, 2) + pow(triangle.v2().y - point.y, 2));

	double totalD = d1 + d2 + d3;

	double totalZ = triangle.v0().depth + triangle.v1().depth + triangle.v2().depth;

	double pixelValue = (100- (((((d1/totalD) * triangle.v0().depth) + ((d2/totalD) * triangle.v1().depth) + ((d3/totalD) * triangle.v2().depth)) / totalZ) * 100));

	if(point.x ==119 && point.y == 144){
		std::cout<<pixelValue<<","<< colour.name<< (depthBuffer[point.x][point.y] < pixelValue) <<std::endl;
		std::cout<<"d1: "<<d1<<" d2: "<<d2<<" d3: "<<d3<< " totalZ: "<< totalZ <<std::endl;
		std::cout<<"z1: "<<triangle.v0().depth<<" z2: "<<triangle.v1().depth<<" z3: "<<triangle.v1().depth <<std::endl;
	}
	if(depthBuffer[point.x][point.y] == 0){
		depthBuffer[point.x][point.y] = pixelValue;
		//std::cout << pixelValue<<std::endl;	
		return true;
	} 
	else if(depthBuffer[point.x][point.y] > pixelValue){
		return false;
	}
	else if(depthBuffer[point.x][point.y] < pixelValue){
		depthBuffer[point.x][point.y] = pixelValue;	
		return true;
	}
}*/

//projects a triangle in model coordinate system onto an image plane or from a certain camera angle.
void projectImage(DrawingWindow &window, std::vector<ModelTriangle> vertices){

	//intialize heightxWidth depth buffer.
	std::vector<std::vector<CanvasPoint>> depthBuffer(WIDTH);

	for ( int i = 0 ; i < WIDTH ; i++ ){
   		depthBuffer[i].resize(HEIGHT);
	}


	//changing coordinate systems//getting projected triangle
	for (auto&element: vertices){
		CanvasPoint P1 = getUV(element.vertices[0]);
		//P1.depth = (element.vertices[0]).z;
		CanvasPoint P2 = getUV(element.vertices[1]);
		//P2.depth = (element.vertices[1]).z;
		CanvasPoint P3 = getUV(element.vertices[2]);
		//P3.depth = (element.vertices[2]).z;

		

		
		drawFilledTriangle(window, CanvasTriangle(P1, P2, P3), element.colour, depthBuffer);

	
		//std::cout<<depthBuffer[163][193].depth<<"," <<std::endl;
			
			
		
	}

}



void draw(DrawingWindow &window) {
	std::vector<ModelTriangle> faces = loadOBJ("cornell-box.obj"); 
	projectImage(window, faces);

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
			//indow.setPixelColour(x, y, colour);
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

	
	//std::vector<Colour> palette = loadOBJmaterial("cornell-box.mtl");


	
		

	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		update(window);
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}

