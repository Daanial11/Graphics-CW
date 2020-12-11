#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <ModelTriangle.h>
#include <sstream>
#include <math.h>
#include <iostream>
#include <TextureMap.h>
#include <map>
#include <RayTriangleIntersection.h>
#include <algorithm>
#include <thread>





#define WIDTH 640
#define HEIGHT 480
float camX = 0.0;
float camY = 0.0;
float camZ = 1.5;
float xAngle = 0.0;
float yAngle = 0.0;
float imageScaler = 100.0;
float t = 0;
int mode = 2;
int counter = 0;
int camCounter = 0;
int camCounter2 = 0;
int camCounter3 = 0;
int camCounter4 = 0;
int camCounter5 = 0;
float lightX = 0.0;
float lightY = 0.25;
float lightZ = 0.5;

bool lookatc = true;
TextureMap texture;

std::vector<float> interpolateSingleFloats(float firstNo, float secondNo, int listSize);
std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues);
bool checkDepth(CanvasPoint point, std::vector<std::vector<float>> &depthBuffer, CanvasTriangle triangle, Colour colour);
std::vector<CanvasPoint> getLineCords(CanvasPoint startPoint, CanvasPoint endPoint);
Colour refract(glm::vec3 surfaceNormal, RayTriangleIntersection incidentRay, glm::vec3 rayDir, glm::vec3 intersectionCord, float OutIndex, float InIndex, std::vector<ModelTriangle>& vertices, glm::vec3 light);


//Draws line between two points, uses depthbuffer if depth values are given in depthbuffer or Canvaspoints.
void drawline(CanvasPoint startPoint, CanvasPoint endPoint, uint32_t colour, DrawingWindow &window, std::vector<std::vector<CanvasPoint>> &depthBuffer ){
	 std::vector<CanvasPoint> cords(getLineCords(startPoint, endPoint));
		
		for(auto&element:cords){
			if(element.x < WIDTH && element.x>=0.0 && element.y>=0 && element.y < HEIGHT){
				if(element.depth >= depthBuffer[element.x][element.y].depth){
					window.setPixelColour(round(element.x), round(element.y), colour);
			}
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

void normalProvidedExtract(std::string line, int& a, int& b, int& c) {
	std::istringstream f;
	f.str(line.substr(2));
	char x;
	int l;
	f >> a;
	f >> x;
	f >> x;
	f >> l;
	f >> b;
	f >> x;
	f >> x;
	f >> l;
	f >> c;
	
}

//loads OBJ model data
std::vector<ModelTriangle> loadOBJ(const std::string &filename){
	std::vector<Colour> palette = loadOBJmaterial("textured-cornell-box.mtl");
	std::ifstream in(filename, std::ios::in);

	if (!in) {
        std::cerr << "Cannot open " << filename << std::endl;
        exit(1);

    }

	std::string line;
	std::vector<glm::vec3> vertexList;
	std::vector<glm::vec3> normalList;
	std::vector<TexturePoint> texturePointList;

	std::vector<ModelTriangle> result;

	std::string currentColour;
	bool textureExist = false;
	bool normalProvided = false;

	while(std::getline(in, line)){
		
		
		
		if (line.substr(0, 3) == "vn ") {
			normalProvided = true;
			std::istringstream c;
			c.str(line.substr(3));

			float x, y, z;
			c >> x;
			c >> y;
			c >> z;

			normalList.push_back(glm::vec3(x, y, z));
		}
		
		else if(line.substr(0,2) == "v "){
			textureExist = false;
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

		else if(line.substr(0, 3) == "vt "){
			std::istringstream c;
			c.str(line.substr(3));

			float x, y;

			c>>x;
			c>>y;

			TexturePoint point(x*texture.width, y*texture.height);

			

			texturePointList.push_back(point);

			textureExist = true;
		}

		else if(line.substr(0, 2) == "f "){
			
			std::istringstream f;
			f.str(line.substr(2));

			int a, b, c;
			int f1, f2, f3;

			if (!normalProvided) {
				f >> a;
				f.ignore(1, '/');
				if (textureExist) {
					f >> f1;
				}
				f >> b;
				f.ignore(1, '/');
				if (textureExist) {
					f >> f2;
				}
				f >> c;

				if (textureExist) {
					f.ignore(1, '/');
					f >> f3;
				}
			}
			else {
				normalProvidedExtract(line, a, b, c);
				
			}
			

			ModelTriangle triangle;
			triangle.vertices = {vertexList[a -1], vertexList[b -1], vertexList[c -1]};
			
			triangle.normal = glm::cross((triangle.vertices[1] - triangle.vertices[0]), (triangle.vertices[2] - triangle.vertices[0]));
			
			if (normalProvided) {
				triangle.vertexNormals = { normalList[a - 1], normalList[b - 1], normalList[c - 1] };
				triangle.vertexNormalsProv = true;
			}
			else triangle.vertexNormalsProv = false;

			if(textureExist){
				triangle.texturePoints = {texturePointList[f1-1], texturePointList[f2-1], texturePointList[f3-1]};
				
			}

			for(auto&element:palette){
				if(currentColour.compare(element.name) == 0){
					triangle.colour=element;
				}

			}

			if (normalProvided) {
				triangle.colour = palette[4];
			}

			result.push_back(triangle);

		
		}

		
	}


	
	return result;

}

bool validateIntersection(float u, float v, float t) {
	

	//std::cout << "u: " << u << " v:" << v << " t:" << t << std::endl;
	if (t > 0.0 && (u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && ((double)u + (double)v <= 1.0)) {
		
		return true;
	}

	return false;
}

RayTriangleIntersection getClosestIntersection(glm::vec3 camPos, glm::vec3 rayDir, std::vector<ModelTriangle>& triangles, int dontInclude) {

	std::vector<RayTriangleIntersection> solutions;

	for (int i = 0; i < triangles.size(); i++) {
		glm::vec3 e0 = triangles[i].vertices[1] - triangles[i].vertices[0];
		glm::vec3 e1 = triangles[i].vertices[2] - triangles[i].vertices[0];
		glm::vec3 startPoint = camPos - triangles[i].vertices[0];
		//glm::vec3 rayDir(glm::normalize(canvasPoint - startPoint));
		glm::mat3 DEMatrix(-rayDir, e0, e1);
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * startPoint;

		
		if (validateIntersection(possibleSolution[1], possibleSolution[2], possibleSolution[0]*-1) && i!=dontInclude) {
			solutions.push_back(RayTriangleIntersection(possibleSolution, possibleSolution[0] * -1, triangles[i], i));
		
			
		}

	}

	

	RayTriangleIntersection closestSolution;
	closestSolution.distanceFromCamera = 1000.0;
	
	
	for (const RayTriangleIntersection& element : solutions) {
		//std::cout << element.intersectionPoint[1] << std::endl;
		if (element.distanceFromCamera < closestSolution.distanceFromCamera) {

			closestSolution = element;
			//std::cout << closestSolution.intersectedTriangle.colour <<closestSolution.distanceFromCamera << std::endl;
		}
	}
	
	
	
	//std::cout << closestSolution.distanceFromCamera << std::endl;
	return closestSolution;
}

bool shadowCheck(RayTriangleIntersection intersection, glm::vec3 light, std::vector<ModelTriangle>& triangles) {

	int numberOfSolutions = 0;

	glm::vec3 v0(intersection.intersectedTriangle.vertices[0]);
	glm::vec3 v1(intersection.intersectedTriangle.vertices[1]);
	glm::vec3 v2(intersection.intersectedTriangle.vertices[2]);

	float u = intersection.intersectionPoint[1];
	float v = intersection.intersectionPoint[2];
	
	//cordinate in 3-d space
	auto cord = v0 + u*(v1 - v0) + v*(v2 - v0);
	
	auto shadowRayDir = glm::normalize(light - cord);
	auto distance = glm::distance(light, cord);

	for (int i = 0; i < triangles.size(); i++) {
		glm::vec3 e0 = triangles[i].vertices[1] - triangles[i].vertices[0];
		glm::vec3 e1 = triangles[i].vertices[2] - triangles[i].vertices[0];
		glm::vec3 startPoint = cord - triangles[i].vertices[0];
		
		glm::mat3 DEMatrix(-shadowRayDir, e0, e1);
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * startPoint;

		

		//shadow ray is currently considering objects past the intersection point, so a point on the box and the point on the box on the other side, fix this
		//check for a known point where you expect a shadow, check the intersection vector, for both shadowray and normal ray, find a way with this to fix above
		if (validateIntersection(possibleSolution[1], possibleSolution[2], possibleSolution[0]) &&  possibleSolution[0] < distance && intersection.triangleIndex != i) {
			numberOfSolutions++;
			
		}

	}



	if (numberOfSolutions > 0) {
		return true;
	}
	return false;
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



glm::mat3 lookAt(glm::vec3 o, glm::vec3 camPos){
	glm::vec3 arbitary(0.0, 1.0, 0.0);



	glm::vec3 forward = glm::normalize(o - camPos);
	glm::vec3 right = glm::normalize(glm::cross(arbitary, forward));
	glm::vec3 up = glm::cross(right, forward);

	glm::mat3 view(right, up, forward);

	return view;
}



//transforms point in 3D space to a point in 2D space which represents the original from a certain camera angle and focallenght
CanvasPoint getUV(glm::vec3 point){

	glm::vec3 cameraPosition(camX, camY, camZ);

	

	
	float focalLength = 2.0;

	//rotate along y axis in the space
	glm::vec3 column0(cos(yAngle), 0.0, sin(yAngle));
	glm::vec3 column1(0.0, 1.0, 0.0);
	glm::vec3 column2(-sin(yAngle), 0.0, cos(yAngle));
	glm::mat3 n(column0, column1, column2); // sets columns of matrix n

	//rotate along x axis in the space
	glm::vec3 column3(1.0, 0.0, 0.0);
	glm::vec3 column4(0.0, cos(xAngle), -sin(xAngle));
	glm::vec3 column5(0.0, sin(xAngle), cos(xAngle));
	glm::mat3 n1(column3, column4, column5); // sets columns of matrix n

	glm::mat3 camOri = n*n1;

	
	

	glm::vec3 adjustedVector = point - cameraPosition;

	adjustedVector = camOri * adjustedVector ;
	
	if (lookatc) {
		adjustedVector = adjustedVector * lookAt(glm::vec3{ 0, 0, 0 }, cameraPosition);
	}
	
	
	float u = imageScaler*(focalLength * (adjustedVector.x / adjustedVector.z)) + (WIDTH/2.0);
    float v = imageScaler*(focalLength * (adjustedVector.y / adjustedVector.z)) + (HEIGHT/2.0);



 
 
	//width-xvalue to flip everything into correct perspective, each canvas point also has a depth 1/z which will be used later
    return CanvasPoint(WIDTH - int(u), int(v), 1/adjustedVector.z);
}

//converts a Colour colour into uint32_t value used by setpixel function
uint32_t colourConvert(Colour colour){

	
	
	uint32_t thisColour = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);
	
	

	return thisColour;

}
// https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
void Barycentric(glm::vec2 p, glm::vec2 a, glm::vec2 b, glm::vec2 c, float& u, float& v, float& w)
{
	auto v0 = b - a, v1 = c - a, v2 = p - a;
	float d00 = glm::dot(v0, v0);
	float d01 = glm::dot(v0, v1);
	float d11 = glm::dot(v1, v1);
	float d20 = glm::dot(v2, v0);
	float d21 = glm::dot(v2, v1);
	float denom = d00 * d11 - d01 * d01;
	v = (d11 * d20 - d01 * d21) / denom;
	w = (d00 * d21 - d01 * d20) / denom;
	u = 1.0f - v - w;
}

uint32_t getTexturePixel(CanvasTriangle triangle, float x, float y, TextureMap &texture){


	float alpha, beta, gamma;

	Barycentric(glm::vec2(x,y), glm::vec2(triangle.v0().x, triangle.v0().y), glm::vec2(triangle.v1().x, triangle.v1().y), glm::vec2(triangle.v2().x, triangle.v2().y), gamma, alpha, beta);
	

	float interpolatedX = (alpha * triangle.v0().texturePoint.x) + (beta * triangle.v1().texturePoint.x) + (gamma * triangle.v2().texturePoint.x);
	float interpolatedY = (alpha * triangle.v0().texturePoint.y) + (beta * triangle.v1().texturePoint.y) + (gamma * triangle.v2().texturePoint.y);
	
	
	int index = ((int)interpolatedY * texture.width) + (int)interpolatedX;
	
	
	uint32_t thisColour = texture.pixels[index];
	//std::cout << "X: "<< interpolatedX <<" Y: "<<interpolatedY << " colour uint32_t: "<< thisColour << std::endl;
	

	

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
std::vector<CanvasPoint> getEdgeDepthBuffer(CanvasTriangle triangle){

	std::vector<CanvasPoint> result;
	//getting depth data for points along edges
	std::vector<CanvasPoint> v0v1(getLineCords(triangle.v0(), triangle.v1()));
	std::vector<CanvasPoint> v1v2(getLineCords(triangle.v1(), triangle.v2()));
	std::vector<CanvasPoint> v2v0(getLineCords(triangle.v2(), triangle.v0()));
	
	//add interpolated edge depths to local depth buffer
	for(auto&element:v0v1){
		if(element.x < WIDTH && element.x>=0.0 && element.y>=0 && element.y < HEIGHT){
			result.push_back(element);
		}
	}

	for(auto&element:v1v2){
		if(element.x < WIDTH && element.x>=0.0 && element.y>=0 && element.y < HEIGHT){
			result.push_back(element);
		}
	}

	for(auto&element:v2v0){
		if(element.x < WIDTH && element.x>=0.0 && element.y>=0 && element.y < HEIGHT){
			result.push_back(element);
		}
	}

	return result;
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

	//possibly store these in hashmap instead to speed up, for getting current row just query the y value for that row to get all the cords in that row <key, list<x values>>
	std::vector<CanvasPoint> localDepthBuffer = getEdgeDepthBuffer(triangle);
	
	
	//go through imaginery square around xDiff and yDiff, interpolate the depth values and draw each pixel in a given row starting from the entry pixel, ending with the exit pixel out of the triangle.
	for (int y = yRange[0].y; y<=yRange[2].y; y++){
		std::vector<CanvasPoint> currentRowInside;
		
	

		for (int x = xRange[0].x; x<=xRange[2].x; x++){
			if(x < WIDTH && x>=0.0 && y>=0 && y < HEIGHT){
				for (auto& element : localDepthBuffer) {
					if (element.x == x && element.y == y) {
						currentRowInside.push_back(element);
					}
				}
				
			}
			



	

		}


		//getting interpolated depth values between entry and exit pixels into triangle for given row
		if(currentRowInside.size()>0){
			currentRowInside = getLineCords(currentRowInside[0], currentRowInside[currentRowInside.size()-1]);
		}
		

		//compare local depth buffer with data of current triangle we'r trying to draw with master depth buffer values. if something is already ahead of it (higher 1/z/closer to camera) then it won't setpixel.
		for(auto&element:currentRowInside){

			uint32_t colourValue;

			if(triangle.v0().texturePoint.x !=0){
				colourValue = getTexturePixel(triangle, element.x, element.y, texture);
				//std::cout<<colourValue<<std::endl;
			} else colourValue = colourConvert(colour);
			
			if(depthBuffer[element.x][element.y].depth == 0){
				depthBuffer[element.x][element.y].depth = element.depth;
				window.setPixelColour(element.x, element.y, colourValue);
			}

			else if(depthBuffer[element.x][element.y].depth < element.depth){
				depthBuffer[element.x][element.y].depth = element.depth;
				window.setPixelColour(element.x, element.y, colourValue);
			}

			/*else if(depthBuffer[element.x][element.y].depth < 0){
				depthBuffer[element.x][element.y].depth = element.depth;
				window.setPixelColour(element.x, element.y, colourValue);
			}*/
		}
	}

	//draw outline of each triangle, done at end as the depthbuffer is fully populated
	drawTriangle(window, triangle, colour, depthBuffer);


}

float proximityLighting(glm::vec3 cord, glm::vec3 light) {

	auto distance = glm::distance(light, cord);


	float lightMultiplier = (1 - (0.25 * 3.1415 * pow(distance, 1.5)));


	

	return lightMultiplier;
}

float clamp(float value, float low, float max) {
	if (value > max) {
		return max;
	}
	else if (value < low) {
		return low;
	}
	else return value;
}

float incidenceLighting(glm::vec3 normal, glm::vec3 cord, glm::vec3 light) {
	
	auto lightDir = glm::normalize(light-cord);
	

	auto angleOfIncidence = glm::dot(lightDir, glm::normalize(normal));

	angleOfIncidence = clamp(angleOfIncidence, 0.0f, 1.0f);

	
	

	


	return angleOfIncidence;
}

float specLighting(glm::vec3 normal, glm::vec3 cord, glm::vec3 light) {


	auto viewDir = glm::normalize(-1.0f * cord);
	auto lightDir = glm::normalize(light-cord);

	auto norm = glm::normalize(normal);
	auto spec = 0.0f;
	auto reflection = glm::normalize(2.0f * glm::dot(norm, lightDir) * norm - lightDir);

	auto cosAngle = clamp(glm::dot(viewDir, reflection),0.0f, 1.0f);
	float specIntensity = 0.8;

	if (cosAngle == 1.0) {
		spec = 1.0;
	}
	else spec = pow(cosAngle, 2.0f) * specIntensity;

	return spec;

}

Colour singlePointColour(glm::vec3 normal, glm::vec3 cord, glm::vec3 light, Colour currentColour, bool shadowCheck) {
	
	
	float prox = proximityLighting(cord, light);
	float incidence = incidenceLighting(normal, cord, light);

	float spec = specLighting(normal, cord, light);
	float baseBrightness = 0.0f;
	float totalMultiplier = clamp((prox + spec) * incidence, 0.2f, 1.0f);

	
	
	currentColour.blue = currentColour.blue * totalMultiplier;

	currentColour.green = currentColour.green * totalMultiplier;
	currentColour.red = currentColour.red * totalMultiplier;
	

	

	return currentColour;
}
Colour reflectiveMaterial(RayTriangleIntersection ray, glm::vec3 rayDir, std::vector<ModelTriangle>& vertices, glm::vec3 startPoint, glm::vec3 light) {

	auto angleOfIncidence = glm::dot(rayDir, glm::normalize(ray.intersectedTriangle.normal));


	auto reflectDir = rayDir - ((2.0f * glm::normalize(ray.intersectedTriangle.normal)) * angleOfIncidence);

	

	auto reflectInter = getClosestIntersection(startPoint, reflectDir, vertices, ray.triangleIndex);

	//adjusting colour with lighting so reflection looks normal
	glm::vec3 normal = reflectInter.intersectedTriangle.normal;
	auto cord = (reflectInter.intersectedTriangle.vertices[0]) + (reflectInter.intersectionPoint[1] * (reflectInter.intersectedTriangle.vertices[1] - reflectInter.intersectedTriangle.vertices[0])) + (reflectInter.intersectionPoint[2] * (reflectInter.intersectedTriangle.vertices[2] - reflectInter.intersectedTriangle.vertices[0]));
	
	if (reflectInter.intersectedTriangle.colour.name.compare("Red") == 0) {
		reflectInter.intersectedTriangle.colour = refract(normal, reflectInter, reflectDir, cord, 1.0, 1.5, vertices, light);
	}
	else {
		float prox = proximityLighting(cord, light);
		float incidence = incidenceLighting(normal, cord, light);
		float spec = specLighting(normal, cord, light);
		float baseBrightness = 0.0f;
		float totalMultiplier = clamp((prox + spec) * incidence, 0.2f, 1.0f);
		reflectInter.intersectedTriangle.colour.blue = reflectInter.intersectedTriangle.colour.blue * totalMultiplier;
		reflectInter.intersectedTriangle.colour.green = reflectInter.intersectedTriangle.colour.green * totalMultiplier;
		reflectInter.intersectedTriangle.colour.red = reflectInter.intersectedTriangle.colour.red * totalMultiplier;
	}


	return reflectInter.intersectedTriangle.colour;

}

Colour gouraudShade(RayTriangleIntersection ray, glm::vec3 light, float u, float v, std::vector<ModelTriangle>& vertices) {

	
	glm::vec3 cam(camX, camY, camZ);
	
	//check if triangle vertexs are in shadow before interpolating for the point

	auto rayDIR0 = glm::normalize((ray.intersectedTriangle.vertices[0]/imageScaler) - cam);
	RayTriangleIntersection ray0(getClosestIntersection(glm::vec3(camX, camY, camZ), rayDIR0, vertices, -1));
	bool shadowV0 = shadowCheck(ray0, light, vertices);

	auto rayDIR1 = glm::normalize((ray.intersectedTriangle.vertices[1] / imageScaler) - cam);
	RayTriangleIntersection ray1(getClosestIntersection(glm::vec3(camX, camY, camZ), rayDIR1, vertices, -1));
	bool shadowV1 = shadowCheck(ray1, light, vertices);

	auto rayDIR2 = glm::normalize((ray.intersectedTriangle.vertices[2] / imageScaler) - cam);
	RayTriangleIntersection ray2(getClosestIntersection(glm::vec3(camX, camY, camZ), rayDIR2, vertices, -1));
	bool shadowV2 = shadowCheck(ray2, light, vertices);
	
	
	//get shaded colour values using normals
	Colour v0Colour = singlePointColour(ray.intersectedTriangle.vertexNormals[0], ray.intersectedTriangle.vertices[0], light, ray.intersectedTriangle.colour, shadowV0);
	Colour v1Colour = singlePointColour(ray.intersectedTriangle.vertexNormals[1], ray.intersectedTriangle.vertices[1], light, ray.intersectedTriangle.colour, shadowV1);
	Colour v2Colour = singlePointColour(ray.intersectedTriangle.vertexNormals[2], ray.intersectedTriangle.vertices[2], light, ray.intersectedTriangle.colour, shadowV2);

	

	Colour interpolatedColour;


	//use barycentric cords
	interpolatedColour.blue = v0Colour.blue + (u * (v1Colour.blue - v0Colour.blue)) + (v * (v2Colour.blue - v0Colour.blue));
	interpolatedColour.red =  v0Colour.red + (u * (v1Colour.red - v0Colour.red)) + (v * (v2Colour.red - v0Colour.red));
	interpolatedColour.green = v0Colour.green + (u * (v1Colour.green - v0Colour.green)) + (v * (v2Colour.green - v0Colour.green));

	

	return interpolatedColour;
}

Colour phongShade(glm::vec3 point, RayTriangleIntersection ray, glm::vec3 light, float u, float v, std::vector<ModelTriangle>& vertices) {

	glm::vec3 pointNormal(ray.intersectedTriangle.vertexNormals[0] + (u * (ray.intersectedTriangle.vertexNormals[1] - ray.intersectedTriangle.vertexNormals[0])) + (v * (ray.intersectedTriangle.vertexNormals[2] - ray.intersectedTriangle.vertexNormals[0])));
	bool shadow = shadowCheck(ray, light, vertices);

	Colour c = singlePointColour(pointNormal, point, light, ray.intersectedTriangle.colour, shadow);

	return c;

}

float getReflectance(glm::vec3 surfaceNormal, glm::vec3 incidentRay, float OutIndex, float InIndex, float cosAngle, float sin2T) {
	float indexRatio = OutIndex / InIndex;

	if (sin2T > 1.0) {
		return 1.0;
	}

	float cosT = sqrt(1.0 - sin2T);
	float rOrtho = ((OutIndex * cosAngle) - (InIndex * cosT)) / ((OutIndex * cosAngle) + (InIndex * cosT));
	float rPara = ((InIndex * cosAngle) - (OutIndex * cosT)) / ((InIndex * cosAngle) + (OutIndex * cosT));



	return (pow(rOrtho, 2.0) + pow(rPara, 2.0))/4.0;

}


Colour refract(glm::vec3 surfaceNormal, RayTriangleIntersection incidentRay, glm::vec3 rayDir, glm::vec3 intersectionCord, float OutIndex, float InIndex, std::vector<ModelTriangle>& vertices, glm::vec3 light){
	
	float indexRatio = OutIndex / InIndex;
	float cosAngle = -glm::dot(rayDir, surfaceNormal);
	float sin2T = sqrt((1 - pow(indexRatio, 2.0)) * (1.0 - pow(cosAngle, 2.0)));


	glm::vec3 refractedRayDir = (indexRatio * rayDir) + (indexRatio * (cosAngle - sin2T)) * surfaceNormal;
	
	RayTriangleIntersection refractedRay = getClosestIntersection(intersectionCord, refractedRayDir, vertices, incidentRay.triangleIndex);

	float indexRatio2 = InIndex / OutIndex;
	float cosAngle2 = -glm::dot(refractedRayDir, refractedRay.intersectedTriangle.normal);
	float sin2T2 = sqrt((1 - pow(indexRatio, 2.0)) * (1.0 - pow(cosAngle, 2.0)));

	glm::vec3 refractedRayDir2 = (indexRatio2 * refractedRayDir) + (indexRatio2 * (cosAngle2 - sin2T2)) * refractedRay.intersectedTriangle.normal;
	
	auto cord = (refractedRay.intersectedTriangle.vertices[0]) + (refractedRay.intersectionPoint[1] * (refractedRay.intersectedTriangle.vertices[1] - refractedRay.intersectedTriangle.vertices[0])) + (refractedRay.intersectionPoint[2] * (refractedRay.intersectedTriangle.vertices[2] - refractedRay.intersectedTriangle.vertices[0]));

	
	RayTriangleIntersection refractedRay2 = getClosestIntersection(cord, refractedRayDir2, vertices, refractedRay.triangleIndex);

	



	auto cord2 = (refractedRay2.intersectedTriangle.vertices[0]) + (refractedRay2.intersectionPoint[1] * (refractedRay2.intersectedTriangle.vertices[1] - refractedRay2.intersectedTriangle.vertices[0])) + (refractedRay2.intersectionPoint[2] * (refractedRay2.intersectedTriangle.vertices[2] - refractedRay2.intersectedTriangle.vertices[0]));
	
	if (refractedRay2.triangleIndex == 31 || refractedRay2.triangleIndex == 26) {
		refractedRay2.intersectedTriangle.colour = reflectiveMaterial(refractedRay2, refractedRayDir2, vertices, cord2, light);
	}
	else {
		float prox = proximityLighting(cord2, light);
		float incidence = incidenceLighting(refractedRay2.intersectedTriangle.normal, cord2, light);
		float spec = specLighting(refractedRay2.intersectedTriangle.normal, cord2, light);
		float baseBrightness = 0.0f;
		float totalMultiplier = clamp((prox + spec) * incidence, 0.2f, 1.0f);

		float fren = clamp(getReflectance(surfaceNormal, rayDir, OutIndex, InIndex, cosAngle, sin2T), 0.6, 1.0f);

		

		if (refractedRay.intersectedTriangle.colour.name.compare("Green") == 0) {
			refractedRay2.intersectedTriangle.colour = refractedRay.intersectedTriangle.colour;
		}
		refractedRay2.intersectedTriangle.colour.blue = refractedRay2.intersectedTriangle.colour.blue * totalMultiplier * fren;
		refractedRay2.intersectedTriangle.colour.green = refractedRay2.intersectedTriangle.colour.green * totalMultiplier * fren;
		refractedRay2.intersectedTriangle.colour.red = refractedRay2.intersectedTriangle.colour.red * totalMultiplier * fren;
	}
	


	

	return refractedRay2.intersectedTriangle.colour;

}




void projectImageRayTrace(DrawingWindow& window, std::vector<ModelTriangle> &vertices) {
	float focalLength = 2.0;

	glm::vec3 column0(cos(yAngle), 0.0, sin(yAngle));
	glm::vec3 column1(0.0, 1.0, 0.0);
	glm::vec3 column2(-sin(yAngle), 0.0, cos(yAngle));
	glm::mat3 n(column0, column1, column2); // sets columns of matrix n

	//rotate along x axis in the space
	glm::vec3 column3(1.0, 0.0, 0.0);
	glm::vec3 column4(0.0, cos(xAngle), -sin(xAngle));
	glm::vec3 column5(0.0, sin(xAngle), cos(xAngle));
	glm::mat3 n1(column3, column4, column5); // sets columns of matrix n

	glm::mat3 camOri = n * n1 ;

	
	
	glm::vec3 light(lightX, lightY, lightZ);

	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			glm::vec3 camSpacePixel(((int)WIDTH / 2 - (int)x) / (imageScaler), ((int)y - (int)HEIGHT / 2) / (imageScaler), focalLength);

			

			auto worldSpacePixel = (camSpacePixel * glm::inverse(camOri)) + glm::vec3(camX, camY, camZ);
			auto rayDIR = glm::normalize(worldSpacePixel - glm::vec3(camX, camY, camZ));
			
			
			RayTriangleIntersection ray(getClosestIntersection(glm::vec3(camX, camY, camZ), rayDIR, vertices, -1));
			
			RayTriangleIntersection shadowRay;

			//1000 is just a default value for no intersection
			if (ray.distanceFromCamera != 1000.0) {
				
				bool shadow = shadowCheck(ray, light, vertices);
				Colour adjustedColour;
				glm::vec3 v0(ray.intersectedTriangle.vertices[0]);
				glm::vec3 v1(ray.intersectedTriangle.vertices[1]);
				glm::vec3 v2(ray.intersectedTriangle.vertices[2]);

				float u = ray.intersectionPoint[1];
				float v = ray.intersectionPoint[2];
				//cordinate in 3-d space
				auto cord = (v0) + (u * (v1 - v0)) + (v * (v2 - v0));
				
				if (ray.intersectedTriangle.vertexNormalsProv) {
					//adjustedColour = gouraudShade(ray, light, u, v, vertices);
					adjustedColour = phongShade(cord, ray, light, u, v, vertices);

				}
				else if (ray.intersectedTriangle.colour.name.compare("Red") == 0) {
					adjustedColour = refract(ray.intersectedTriangle.normal, ray, rayDIR, cord, 1.0f, 1.5f, vertices, light);
				}
				else if (ray.triangleIndex == 31 || ray.triangleIndex == 26) {
					adjustedColour = reflectiveMaterial(ray, rayDIR, vertices, cord, light);
				}

				else if (shadow) {
					//if shadow use 1/10th of colour value, flat shadow
					ray.intersectedTriangle.colour.blue = ray.intersectedTriangle.colour.blue * 0.1;
					ray.intersectedTriangle.colour.red = ray.intersectedTriangle.colour.red * 0.1;
					ray.intersectedTriangle.colour.green = ray.intersectedTriangle.colour.green * 0.1;
					adjustedColour = ray.intersectedTriangle.colour;
				}
				else {

					
					glm::vec3 normal = ray.intersectedTriangle.normal;

					

					float prox = proximityLighting(cord, light);
					float incidence = incidenceLighting(normal, cord, light);

					//light intersecting with ceiling for some reason, precision issues
					if (ray.triangleIndex < 4) {
						incidence = 0.7 * prox;
					}
					float spec = specLighting(normal, cord, light);
					float baseBrightness = 0.0f;
					float totalMultiplier = clamp((prox + spec) * incidence, 0.3f, 1.0f);

					

					

					ray.intersectedTriangle.colour.blue = ray.intersectedTriangle.colour.blue * totalMultiplier;
					ray.intersectedTriangle.colour.green = ray.intersectedTriangle.colour.green * totalMultiplier;
					ray.intersectedTriangle.colour.red = ray.intersectedTriangle.colour.red * totalMultiplier;
					
					adjustedColour = ray.intersectedTriangle.colour;
				}

				window.setPixelColour(x, y, colourConvert(adjustedColour));

			}

			
		}
	}
}

//projects a triangle in model coordinate system onto an image plane or from a certain camera angle.
void projectImageRaster(DrawingWindow &window, std::vector<ModelTriangle> vertices){

	//intialize heightxWidth depth buffer. faster to just create it once and clear it everytime project image is called.
	std::vector<std::vector<CanvasPoint>> depthBuffer(WIDTH);
	
	for ( int i = 0 ; i < WIDTH ; i++ ){
   		depthBuffer[i].resize(HEIGHT);
	}
	
	

	//changing coordinate systems//getting projected triangle
	for (auto&element: vertices){
		CanvasPoint P1 = getUV(element.vertices[0]);
		P1.texturePoint = element.texturePoints[0];
		//P1.depth = (element.vertices[0]).z;
		CanvasPoint P2 = getUV(element.vertices[1]);
		P2.texturePoint = element.texturePoints[1];
		//P2.depth = (element.vertices[1]).z;
		CanvasPoint P3 = getUV(element.vertices[2]);
		P3.texturePoint = element.texturePoints[2];

		
		//P3.depth = (element.vertices[2]).z;

		

		
		drawFilledTriangle(window, CanvasTriangle(P1, P2, P3), element.colour, depthBuffer);

	
		//std::cout<<depthBuffer[163][193].depth<<"," <<std::endl;
			
			
		
	}

}

void wireframeMode(DrawingWindow& window, std::vector<ModelTriangle> vertices) {

	std::vector<std::vector<CanvasPoint>> depthBuffer(WIDTH);

	for (int i = 0; i < WIDTH; i++) {
		depthBuffer[i].resize(HEIGHT);
	}

	for (auto& element : vertices) {
		CanvasPoint P1 = getUV(element.vertices[0]);


		CanvasPoint P2 = getUV(element.vertices[1]);


		CanvasPoint P3 = getUV(element.vertices[2]);








		drawTriangle(window, CanvasTriangle(P1, P2, P3), Colour(255, 255, 255), depthBuffer);
	}
}


void draw(DrawingWindow &window) {
	window.clearPixels();
	std::vector<ModelTriangle> faces = loadOBJ("textured-cornell-box.obj");
	std::vector<ModelTriangle> sphere = loadOBJ("sphere.obj");

	faces.insert(faces.end(), sphere.begin(), sphere.end());

	if (mode == 1) {
		projectImageRayTrace(window, faces);
	}
	else if (mode == 2) projectImageRaster(window, faces);
	else if (mode == 3) wireframeMode(window, faces);



	
}

void update(DrawingWindow &window) {
	// Function for performing animation (shifting artifacts or moving the camera)
	
	if (counter < 31) {
		t += 0.05;
		if (t > 2 * 3.14) {
			t = 0;
		}
		camX = 1.5 * cos(t);
		camZ = 1.5 * sin(t);
	}
	else if (counter >= 31 && camCounter<32) {

		glm::vec3 startPos{ 0.0, 0.0, 1.5 };
		glm::vec3 endPos{ 0.3, -0.1, 0.9 };

		auto xPos = interpolateSingleFloats(0.0, 0.3, 32);
		auto yPos = interpolateSingleFloats(0.0, -0.1, 32);
		auto zPos = interpolateSingleFloats(1.5, 0.9, 32);

		camX = xPos[camCounter];
		camY = yPos[camCounter];
		camZ = zPos[camCounter];

		camCounter++;
	}
	else if (counter > 62 && camCounter2 < 45) {
		//switch to ray tracer, work towards centre
		//lookatc = false;
		mode = 1;
		glm::vec3 startPos{ 0.3, -0.1, 0.9 };
		glm::vec3 endPos{ 0.0, 0.1, 0.9 };

		auto xPos = interpolateSingleFloats(0.3, 0.0, 45);
		auto yPos = interpolateSingleFloats(-0.1, 0.1, 45);
		auto zPos = interpolateSingleFloats(0.9, 0.9, 45);

		camX = xPos[camCounter2];
		camY = yPos[camCounter2];
		camZ = zPos[camCounter2];
		
		camCounter2++;
	}

	else if (counter >= 107 && camCounter3 < 50) {
		//switch to ray tracer, work towards centre
		//lookatc = false;
		mode = 1;
		glm::vec3 startPos{ 0.0, 0.1, 0.9 };
		glm::vec3 endPos{ -0.3, 0.3, 0.9 };

		auto xPos = interpolateSingleFloats(0.0, -0.3, 50);
		auto yPos = interpolateSingleFloats(0.1, 0.3, 50);
		auto zPos = interpolateSingleFloats(0.9, 0.9, 50);

		camX = xPos[camCounter3];
		camY = yPos[camCounter3];
		camZ = zPos[camCounter3];

		camCounter3++;
	}

	else if (counter >= 157 && camCounter4 < 59) {
		//switch to ray tracer, work towards centre
		//lookatc = false;
		mode = 1;
		glm::vec3 startPos{ -0.3, 0.3, 0.9 };
		glm::vec3 endPos{ -0.4, 0.3, 0.05 };

		auto xPos1 = interpolateSingleFloats(-0.3, -0.4, 59);
		auto yPos1 = interpolateSingleFloats(0.3, 0.3, 59);
		auto zPos1 = interpolateSingleFloats(0.9, 0.05, 59);

		auto xangleI = interpolateSingleFloats(0.0, 0.5, 59);
		auto yangleI = interpolateSingleFloats(0.0, 1.4, 59);
		xAngle = xangleI[camCounter4];
		yAngle = yangleI[camCounter4];

		camX = xPos1[camCounter4];
		camY = yPos1[camCounter4];
		camZ = zPos1[camCounter4];

		
		camCounter4++;
	}

	else if (counter >= 190 && camCounter5 < 33) {
		//switch to ray tracer, work towards centre
		//lookatc = false;
		mode = 1;
		glm::vec3 startPos{ -0.4, 0.3, 0.05 };
		glm::vec3 endPos{ -0.1, 0.3, 0.05 };

		auto xPos = interpolateSingleFloats(-0.4, -0.1, 33);
		auto yPos = interpolateSingleFloats(0.3, 0.3, 33);
		auto zPos = interpolateSingleFloats(0.05, 0.05, 33);

		
		auto lightMovementZ = interpolateSingleFloats(0.5, 0.19, 33);

		lightZ = lightMovementZ[camCounter5];
		camX = xPos[camCounter5];
		camY = yPos[camCounter5];
		camZ = zPos[camCounter5];


		camCounter5++;
	}




	std::cout << " camX: " << camX << " camY: " << camY << " camZ: " << camZ << std::endl;
	auto name = "img" + std::to_string(counter) + ".bmp";
	counter++;
	window.saveBMP(name);
	std::cout << " camX: " << camX << " camY: " << camY << " camZ: " << camZ << " xangle: "<<xAngle<<" yangle: "<<yAngle<< std::endl;
	std::cout << " C: " << camCounter << " C2: " << camCounter2 << " C3: " << camCounter3 << " C4: " << camCounter4 << " C5: " << camCounter5 << " CM: " << counter << std::endl;
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT){ 
			std::cout << "LEFT" << std::endl;
			camX = camX - 0.1;
		}
		else if (event.key.keysym.sym == SDLK_RIGHT){
			 std::cout << "RIGHT" << std::endl;
			 camX = camX + 0.1;
		}
		else if (event.key.keysym.sym == SDLK_UP){
			 std::cout << "UP" << std::endl;
			 camY = camY + 0.1;
		}
		else if (event.key.keysym.sym == SDLK_s){
			 std::cout << "rotate along x" << std::endl;
			 xAngle = xAngle + 0.1;
		}
		else if (event.key.keysym.sym == SDLK_DOWN){
			 std::cout << "DOWN" << std::endl;
			 camY = camY - 0.1;
		}
		else if (event.key.keysym.sym == SDLK_w){
			 std::cout << "rotate along x" << std::endl;
			 xAngle = xAngle - 0.1;
		}

		
		else if (event.key.keysym.sym == SDLK_e){
			 std::cout << "rotate along y" << std::endl;
			 yAngle = yAngle - 0.1;
		}
		else if (event.key.keysym.sym == SDLK_q){
			 std::cout << "rotate along y" << std::endl;
			 yAngle = yAngle + 0.1;
		}

		else if (event.key.keysym.sym == SDLK_1) {
			
			std::cout << "raytracing mode" << std::endl;
			mode = 1;
		
		}
		else if (event.key.keysym.sym == SDLK_2) {

			std::cout << "raster mode" << std::endl;
			mode = 2;

		}
		else if (event.key.keysym.sym == SDLK_3) {

			std::cout << "wireframe mode" << std::endl;
			mode = 3;

		}

		
	} 
	else if (event.type == SDL_MOUSEBUTTONDOWN) window.savePPM("output.ppm");

	else if (event.type == SDL_MOUSEWHEEL){
		if(event.wheel.y > 0){
			std::cout << "forward scroll" << std::endl;
			camZ = camZ - 0.1;
		}

		else if(event.wheel.y < 0){
			std::cout << "backward scroll" << std::endl;
			camZ = camZ + 0.1;
		}
	}
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
	texture = TextureMap("texture.ppm");
	
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

