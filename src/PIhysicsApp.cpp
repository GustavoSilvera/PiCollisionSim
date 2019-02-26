#include "cinder/app/AppNative.h"
#include "cinder/gl/gl.h"
#include "cinder/ImageIo.h"
#include "cinder/gl/Texture.h"
#include "cinder/Vector.h"
#include "cinder/Text.h"
#include "cinder/Font.h"
#include "cinder/audio/Voice.h"
#include "cinder/gl/TextureFont.h"
#include "cinder/Utilities.h"
#include <ostream>
#include <fstream>
#include <vector>
#include <filesystem>
#include <utility>
//my own headers

//declaration for the main things, simulation and whatnot
using namespace ci;
using namespace ci::app;
using namespace std;

const double ppm = 100;//scale of how many pixels are in 1 SI Meter

Font mFont;//custom font for optimized drawing
gl::TextureFontRef mTextureFont;//custom opengl::texture
static double kFreq = 60;//refresh rate
static double sqr(double x) { return x * x; }

class vec3 {
public:
	vec3(double x = 0.0, double y = 0.0, double z = 0.0) : X(x), Y(y), Z(z) {}
	double X, Y, Z;
	vec3 times(double f) { return vec3(X*f, Y*f, Z*f); }
	double distance(vec3 v) { return sqrt(sqr(v.X - X) + sqr(v.Y - Y)); }
	double distanceV3(vec3 v) { return sqrt(sqr(v.X - X) + sqr(v.Y - Y) + sqr(v.Z - Z)); }
	vec3 operator+(vec3 v) { return vec3(X + v.X, Y + v.Y, Z + v.Z); }
	vec3 operator*(double f) { return vec3(f*X, f*Y, f*Z); }
};

class Polynomial {
public:
	Polynomial()  {}
	explicit Polynomial(double c0) {
		coeffs.push_back(c0);
	}
	Polynomial(double c1, double c0) {
		coeffs.push_back(c0);
		coeffs.push_back(c1);
	}
	Polynomial(double c2, double c1, double c0) {
		coeffs.push_back(c0);
		coeffs.push_back(c1);
		coeffs.push_back(c2);
	}
	Polynomial(double c, int i) { // c x ^ i
		coeffs.resize(i);
		coeffs.push_back(c);
	}
	
	double eval(double x) const {
		double result = 0.0;
		double xi = 1.0;
		for (int i = 0; i < coeffs.size(); ++i) {
			result += coeffs[i] * xi;
			xi *= x;
		}
		return result;
	}

	int size() const { return coeffs.size(); }
	double operator[]  (int index) const {
		return coeffs[index];
	}

	Polynomial & operator+=(const Polynomial &p) {
		coeffs.resize(std::max(size(), p.size()));
		for	(int i = 0; i < p.size(); ++i) {
			coeffs[i] += p[i];
		}
		return *this;
	}
	Polynomial & operator-=(const Polynomial &p) {
		coeffs.resize(std::max(size(), p.size()));
		for (int i = 0; i < p.size(); ++i) {
			coeffs[i] -= p[i];
		}
		return *this;
	}
	Polynomial & operator*=(const double c) {
		for (int i = 0; i < coeffs.size(); ++i) {
			coeffs[i] *= c;
		}
		return *this;
	}
private:
	// coeff[i] is term for x^[i]
	std::vector<double> coeffs;
};

std::pair<double, double> Quadratic(double a, double b, double c) {
	double bsm4ac = sqr(b) - 4 *a*c;
	if (bsm4ac >= 0) {
		return { (-b + sqrt(bsm4ac)) / (2 * a),
		 		 (-b - sqrt(bsm4ac)) / (2 * a) };
	}
	return { std::numeric_limits<double>::quiet_NaN(),
		std::numeric_limits<double>::quiet_NaN() };
}

Polynomial operator*(const Polynomial &p0, const Polynomial &p1) {
	Polynomial result;
	for (int i = 0; i < p0.size(); ++i) {
		double pi = p0[i];
		for (int j = 0; j < p1.size(); ++j) {
			double pj = p1[j];
			result += Polynomial(pi * pj, i+j);
		}
	}
	return result;
}

Polynomial operator*(double c, const Polynomial &p) {
	Polynomial result = p;
	result *= c;
	return result;
}

Polynomial operator-(const Polynomial &p0, const Polynomial &p1) {
	Polynomial result = p0;
	result -= p1;
	return result;
}

Polynomial sqr(const Polynomial &p) {
	return p * p;
}

std::pair<double, double> Quadratic(const Polynomial &a, const Polynomial &b, const Polynomial &c, double mr) {
	double bsm4ac = (sqr(b) - 4 * a * c).eval(mr);
	if (bsm4ac >= 0) {
		double a2 = (2 * a).eval(mr);
		return { -b.eval(mr) / a2 + sqrt(bsm4ac) / a2,
			     -b.eval(mr) / a2 - sqrt(bsm4ac) / a2 };
	}
	return { std::numeric_limits<double>::quiet_NaN(),
		std::numeric_limits<double>::quiet_NaN() };
}

/////////////////////scene.h
using Coordinates = std::pair<Vec2f, Vec2f>;
class Box_properties {
public:
	double velocity;
	double mass;
	double kinE;
	double momentum;
	double posX;
};
class Box {
public:
	Box(double distance_to_wall, double size):
		distance_to_wall{ distance_to_wall },
		size{ size } {	}
	Coordinates Coords() const {
		return { Vec2f{ static_cast<float>(distance_to_wall),0.0f },//leftmost edge	
			Vec2f{ static_cast<float>(distance_to_wall + size), static_cast<float>(size) } };//rightmost edge
	}
	double DistanceToWall() const {
		return distance_to_wall;
	}
	double Size() const { return size; }
private:
	double distance_to_wall;
	double size;
};

void EnergyTransfer(double& small_speed, double small_mass, double& large_speed, double large_mass) { 
	//std::swap(*small_speed, *large_speed); // full energy transfer
	const double mr = (large_mass / small_mass);//mass ration m2/m1
	//divide everything by mr
    //	double a = (sqr(mr) + mr);
	Polynomial a{ 1, 1, 0 };  // mr^2 + mr + 0
	//	double b = -2 * (sqr(mr)*large_speed + mr * small_speed);
	Polynomial b{ -2 * large_speed, -2 * small_speed, 0 }; // -2(mr^2 * ls + mr * ss)
	//	double c = sqr(mr)*sqr(large_speed) + 2 * small_speed*large_speed*mr - mr * sqr(large_speed);
	Polynomial c{ sqr(large_speed), 2 * small_speed * large_speed - sqr(large_speed), 0 };
	//DOSENT WORK FOR 1e+8
	/*double a = 1 + mr;
	double b = -2 * (mr * large_speed + small_speed);
	double c = mr * sqr(large_speed) + 2 * small_speed * large_speed - sqr(large_speed);*/
	std::pair<double, double> velZeros= Quadratic(a, b, c, mr);
	double oldLargeVelocity = large_speed;//v0 of other block
	if (std::abs(large_speed - velZeros.first) > std::abs(large_speed - velZeros.second)) {//takes the one further away
		large_speed = velZeros.first;//vF = other block
	}
	else large_speed = velZeros.second;
	small_speed = small_speed + mr * oldLargeVelocity - mr * large_speed;//vF of current block
}

class State {
public:
	State(double small_x, double large_x, double time) ://snapshot of time
		small{ small_x, 1.0},
		large{ large_x, 1.5},
		time{ time }
	{}
	Box Small() const { return small; }
	Box Large() const { return large; }
	double Time() const { return time; }

	State NextCollision(double small_speed, double small_mass, double large_speed, double large_mass) const {
		if (small.DistanceToWall() == 0.0) {
			// Bounce small off the wall
			small_speed = -small_speed;//perfectly elastic collision
			//large_speed = large_speed;
		} else {
			// energy/momentum transfer
			EnergyTransfer(small_speed, small_mass, large_speed, large_mass);
		}
		double time_to_collision = 60 * 60 * 24 * 365.25 * 100; // 100 years (forever)
		if (small_speed > large_speed ) {
			// Going towards each other.
			double distance = large.DistanceToWall() - small.DistanceToWall() - small.Size();
			time_to_collision = distance / (small_speed - large_speed);
		}
		if (small_speed < 0.0) {
			double time_to_wall = small.DistanceToWall() / - small_speed;
			if (time_to_wall < time_to_collision) {
				return State{ 0.0,
					large.DistanceToWall() + time_to_wall * large_speed,
					time + time_to_wall };
			}
		}
		return State{	small.DistanceToWall() + time_to_collision * small_speed,
						large.DistanceToWall() + time_to_collision * large_speed,
						time + time_to_collision 
		};
	}
private:
	Box small, large;
	double time;
};

double Interpolate(double a, double b, double ratio) {
	return a + (b - a)*ratio;
}

State Interpolate(const State &from, const State &to, double time) {
	double ratio = (time - from.Time()) / (to.Time() - from.Time());
	return State{ Interpolate(from.Small().DistanceToWall(), to.Small().DistanceToWall(), ratio),
		Interpolate(from.Large().DistanceToWall(), to.Large().DistanceToWall(), ratio),
		time };
}

class Scene {
public:
	Scene() :
		current{1.0, 5.0, 0.0},
		previous_clack{ current },
		next_clack{ 1.0, 2.0, 3.0 }
	{	
		UpdateBoxProperties();
	}
	// TODO: remove small_mass/large_mass and make getters for Box_properties
	double small_mass = 1, large_mass = 1;
	Box_properties small_prop, large_prop;
	bool Update(double new_time) {
		bool clacked = false;
		while (next_clack.Time() < new_time) {
			clacked = true;
			double time_delta = next_clack.Time() - previous_clack.Time();
			double small_speed = (next_clack.Small().DistanceToWall() - previous_clack.Small().DistanceToWall()) / time_delta;
			double large_speed = (next_clack.Large().DistanceToWall() - previous_clack.Large().DistanceToWall()) / time_delta;
			previous_clack = next_clack;
			next_clack = next_clack.NextCollision(small_speed, small_mass, large_speed, large_mass);
			++num_clacks;
		}
		if (clacked) UpdateBoxProperties();
		current = Interpolate(previous_clack, next_clack, new_time);
		//update box properties (for drawing)
		// TODO: REMOVE THESE HORRIBLE HACKS!
		small_prop.mass = small_mass;
		large_prop.mass = large_mass;
		small_prop.posX = current.Small().DistanceToWall();
		large_prop.posX = current.Large().DistanceToWall();
		return clacked;
	}
	Coordinates SmallBoxCoords() const {
		return current.Small().Coords();
	}
	Coordinates LargeBoxCoords() const {
		return current.Large().Coords();
	}
	Box SmallBox() const {
		return current.Small();
	}
	Box LargeBox() const {
		return current.Large();
	}
	int TotalClacks() const {
		return num_clacks;
	}
private:
	void UpdateBoxProperties() {
		double time_delta = next_clack.Time() - previous_clack.Time();
		double small_speed = (next_clack.Small().DistanceToWall() - previous_clack.Small().DistanceToWall()) / time_delta;
		double large_speed = (next_clack.Large().DistanceToWall() - previous_clack.Large().DistanceToWall()) / time_delta;
		small_prop.velocity = small_speed;
		large_prop.velocity = large_speed;
	}
	State current, previous_clack, next_clack;
	int num_clacks = 0;
};
/////////////////////scene.h

//begin
int tX = 1200;
class PIhysicsApp : public AppNative {
public:
	PIhysicsApp()  {}

	void prepareSettings(Settings *settings);
	void setup();
	void mouseDown(MouseEvent event);
	void mouseUp(MouseEvent event);
	void mouseMove(MouseEvent event);
	void keyDown(KeyEvent event);
	void keyUp(KeyEvent event);
	void update() override {
		if(kFreq != 0) current_time += 1.0 / kFreq;
		if (scene.Update(current_time)) clack();
	}

	void textDraw(class Box_properties s);
	static void drawFontText(float text, vec3 pos);
	audio::VoiceRef clackSound;
	void clack();
	//	struct text {
	//		string s;
	//	};
	void draw() override;
	Area MapArea(Coordinates box) {
		const Vec2f ref = { 2, 10 }; // 2m right, 8 down

		const double kHeight = 1.0;
		return Area(ppm*(ref + Vec2f(box.first.x, -box.first.y)),
			ppm*(ref + Vec2f(box.second.x, -box.second.y)));
	}
private:
	// Change screen resolution
	int mScreenWidth, mScreenHeight;
	float initWidth, initHeight;
	void getScreenResolution(int& width, int& height);

	double current_time = 0.0;
	class Scene scene;
};

void PIhysicsApp::prepareSettings(Settings *settings) {
	//! setup our window
	settings->setTitle("PIhysics");
	settings->setFrameRate(60);//60fps
	gl::enableVerticalSync();//vsync
	getScreenResolution(mScreenWidth, mScreenHeight);//getss resolution relative to monitor
	settings->setWindowPos(mScreenWidth / 6, mScreenHeight / 6);
	int aspectRatio = mScreenWidth / 7;//using 4/7ths of monitor resolution
	initWidth = aspectRatio * 4;
	initHeight = aspectRatio * 3;
	//winConst = 0.000260417*(mScreenWidth)+0.2;
	settings->setWindowSize(initWidth, initHeight);//maintains 4:3 aspect ratio
}
void PIhysicsApp::getScreenResolution(int& width, int& height) {
	// Get a handle to the desktop window
	const HWND hDesktop = GetDesktopWindow();

	// Get the size of the screen into a rectangle
	RECT rDesktop;
	GetWindowRect(hDesktop, &rDesktop);

	// The top left corner will have coordinates (0, 0)
	// and the bottom right corner will have coordinates
	// (width, height)
	width = rDesktop.right;
	height = rDesktop.bottom;

}
//initial setup for all the variables and stuff
void PIhysicsApp::setup() {
	srand(time(NULL));//seeds random number generator
	mFont = Font("Times", 45);//fixed custom font
	mTextureFont = gl::TextureFont::create(mFont);
	audio::SourceFileRef sourceFile = audio::load(app::loadAsset("clack.mp3"));
	clackSound = audio::Voice::create(sourceFile);
	//clack2 = audio::Voice::create(sourceFile);
}
//setup for on screen buttons
//when mouse is clicked
void PIhysicsApp::mouseDown(MouseEvent event) {
	if (event.isLeft()) {

	}
}
//when mouse is released
void PIhysicsApp::mouseUp(MouseEvent event) {

}
//when mouse is moved (used with joystick)
void PIhysicsApp::mouseMove(MouseEvent event) {

}
//what to do when keyboard key is pressed
void PIhysicsApp::keyDown(KeyEvent event) {
	//base
	/*if (event.getCode() == KeyEvent::KEY_UP || event.getChar() == 'w' || event.getChar() == 'W')	v.r[0].ctrl.KeyUp = true;
	if (event.getCode() == KeyEvent::KEY_DOWN || event.getChar() == 's' || event.getChar() == 'S')	v.r[0].ctrl.KeyDown = true;
	if (event.getCode() == KeyEvent::KEY_LEFT || event.getChar() == 'a' || event.getChar() == 'A')	v.r[0].ctrl.KeyLeft = true;
	if (event.getCode() == KeyEvent::KEY_RIGHT || event.getChar() == 'd' || event.getChar() == 'D') v.r[0].ctrl.KeyRight = true;*/
}
//what to do when keyboard key is released
void PIhysicsApp::keyUp(KeyEvent event) {
	//base
	if (event.getCode() == KeyEvent::KEY_UP) scene.large_mass *= 100;
	if (event.getCode() == KeyEvent::KEY_DOWN) if (scene.large_mass >= 100) scene.large_mass /= 100;
	if (event.getChar() == 'R' || event.getChar() == 'r') {
		scene = Scene();//reset
		current_time = 0;
	}
	if (event.getCode() == KeyEvent::KEY_SPACE) {
		if (kFreq == 60) kFreq = 0;
		else kFreq = 60;
	}

	/*if (event.getCode() == KeyEvent::KEY_DOWN || event.getChar() == 's' || event.getChar() == 'S') v.r[0].ctrl.KeyDown = false;
	if (event.getCode() == KeyEvent::KEY_UP || event.getChar() == 'w' || event.getChar() == 'W') v.r[0].ctrl.KeyUp = false;
	if (event.getCode() == KeyEvent::KEY_RIGHT || event.getChar() == 'd' || event.getChar() == 'D') v.r[0].ctrl.KeyRight = false;
	if (event.getCode() == KeyEvent::KEY_LEFT || event.getChar() == 'a' || event.getChar() == 'A') v.r[0].ctrl.KeyLeft = false;
	*/
}
//overall application update function
void PIhysicsApp::clack() {
	clackSound->stop();
	clackSound->start();
}

//drawing the text used for debugging or just other details
void PIhysicsApp::textDraw(class Box_properties s) {//function for drawing the buttons 
										 //(	WARNING: RESOURCE HOG!!!!!!!!!!!)
	struct text {
		string s;
		double f;
	};
	text t[] = {
		{ "mass:", (s.mass) },
		{ "vel:", (s.velocity) }
	//{ "momt", (s.momentum) },
	//{ "kinE", (s.kinE) }
	};
	int i = 0;
	const int refY = 9, refX = 2.5;
	
	for (text& ti : t) {
		int tY = (i + 4);//increment y position for each button based off index
		//double posX = (s.Coords().first + s.Coords().second) / 2;//middle of two edges
		gl::drawString(ti.s, ppm*Vec2f(s.posX + refX, refY - tY), Color(1, 1, 1), Font("Times", 45));
		drawFontText(ti.f, vec3(s.posX + refX + 1, refY - tY).times(ppm));
		++i;
	}
	gl::color(1, 1, 1);
}

void PIhysicsApp::drawFontText(float text, vec3 pos) {
	std::stringstream dummyText;
	std::string PRINT;
	dummyText << text;
	dummyText >> PRINT;
	gl::color(1, 1, 1);
	mTextureFont->drawString(PRINT, Vec2f(pos.X, pos.Y + 20));
	gl::color(1, 1, 1);
}

void PIhysicsApp::draw() {
	gl::enableAlphaBlending();//good for transparent images
	gl::clear(Color(0, 0, 0));

	gl::color(1, 1, 1);
	gl::drawSolidRect(MapArea(scene.SmallBoxCoords()));
	gl::drawSolidRect(MapArea(scene.LargeBoxCoords()));
	gl::drawSolidRect(MapArea(Coordinates{{-0.1f, 0.0f }, { 0.0f, 5.0f }})); // Wall
	gl::drawSolidRect(MapArea(Coordinates{{-0.1f, 0.0f },{ 10.0f, -0.1f } })); // Floor
	textDraw(scene.small_prop);
	textDraw(scene.large_prop);

	gl::color(1, 1, 1);

	gl::drawString("FPS: ", Vec2f(getWindowWidth() - 250, 10), Color(0, 1, 0), Font("Arial", 45));
	drawFontText(getAverageFps(), vec3(getWindowWidth() - 130, 10));

	gl::drawString("Clacks: ", Vec2f(getWindowWidth() - 250, 100), Color(0, 1, 0), Font("Times", 45));
	drawFontText(scene.TotalClacks(), vec3(getWindowWidth() - 130, 100));

	gl::drawString("Refresh: ", Vec2f(getWindowWidth() - 250, 190), Color(0, 1, 0), Font("Times", 45));
	drawFontText(kFreq, vec3(getWindowWidth() - 130, 190));
}
//awesomesause
CINDER_APP_NATIVE(PIhysicsApp, RendererGl)
