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
//my own headers

//declaration for the main things, simulation and whatnot
using namespace ci;
using namespace ci::app;
using namespace std;

const double ppm = 100;//scale of how many pixels are in 1 SI Meter
const double refY = 10;//8m down from top of page
int numClacks = 0;
Font mFont;//custom font for optimized drawing
gl::TextureFontRef mTextureFont;//custom opengl::texture
static double kFreq = 60;//refresh rate
static double sqr(double x) { return x * x; }
double oldKfreq = 60;
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
class quadratic {
public:
	quadratic(double aparam, double bparam, double cparam) : a(aparam), b(bparam), c(cparam) {}
	double a, b, c;
	double solveAdd() {
		double ans = 0;
		try {
			ans = ((-b + sqrt(sqr(b) - 4 * a*c)) / (2 * a));
		}
		catch (const std::exception &e) {
			return 0;
		}
		return(ans);
	}
	double solveSub() {
		double ans = 0;
		try {
			ans = ((-b - sqrt(sqr(b) - 4 * a*c)) / (2 * a));
		}
		catch (const std::exception &e) {
			return 0;
		}
		return(ans);
	}
	vec3 solve() {
		vec3 ans;
		ans.X = solveAdd();
		ans.Y = solveSub();
		ans.Z = 0;
		return ans;
	}
};
class wall {
public:
	wall(vec3 p, vec3 dimX, vec3 dimY) : pos(p), dimensX(dimX), dimensY(dimY) {}
	vec3 pos;
	vec3 dimensX, dimensY;
	Area drawArea() {
		double L = dimensX.X;
		double R = dimensX.Y;
		double U = dimensY.X;
		double D = dimensY.Y;
		return (Area(ppm*(Vec2f(pos.X - L, (pos.Y - U) + refY)), ppm*(Vec2f(pos.X + R, (pos.Y + D) + refY))));
	}
};
wall w = wall(vec3(1, 0), vec3(0.2, 0.0), vec3(6.0, 0.0));
wall f = wall(vec3(1, 0), vec3(0.2, 15.0), vec3(0.0, 0.2));

class square {//ALL SI UNITS
public:
	square(float m, double v, vec3 p, double s) : mass(m), velocity(v), pos(p), size(s) {}
	vec3 pos;
	double mass;
	double size;
	double velocity;
	double momentum;
	double kinE;

	Area drawArea() {
		double r = size / 2;
		return (Area(ppm*(Vec2f(pos.X - r, (pos.Y) + refY)), ppm*(Vec2f(pos.X + r, (pos.Y - 2 * r) + refY))));
	}
	bool update(square other) {
		if (kFreq != 0) pos.X -= velocity / kFreq;
		else pos.X += 0;
		//pos.Y += velocity / 60;//NO VELOCITY IN Y UNIMPORTANNT
		kinE = 0.5 * (mass)* sqr(velocity);//k=1/2mv^2
		momentum = mass * velocity;
		return collideWall(w);
	}

	bool collideWall(class wall w) {//perfectly elastic collision
		double r = size / 2;
		if (pos.X - r < w.pos.X) {//checks both position and velocity direction
								  //kinematics time!!
			double timeA, timeB;//times within collision
			double oldVelocity = velocity;//v0 of current
			double distance = ((pos.X - r) - w.pos.X);//distance between edges of each squares;
			timeA = distance / oldVelocity;
			//NON ABS() BC REFRESH IS AFTER THE COLLISION PASSES THROUGH
			timeB = 1 / kFreq - fabs(timeA);//rest of the time
			velocity *= -1;//full elastic collision, 100% perfect
			pos.X = pos.X - (oldVelocity * timeA + velocity * timeB);//updates position of current block
			return true;
		}
		return false;
	}
 	bool collideSquare(class square *s) {//perfectly elastic collision
		double r = size / 2;
		if (fabs(pos.X - s->pos.X) <= r + (s->size / 2)) {//use kinematics
			double mr = (s->mass / mass);//mass ration m2/m1
			double a = sqr(mr) + mr;
			double b = -2 * (sqr(mr)*s->velocity + mr * velocity);
			double c = sqr(mr)*sqr(s->velocity) + 2 * velocity*s->velocity*mr - mr * sqr(s->velocity);
			quadratic newVel2 = quadratic(a, b, c);
			vec3 velZeros = newVel2.solve();
			//kinematics time!!
			double timeA, timeB;//times within collision
			double oldVelocity = velocity;//v0 of current
			double oldOtherVelocity = s->velocity;//v0 of other block
			double distance = ((pos.X + r) - (s->pos.X - s->size / 2));//distance between edges of each squares;
			timeA = distance / (oldOtherVelocity + oldVelocity);
			timeB = 1 / kFreq - fabs(timeA);//rest of the time
			if (fabs(s->velocity - velZeros.X) > 0.0001) {
				s->velocity = velZeros.X;//vF = other block
			}
			else s->velocity = velZeros.Y;
			velocity = velocity + mr * (oldOtherVelocity - s->velocity);//vF of current block
			if((s->mass == 100 || s->mass == 10000) && numClacks == 0){
				timeA = timeB;
				timeB = 1 / kFreq - fabs(timeA);
			}
			pos.X = pos.X - (oldVelocity * timeA + velocity * timeB);//updates position of current block
			s->pos.X = s->pos.X - (oldOtherVelocity * timeA + s->velocity * timeB);//updates position of new block
																				   //pos.X = s->pos.X - s->size / 2 - r;//resets position to be just before touching the other block
			return true;
		}
		return false;
	}
};
square sqr1 = square(1, 0, vec3(4, 0), 1);
square sqr2 = square(100, 1, vec3(7.5, 0), 1.5);

//begin
int tX = 1200;
bool pause = false;
class PIhysicsApp : public AppNative {
public:
	PIhysicsApp() {}
	void prepareSettings(Settings *settings);
	void setup();
	void mouseDown(MouseEvent event);
	void mouseUp(MouseEvent event);
	void mouseMove(MouseEvent event);
	void keyDown(KeyEvent event);
	void keyUp(KeyEvent event);
	void update();
	void textDraw(square s);
	static void drawFontText(float text, vec3 pos);
	audio::VoiceRef clackSound;
	void clack();
	struct text {
		string s;
	};
	void draw();
private:
	// Change screen resolution
	int mScreenWidth, mScreenHeight;
	float initWidth, initHeight;
	void getScreenResolution(int& width, int& height);
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
	if (event.getCode() == KeyEvent::KEY_UP) sqr2.mass *= 100;
	if (event.getCode() == KeyEvent::KEY_DOWN) if (sqr2.mass >= 100) sqr2.mass /= 100;
	if (event.getCode() == KeyEvent::KEY_RIGHT) kFreq += 60;
	if (event.getCode() == KeyEvent::KEY_LEFT) if (kFreq >= 60) kFreq -= 60;
	if (event.getChar() == 'R' || event.getChar() == 'r') {
		sqr1 = square(1, 0, vec3(4, 0), 1);
		sqr2 = square(100, 1, vec3(7.5, 0), 1.5);
		numClacks = 0;
		kFreq = 60;
	}
	if (event.getCode() == KeyEvent::KEY_SPACE) {
		pause = !pause;
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
	numClacks++;
}
void PIhysicsApp::update() {
	if (sqr1.update(sqr2)) clack();
	if (sqr2.update(sqr1)) clack();
	if (sqr1.collideSquare(&sqr2)) clack();
	if (!pause) {
		kFreq = oldKfreq;
		if (fabs(sqr1.velocity) > 5 && 1.5 * sqr(sqr1.velocity) < 180000) {
			kFreq = 1.5 * sqr(sqr1.velocity);
		}
		oldKfreq = kFreq;
	}
	else kFreq = 0;
}
//drawing the text used for debugging or just other details
void PIhysicsApp::textDraw(square s) {//function for drawing the buttons 
									  //(	WARNING: RESOURCE HOG!!!!!!!!!!!)
	struct text {
		string s;
		double f;
	};
	text t[] = {

		{ "mass", (s.mass) },
	{ "vel", (s.velocity) }
	//{ "momt", (s.momentum) },
	//{ "kinE", (s.kinE) }
	//velocity and acceleration measured with drawDials
	};
	int i = 0;
	for (text& ti : t) {
		int tY = (i + 4*s.size);//increment y position for each button based off index
		gl::drawString(ti.s, ppm*Vec2f(s.pos.X - s.size / 2, (s.pos.Y + refY) - tY), Color(1, 1, 1), Font("Times", 45));
		drawFontText(ti.f, vec3(s.pos.X + s.size / 2, (s.pos.Y + refY) - tY).times(ppm));
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
	gl::drawSolidRect(sqr1.drawArea());
	gl::drawSolidRect(sqr2.drawArea());
	gl::drawSolidRect(w.drawArea());
	gl::drawSolidRect(f.drawArea());
	textDraw(sqr1);
	textDraw(sqr2);

	gl::color(1, 1, 1);

	gl::drawString("FPS: ", Vec2f(getWindowWidth() - 250, 10), Color(0, 1, 0), Font("Arial", 45));
	drawFontText(getAverageFps(), vec3(getWindowWidth() - 130, 10));

	gl::drawString("Clacks: ", Vec2f(getWindowWidth() - 250, 100), Color(0, 1, 0), Font("Times", 45));
	drawFontText(numClacks, vec3(getWindowWidth() - 130, 100));

	gl::drawString("Refresh: ", Vec2f(getWindowWidth() - 250, 190), Color(0, 1, 0), Font("Times", 45));
	drawFontText(kFreq, vec3(getWindowWidth() - 130, 190));


}
//awesomesause
CINDER_APP_NATIVE(PIhysicsApp, RendererGl)