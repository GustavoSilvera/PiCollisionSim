#include "cinder/app/AppNative.h"
#include "cinder/gl/gl.h"

using namespace ci;
using namespace ci::app;
using namespace std;

class PIhysicsApp : public AppNative {
  public:
	void setup();
	void mouseDown( MouseEvent event );	
	void update();
	void draw();
};

void PIhysicsApp::setup()
{
}

void PIhysicsApp::mouseDown( MouseEvent event )
{
}

void PIhysicsApp::update()
{
}

void PIhysicsApp::draw()
{
	// clear out the window with black
	gl::clear( Color( 0, 0, 0 ) ); 
}

CINDER_APP_NATIVE( PIhysicsApp, RendererGl )
