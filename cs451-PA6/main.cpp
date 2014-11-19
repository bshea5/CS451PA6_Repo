//------------------------------------------------------------------------------
//  Copyright 2007-2014 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------


#include "main.h"
#include "raytracer.h"

//
//
//
//
//   The MAIN function
//
//
//
//

int main(int argc, char ** argv)
{
    //
    if(!parseArg(argc,argv)){
        printUsage(argv[0]);
        return 1;
    }

    //read models
	initialize(rt_filename);
    //



    /////////////////////////////////////////////////////////////////
	//
    // setup glut/gli
	//
	/////////////////////////////////////////////////////////////////
	unsigned int image_w = 680, image_h = 680;
    glutInit( &argc, argv );
    glutInitDisplayMode( GLUT_RGB|GLUT_DOUBLE|GLUT_DEPTH );
	glutInitWindowSize(image_w, image_h);
    glutInitWindowPosition( 50, 50 );
	string title = "CS451 PA6 Ray Casting: ";
	glutCreateWindow(title.c_str());

    InitGL();
    gli::gliInit();
    gli::gliDisplayFunc(Display);
    glutReshapeFunc(Reshape);
    glutKeyboardFunc(Keyboard);

    //set camera position
    gli::setCameraPosZ(R*2.1f);

	/////////////////////////////////////////////////////////////////
	//
	// Ray trace rendering here
	//
	/////////////////////////////////////////////////////////////////

	unsigned int n_ray = 1;
	RayTracer rt(models);
	rt.render(image_w, image_h, n_ray);
	string imagename = "CS451_PA6_Ray_Casting";
	imagename = imagename +"_"+ to_string(time(NULL)) + "_.ppm";
	cout << "- Saving image (" << imagename << ")" << endl;
	rt.save2file(imagename);

    /////////////////////////////////////////////////////////////
    gli::gliMainLoop();

    return 0;
}

