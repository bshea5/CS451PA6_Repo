#include "raytracer.h"

//some helper functions
inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; }

inline int toInt(double x){ return int( clamp(x)*255 + .5); }


//constructor
RayTracer::RayTracer(list<model>& models) : m_models(models)
{	
	//must have window and camera set up first 
	camera_pos = gli::getCameraPos();
	screenWidth = glutGet(GLUT_WINDOW_WIDTH);
	screenHeight = glutGet(GLUT_WINDOW_HEIGHT);

	glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
	glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
	glGetIntegerv(GL_VIEWPORT, viewport);
}


//render an image with "nray_per_pixel" rays created per pixel
void RayTracer::render(unsigned int img_w, unsigned int img_h, unsigned int nray_per_pixel)
{
	//create a black image of size img_w X img_h
	m_img = vector< vector< Vector3d > >(img_h, vector< Vector3d >(img_w, Vector3d(0,0,0) ));

	//generate rays
	for (int y = 0; y < img_h; y++)
	{

		for (int x = 0; x < img_w; x++)
		{
			Vector3d color;
			for (int n = 0; n < nray_per_pixel; n++)
			{
				Ray r=create_a_random_ray(x,y);
				pair<model *, triangle *> X = intersect(r);

				//determine the color of this ray 
				Vector3d rc = raycolor(*X.first, X.second, r);
				color = rc + color;
			}

			m_img[y][x] = color / nray_per_pixel;
		}//end of x

	}//end of y
}

// render an image with "nray_per_pixel" rays created per pixel
// create this ray randomly in the given pixel (x,y)
Ray RayTracer::create_a_random_ray(unsigned int x, unsigned int y)
{
	Ray r;
	Vector3d D;	//distance from camera to origin of the ray on the near plane

	GLdouble winX = x + (GLdouble)(rand() % 10) / 10.0; //get random decimal value [0,1] for offset
	GLdouble winY = y + (GLdouble)(rand() % 10) / 10.0;
	GLdouble objX, objY, objZ;

	gluUnProject(	winX,
					winY,
				 	0,				//get the near plane
				 	modelMatrix,	//model
				 	projMatrix, 	//projection
				 	viewport,		//viewport
				 	&objX,
				 	&objY,
				 	&objZ	);

	r.o = Point3d(objX, objY, objZ);	//set as ray origin

	// origin - camera position to get ray vector
	r.v = Vector3d(r.o[0], r.o[1], r.o[2]).normalize() - 
		  Vector3d(camera_pos[0], camera_pos[1], camera_pos[2]).normalize();

		  //check if rays are generating
	// int static count = 0;
	// count++;
	// if (count > 462300)
	// 	std::cout << "ray vector: " << r.v << std::endl;
	// if ( cound % 1000 == 0)
	//	draw line

	//TODO: implement this
	//hint: see slides on generating rays for perspective views
	// (x / number of pixels) * screenWidth + random[0,1]

	return r;
}

//returns a model and triangle that intersect ray r
//return pair<NULL,NULL> if no intersection is found
pair<model* , triangle* > RayTracer::intersect(Ray r)
{
	for (list<model>::iterator i = m_models.begin(); i != m_models.end(); i++)
	{
		triangle * t = intersect(*i, r);
		if (t != NULL) return make_pair(&(*i),t);
	}
	return make_pair((model *)NULL, (triangle *)NULL);
}

//returns a triangle in model m that intersects ray r
//return NULL if no intersection is found
triangle* RayTracer::intersect(model& m, Ray r)
{
	for (int i = 0; i < m.t_size; i++)
	{
		Point3d x;
		if ( intersect(m, &m.tris[i], r, x) )
			return &m.tris[i];
	}

	return NULL;
}

// determine if there is an intersection between triangle t and ray r
// return true if there is an intersection and store the intersection in x
// return false otherwise and x is undefined in this case
bool RayTracer::intersect(model& m, triangle* t, Ray r, Point3d& x)
{
	//TODO: implement this
	return false;
}

//
// determine the color of ray r, by analizing the intersection between t and r 
// 
Vector3d RayTracer::raycolor(model& m, triangle* t, const Ray& r)
{
	Vector3d color;

	//TODO: implement this

	//
	//hint: (1) what is the normal direction at the intersection between t and r
	//      (2) what is the color of this triangle? 
	//      (3) where is the light?
	//      (4) is the intersection in shadow? (call inshadow(const Point3d& p))

	return color;
}

//check if a point p is in shadow
bool RayTracer::inshadow(const Point3d& p)
{
	//TODO: implement this
	return false;
}

//save rendered image to file
bool RayTracer::save2file(const string& filename)
{
	FILE *f = fopen(filename.c_str(), "w");         // Write image to PPM file.

	int h = m_img.size();
	if (h == 0) return true; //nothing to save...

	int w = m_img.front().size();
	
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);

	for (int i = 0; i < h; i++) 
		for (int j = 0; j < w; j++) 
			fprintf(f, "%d %d %d ", toInt(m_img[i][j][0]), toInt(m_img[i][j][1]), toInt(m_img[i][j][2]));

	fclose(f);
}


