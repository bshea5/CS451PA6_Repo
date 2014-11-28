#include "raytracer.h"
#include <algorithm>
#include <cfloat>

//NEW in PA7
#define MAX_RAY_DEPTH 5
#define SHADOW_SAMPLE_SIZE 10

//some helper functions
inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; }

inline int toInt(double x){ return int( clamp(x)*255 + .5); }

//global variable vector, :( , so dirty
Point3d current_P;

//constructor
RayTracer::RayTracer(list<object3D*>& models) : m_models(models)
{
	//find out where lights are
	m_lights.clear();
	for (list<object3D*>::iterator i = models.begin(); i != models.end(); i++)
	{
		object3D* obj = *i;
		if (obj->mat_emission.normsqr() != 0) m_lights.push_back(*i);
	}

	if (m_lights.empty())
	{
		cerr << "! Error: There no light defined. Make sure that some object has has non-zero emission." << endl;
		exit(0);
	}

	//must have window and camera set up first 
	camera_pos = gli::getCameraPos();
	screenWidth = glutGet(GLUT_WINDOW_WIDTH);
	screenHeight = glutGet(GLUT_WINDOW_HEIGHT);

	glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
	glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
	glGetIntegerv(GL_VIEWPORT, viewport);
}

inline void show_progress_bar(double ratio)
{
	// Show the load bar.
	static const int w = 50;
	int   c = (int)(ratio * w);

	cout << setw(3) << (int)(ratio * 100) << "% [";
	for (int x = 0; x<c; x++) cout << "=";
	for (int x = c; x<w; x++) cout << " ";
	cout << "]\r" << flush;
}

//render an image with "nray_per_pixel" rays created per pixel
void RayTracer::render(unsigned int img_w, unsigned int img_h, unsigned int nray_per_pixel)
{
	//create a black image of size img_w X img_h
	m_img = vector< vector< Vector3d > >(img_h, vector< Vector3d >(img_w, Vector3d(0,0,0) ));
	 
	//initialize the bar
	show_progress_bar(0);

	//generate rays
	for (unsigned int y = 0; y < img_h; y++)
	{

		for (unsigned int x = 0; x < img_w; x++)
		{
			Vector3d color;
			for (unsigned int n = 0; n < nray_per_pixel; n++)
			{
				Ray r=create_a_random_ray(x,y);
				Vector3d rc = raycolor(r, 0);
				color = rc + color;
			}

			m_img[y][x] = color / nray_per_pixel;
		}//end of x

		show_progress_bar(y*1.0 / img_h);

	}//end of y
	show_progress_bar(1.0f);
	cout << endl;
}

// render an image with "nray_per_pixel" rays created per pixel
// create this ray randomly in the given pixel (x,y)
Ray RayTracer::create_a_random_ray(unsigned int x, unsigned int y)
{
	Ray r;
	Vector3d n, f;	//near and far plane

	GLdouble pixelWidth = viewport[2] / m_img.front().size();
	GLdouble pixelHeight = viewport[3] / m_img.size(); 
	GLdouble winX = x + (GLdouble)(rand() % 10) / 10.0; //get random decimal value [0,1] for offset
	GLdouble winY = (viewport[3] - y) + (GLdouble)(rand() % 10) / 10.0;
	winX *= pixelWidth;
	winY *= pixelHeight;

	//get near
	gluUnProject(	winX,
					winY,
				 	0,				//get the near plane
				 	modelMatrix,	//model
				 	projMatrix, 	//projection
				 	viewport,		//viewport
				 	&n[0],
				 	&n[1],
				 	&n[2]	);

	//get far
	gluUnProject(	winX,
					winY,
				 	1,				//get the far plane
				 	modelMatrix,	
				 	projMatrix, 	
				 	viewport,		
				 	&f[0],
				 	&f[1],
				 	&f[2]		);

	//ray vector
	//far - near
	r.v = (f - n).normalize();

	//ray origin
	//add camera position to origin
	r.o = Point3d(	camera_pos[0] + n[0], 
					camera_pos[1] + n[1], 
					camera_pos[2] + n[2]	);

	//all_rays.push_back(r);
	return r;
}

//returns a model and triangle that intersect ray r
//return pair<NULL,NULL> if no intersection is found
pair<object3D *, triangle *> RayTracer::intersect(Ray r)
{
	double min_dist = FLT_MAX;
	triangle * closest_T = NULL;
	object3D * closest_M = NULL;

	for (list<object3D*>::iterator i = m_models.begin(); i != m_models.end(); i++)
	{
		object3D* obj = *i;

		if (dynamic_cast<model*>(obj)) //this object is a mesh
		{
			model* mesh = dynamic_cast<model*>(obj);
			triangle * t = closest_intersect(*mesh, r);
			if (t != NULL)
			{
				Point3d x;
				intersect(*mesh, t, r, x);
				double dist = (x - r.o).normsqr();
				if (dist < min_dist)
				{
					min_dist = dist;
					closest_T = t;
					closest_M = mesh;
				}
			}
		}//--------------------------------------------------------------------------
		else if (dynamic_cast<sphere*>(obj)) //this object is a sphere
		{
			sphere* ball = dynamic_cast<sphere*>(obj);
			Point3d x;
			if (intersect(*ball, r, x))
			{
				double dist = (x - r.o).normsqr();
				if (dist < min_dist)
				{
					min_dist = dist;
					closest_T = NULL;
					closest_M = ball;
				}
			}
		}//--------------------------------------------------------------------------
	}

	return make_pair(closest_M, closest_T);
}

//
//returns true if the sphere intersects ray r
//return false if no intersection is found
//x is the location of intersection if true is returned.
//
bool RayTracer::intersect(sphere& s, Ray r, Point3d& x)
{
	//
	// PA7 TODO: implement this
	// 

	return false;
}


//returns a triangle in model m that intersect ray r
//return NULL if no intersection is found
triangle * RayTracer::intersect(model& m, Ray r)
{
	for (unsigned int i = 0; i < m.t_size; i++)
	{
		Point3d x;
		if ( intersect(m, &m.tris[i], r, x) )
			return &m.tris[i];
	}

	return NULL;
}

//returns a triangle in model m that make closest intersection with ray r
//return NULL if no intersection is found
triangle * RayTracer::closest_intersect(model& m, Ray r)
{
	double min_dist = FLT_MAX;
	triangle * closest_T=NULL;
	for (unsigned int i = 0; i < m.t_size; i++)
	{
		Point3d x;
		if (intersect(m, &m.tris[i], r, x))
		{
			double dist = (x - r.o).normsqr();
			if (dist < min_dist)
			{
				closest_T = &m.tris[i];
				min_dist = dist;
			}
		}//end if
	}//end for i

	return closest_T;
}

// determine if there is an intersection between triangle t and ray r
// return true if there is an intersection and store the intersection in x
// return false otherwise and x is undefined in this case
bool RayTracer::intersect(model& m, triangle* tri, Ray r, Point3d& x)
{
	Vector3d e1, e2, h, s, q;
	Vector3d v0, v1, v2;	
	double a, f;
	double t, u, v;	//B

	//triangle vertices, save points as vectors for easy calculations
	v0 = Vector3d(m.vertices[tri->v[0]].p[0],m.vertices[tri->v[0]].p[1],m.vertices[tri->v[0]].p[2]);
	v1 = Vector3d(m.vertices[tri->v[1]].p[0],m.vertices[tri->v[1]].p[1],m.vertices[tri->v[1]].p[2]);
	v2 = Vector3d(m.vertices[tri->v[2]].p[0],m.vertices[tri->v[2]].p[1],m.vertices[tri->v[2]].p[2]);

	e1 = v1 - v0;
	e2 = v2 - v0;

	h = r.v % e2;	//cross product
	a = e1 * h;		//dot product

	if (a > -0.00001 && a < 0.00001)
		return(false);

	f = 1/a;
	s = r.o - v0;		//distance from vertex0 to ray origin
	u = f * (s * h);

	if (u < 0.0 || u > 1.0)
		return(false);

	q = s % e1;
	v = f * (r.v * q);

	if (v < 0.0 || u + v > 1.0)
		return(false);

	// at this stage we can compute t to find out where
	// the intersection point is on the line
	t = f * (e2 * q);

	if (t > 0.00001) // ray intersection
	{	
		x = r.o + (r.v * t);	//get point of intersection
		return(true);
	}

	else // this means that there is a line intersection
		 // but not a ray intersection
		 return (false);
}


//
// determine the color of ray r
// 
Vector3d RayTracer::raycolor(const Ray& r, int depth)
{
	if (depth >= MAX_RAY_DEPTH) return Vector3d(0,0,0);

	Vector3d color;
	pair<object3D *, triangle *> X = intersect(r);

	//determine the color of this ray 
	if (X.first != NULL && X.second != NULL) //this is mesh intersection
	{
		Vector3d rc = raycolor(*((model*)X.first), X.second, r,depth);
		color = rc + color;
	}
	else if (X.first != NULL) //this is spehere intersection
	{
		Vector3d rc = raycolor(*((sphere*)X.first), r, depth);
		color = rc + color;
	}

	return color;
}

//
// determine the color of ray r, by analizing the intersection between t and r 
// 
Vector3d RayTracer::raycolor(sphere& s, const Ray& r, int depth)
{
	Vector3d color;

	//PA7 TODO: determine the direct color
	//
	//hint: (1) what is the normal direction at the intersection between s and r
	//      (2) what is the color of this sphere? 
	//      (3) where is the light?
	//      (4) is the intersection in shadow? (call inshadow(const Point3d& p))

	Vector3d reflection_color;
	if (s.reflectiveness.normsqr() > 0)
	{
		//PA7 TODO: determine the reflection color
		//
		//      (1) generate relection ray 
		//      (2) determine the color of the ray (remember of increase depth by 1)
	}

	Vector3d refraction_color;
	if (s.transparency.normsqr() > 0)
	{
		//PA7 TODO: determine the refraction color
		//
		//      (1) generate refraction ray 
		//      (2) determine the color of the ray (remember of increase depth by 1)
	}

	//finally gather all the colors
	for (int i = 0; i < 3; i++)
	{
		color[i] = clamp(color[i] + s.reflectiveness[i] * reflection_color[i] + s.transparency[i] * refraction_color[i]);
	}

	return color;
}


//
// determine the color of ray r, by analizing the intersection between t and r 
// 
Vector3d RayTracer::raycolor(model& m, triangle* tri, const Ray& r, int depth)
{
	Vector3d color; 			//color of point to return
	Vector3d a, b, c;			//triangle points
	Vector3d v0, v1, v2;
	Vector3d nx, na, nb, nc;	//normal directions
	Vector3d diffuse, ambient, specular;
	Vector3d light_direction, light_reflected;
	Point3d light_position;
	double ac, dc, sc;			//coefficients
	double u, v, w; 			//Barycentric coords
	double d00, d01, d11, d20, d21, denom;

	//triangle points in vector form
	a = Vector3d(m.vertices[tri->v[0]].p[0],m.vertices[tri->v[0]].p[1],m.vertices[tri->v[0]].p[2]);
	b = Vector3d(m.vertices[tri->v[1]].p[0],m.vertices[tri->v[1]].p[1],m.vertices[tri->v[1]].p[2]);
	c = Vector3d(m.vertices[tri->v[2]].p[0],m.vertices[tri->v[2]].p[1],m.vertices[tri->v[2]].p[2]);

	//normal directions for each point in the triangle
	na = m.vertices[tri->v[0]].n;
	nb = m.vertices[tri->v[1]].n;
	nc = m.vertices[tri->v[2]].n;

	//(1) get normal direction at intersection
	v0 = b - a; 
	v1 = c - a; 
	v2 = current_P - a;
    d00 = v0 * v0;
    d01 = v0 * v1;
    d11 = v1 * v1;
    d20 = v2 * v0;
    d21 = v2 * v1;
    denom = d00 * d11 - d01 * d01;

    //assign Barycentric coords
    v = (d11 * d20 - d01 * d21) / denom;
    w = (d00 * d21 - d01 * d20) / denom;
    u = 1.0f - v - w;

    //get normal direction of point using Barycentric coords
    nx = u * na + v * nb + w * nc;

	//(2) get the color of model, lights will modify this value
	color = m.mat_color;

    //dot product of light and normal 
    //c
    // \  |
    //  \ | 
    //	 \|
    // -----------
	//(3) get light position/direction/relflected
    light_position = Point3d(	light0_position[0], 
     							light0_position[1], 
     							light0_position[2]	);
	light_direction = (light_position - current_P).normalize();
	light_reflected = 2 * (nx * light_direction) * nx - light_direction;

	//set up ambient, diffuse, and specular affects on color
	//using their calculated coefficients, color, and material
	dc = clamp(light_direction * nx);
	diffuse = dc * color;

	ac = 0.255; 
	ambient = color * ac;	

	sc = pow(std::max((light_reflected * light_direction), 0.0), m.mat_shininess);
	specular = sc * m.mat_specular;

	//set up color of point
	color = ambient + diffuse + specular;

	if (inshadow(nx))		//if in shadow, reduce the color
	{
		color = color / 3; //ugh >_<; magic number
	}

	return color;
}

//check if a point p is in shadow
//return 0 is p is not in the shadow
//a value between 0 and 1 to indicate the "probability"
//that p is in full darkness
double RayTracer::inshadow(const Point3d& p)
{
	//
	// PA7 TODO: implement this
	//           note: this is different from your PA6
	//
	// hint: sample N rays to the light (using m_lights)
	//       let n of N rays be visible betwen the light and p
	//       return n/N
	//
	//       in this file N is SHADOW_SAMPLE_SIZE (defined above)
	//

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
	{
		for (int j = 0; j < w; j++)
			fprintf(f, "%d %d %d ", toInt(m_img[i][j][0]), toInt(m_img[i][j][1]), toInt(m_img[i][j][2]));

		show_progress_bar(i*1.0 / h);
	}
	show_progress_bar(1.0f);
	fclose(f);
}


