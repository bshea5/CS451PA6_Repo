#include "raytracer.h"
#include <algorithm>
#include <cfloat>

//some helper functions
inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; }

inline int toInt(double x){ return int( clamp(x)*255 + .5); }

//global variable vector, :( , so dirty
Point3d current_P;

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

inline void show_progress_bar(double ratio)
{
	// Show the load bar.
	static const int w = 50;
	int   c = ratio * w;

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
				if (X.first != NULL && X.second != NULL)
				{
					Vector3d rc = raycolor(*X.first, X.second, r);
					color = rc + color;
				}
			}

			m_img[y][x] = color / nray_per_pixel;
		}//end of x

		show_progress_bar(y*1.0 / img_h);

	}//end of y

	cout << "\nRendering is done." << endl;
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
pair<model*, triangle*> RayTracer::intersect(Ray r)
{
	double min_dist = FLT_MAX;
	triangle* closest_T = NULL;
	model* closest_M = NULL;

	for (list<model>::iterator i = m_models.begin(); i != m_models.end(); i++)
	{
		triangle * t = closest_intersect(*i, r);
		if (t != NULL)
		{
			Point3d x;
			intersect(*i, t, r, x);
			double dist = (x - r.o).normsqr();
			if (dist < min_dist)
			{
				min_dist = dist;
				closest_T = t;
				closest_M = &*i;
				current_P = x;
			}
		}
	}
	return make_pair(closest_M, closest_T);
}

//returns a triangle in model m that intersect ray r
//return NULL if no intersection is found
triangle* RayTracer::intersect(model& m, Ray r)
{
	for (int i = 0; i < m.t_size; i++)
	{
		Point3d x;
		if ( intersect(m, &m.tris[i], r, x) )
		{
			return &m.tris[i];
		}
	}
	return NULL;
}

//returns a triangle in model m that make closest intersection with ray r
//return NULL if no intersection is found
triangle* RayTracer::closest_intersect(model& m, Ray r)
{
	double min_dist = FLT_MAX;
	triangle* closest_T = NULL;
	for (int i = 0; i < m.t_size; i++)
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
// determine the color of ray r, by analizing the intersection between t and r 
// 
Vector3d RayTracer::raycolor(model& m, triangle* tri, const Ray& r)
{
	Vector3d color; 			//color of point to return
	Vector3d a, b, c;			//triangle points
	Vector3d v0, v1, v2;
	Vector3d nx, na, nb, nc;	//normal directions
	Vector3d al, dl, sl; 		//ambient, diffuse, and specular lights
	double ac, dc, sc;			//coefficients
	double shade;
	Vector3d diffuse, ambient, specular;
	Vector3d light_direction, light_reflected;
	Point3d light_position;
	double u, v, w; 			//Barycentric coords

	//triangle points in vector form
	a = Vector3d(m.vertices[tri->v[0]].p[0],m.vertices[tri->v[0]].p[1],m.vertices[tri->v[0]].p[2]);
	b = Vector3d(m.vertices[tri->v[1]].p[0],m.vertices[tri->v[1]].p[1],m.vertices[tri->v[1]].p[2]);
	c = Vector3d(m.vertices[tri->v[2]].p[0],m.vertices[tri->v[2]].p[1],m.vertices[tri->v[2]].p[2]);

	//normal directions for each point in the triangle
	na = m.vertices[tri->v[0]].n;
	nb = m.vertices[tri->v[1]].n;
	nc = m.vertices[tri->v[2]].n;

	//hint: (1) what is the normal direction at the intersection between t and r
	// 				nx = w1n1 + w2n2 + w3n3, where the w's are barycentric coords
	//      (2) what is the color of this triangle? 
	//      (3) where is the light?
	//      (4) is the intersection in shadow? (call inshadow(const Point3d& p))

	//(1) get normal direction at intersection
	v0 = b - a; 
	v1 = c - a; 
	v2 = current_P - a;
    double d00 = v0 * v0;
    double d01 = v0 * v1;
    double d11 = v1 * v1;
    double d20 = v2 * v0;
    double d21 = v2 * v1;
    double denom = d00 * d11 - d01 * d01;

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

	//get ambient, diffuse, specular lights from gliLight.h
	// al = Vector3d(	light0_ambient[0],	light0_ambient[1],	light0_ambient[2]	);
	// dl = Vector3d(	light0_diffuse[0],	light0_diffuse[1],	light0_diffuse[2]	);
	// sl = Vector3d(	light0_specular[0], light0_specular[1], light0_specular[2]	);

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
		//return color;
	}

	return color;
}

//check if a point p is in shadow
bool RayTracer::inshadow(const Point3d& p)
{
	//TODO: implement this
	//shoot a ray from light to the point, if it intersects with a model, inShadow is true
	Ray shadow_ray;					//ray to shoot back at light
	pair<model*, triangle*> X;		//holds model the ray intersects with, <NULL, NULL> if nothing
	Point3d light_position;
	double dist_to_light, dist_to_x;

	//get light position as a point3d
    light_position = Point3d(	light0_position[0], 
 								light0_position[1], 
 								light0_position[2]	);

	//set up shadow ray
	shadow_ray.v = (light_position - current_P).normalize();
	shadow_ray.o = current_P + shadow_ray.v * 0.001;

	dist_to_light = (light_position - current_P).norm();

	//std::cout << "before: " << current_P << "\nafter: " << current_P << std::endl;
	//loop through all models, for each model, intersect(model, ray)
	X = intersect(shadow_ray);
	if (X.first != NULL && X.second != NULL)
	{
		//std::cout << "stuff" << std::endl;
		dist_to_x = (current_P - shadow_ray.o).norm();
		if (dist_to_x < dist_to_light)
		{
			//std::cout << "shadowed" << std::endl;
			return true;
		}
	}
	return false;

}

//save rendered image to file
bool RayTracer::save2file(const string& filename)
{
	FILE* f = fopen(filename.c_str(), "w");         // Write image to PPM file.

	int h = m_img.size();
	if (h == 0) return true; //nothing to save...

	int w = m_img.front().size();
	
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);

	for (int i = 0; i < h; i++) 
		for (int j = 0; j < w; j++) 
			fprintf(f, "%d %d %d ", toInt(m_img[i][j][0]), toInt(m_img[i][j][1]), toInt(m_img[i][j][2]));

	fclose(f);
}


