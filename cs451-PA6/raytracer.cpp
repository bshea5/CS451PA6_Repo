#include "raytracer.h"
#include <algorithm>
#include <cfloat>

//some helper functions
inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; }

inline int toInt(double x){ return int( clamp(x)*255 + .5); }

//global variable vector, i hate myself
// std::vector<Vector3d> lines;

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
	Vector3d near, far;

	GLdouble pixelWidth = viewport[2] / m_img.front().size();
	GLdouble pixelHeight = viewport[3] / m_img.size(); 
	GLdouble winX = x + (GLdouble)(rand() % 10) / 10.0; //get random decimal value [0,1] for offset
	GLdouble winY = y + (GLdouble)(rand() % 10) / 10.0;
	winX *= pixelWidth;
	winY *= pixelHeight;

	//get near
	gluUnProject(	winX,
					winY,
				 	0,				//get the near plane
				 	modelMatrix,	//model
				 	projMatrix, 	//projection
				 	viewport,		//viewport
				 	&near[0],
				 	&near[1],
				 	&near[2]	);

	//get far
	gluUnProject(	winX,
					winY,
				 	1,				//get the far plane
				 	modelMatrix,	
				 	projMatrix, 	
				 	viewport,		
				 	&far[0],
				 	&far[1],
				 	&far[2]		);

	//ray vector
	//far - near
	r.v = (far - near).normalize();

	//ray origin
	//add camera position to origin
	r.o = Point3d(	camera_pos[0] + near[0], 
					camera_pos[1] + near[1], 
					camera_pos[2] + near[2]	);

	all_rays.push_back(r);
	return r;
}

//returns a model and triangle that intersect ray r
//return pair<NULL,NULL> if no intersection is found
pair<model*, triangle*> RayTracer::intersect(Ray r)
{
	double min_dist = FLT_MAX;
	triangle * closest_T = NULL;
	model * closest_M = NULL;

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
	triangle * closest_T=NULL;
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
	//TODO:Implement me!
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
		intersection_points.push_back(x);
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
	Vector3d light_position;	//
	double u, v, w; 			//Barycentric coords
	double scalar;

	//triangle points in vector form
	a = Vector3d(m.vertices[tri->v[0]].p[0],m.vertices[tri->v[0]].p[1],m.vertices[tri->v[0]].p[2]);
	b = Vector3d(m.vertices[tri->v[1]].p[0],m.vertices[tri->v[1]].p[1],m.vertices[tri->v[1]].p[2]);
	c = Vector3d(m.vertices[tri->v[2]].p[0],m.vertices[tri->v[2]].p[1],m.vertices[tri->v[2]].p[2]);

	//normal directions for each point in the triangle
	na = m.vertices[tri->v[0]].n;
	nb = m.vertices[tri->v[1]].n;
	nc = m.vertices[tri->v[2]].n;

	//TODO: implement this

	//
	//hint: (1) what is the normal direction at the intersection between t and r
	// nx = w1n1 + w2n2 + w3n3, where the w's are barycentric coords
	//      (2) what is the color of this triangle? 
	//      (3) where is the light?
	//      (4) is the intersection in shadow? (call inshadow(const Point3d& p))
	//(1)
	v0 = b - a; 
	v1 = c - a; 
	v2 = r.o - a;
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
    nx = v * na + w * nb + u * nc;

    //dot product of light and normal 
    // \  |
    //  \ |
    //	 \|
    // -----------
    // light_position = Vector3d(	gliLight.light0_position[0], 
    // 							gliLight.light0_position[1], 
    // 							gliLight.light0_position[2]	);
    //scalar = nx * 

	return color;
	//return m.mat_color;
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


