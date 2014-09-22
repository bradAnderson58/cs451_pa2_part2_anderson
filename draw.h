//------------------------------------------------------------------------------
//  Copyright 2007-2014 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#pragma once

#include <GL/gli.h>
#include <time.h>

#include <map>
using namespace std;

#include "Basic.h"
#include "model.h"
using namespace mathtool;

//-----------------------------------------------------------------------------
//variables used in rendering

bool showWire=false; //on/off wireframe
bool showSharpEdgeOnly=false;
bool randomColor=false;
bool background=true; //on/off background
bool light=true; //on/off background

//Store display IDs and model colors
map<model*,int> model_solid_gids;
map<model*,int> model_wire_gids;
map<model*,Vector3d> model_colors;

//For use in 3d picking possibly
GLdouble* modelM = new GLdouble[16];	//store the modelview matrix in this
GLdouble* proj = new GLdouble[16];		//store the projection matrix in this
GLint* view = new GLint[4];				//store the viewport matrix in this
GLdouble* winx = new GLdouble;			//store the values of the screen coordinates in these
GLdouble* winy = new GLdouble;
GLdouble* winz = new GLdouble;

GLdouble* objx = new GLdouble;			//store the values of the 3d points in these
GLdouble* objy = new GLdouble;
GLdouble* objz = new GLdouble;

Point3d* selected;						//a pointer to a selected node in FFD_lattice
bool sel = false;						//is a node selected or no?

GLdouble projectedZ;					//z coordinate for unproject
Point3d oldPoint;						//old point for optimizing deformation code
int selIndex;							//index of point for optimizing deformation code

//-----------------------------------------------------------------------------
//this defines the size of the lattice
extern unsigned int lattice_nx, lattice_ny, lattice_nz;
extern vector<Point3d> FFD_lattice; //This stores all lattice nodes, FFD_lattice has size = (lattice_nx X lattice_ny X lattice_nz)
extern vector<vector<double>> FFD_parameterization; //This stores all parameterized coordinates of all vertices from the model
                                             //FFD_parameterization has size = models.front().v_size

//-----------------------------------------------------------------------------

inline void DisplayLattice()
{
	//glDisable(GL_LIGHTING); //disable lighting

	//TODO: draw lattice nodes using FFD_lattice *done
	
	//For each point in FFD_lattice, draw a point in world
	glColor3f(1, 0, 0);
	glPointSize(6.0);
	glBegin(GL_POINTS);
	for (Point3d spot : FFD_lattice){
		glVertex3f(spot[0], spot[1], spot[2]);       
	}
	glEnd();

	//If a point is selected, draw the selected node larger and different color
	if (sel){  
		Point3d t = *selected;
		glColor3f(0, 0, 1);
		glPointSize(10.0);
		glBegin(GL_POINTS);
		glVertex3f(t[0], t[1], t[2]);
		glEnd();
	}

	//TODO: draw lattice edges using FFD_lattice *done
	glColor3f(0, 1, 0);
	int ctr = 0;
	Point3d cur;
	//I did this because math is hard----------------
	vector<vector<vector<Point3d>>> sanity;
	for (int x = 0; x < lattice_nx; ++x){
		vector<vector<Point3d>> temp1;
		for (int y = 0; y < lattice_ny; ++y){
			vector<Point3d> temp2;
			for (int z = 0; z < lattice_nz; ++z){
				temp2.push_back(FFD_lattice[ctr]);
				++ctr;
			}//z
			temp1.push_back(temp2);
		}//y
		sanity.push_back(temp1);
	}//x
	//------------------------------------------------
	//For each point, draw lines to it's neighbors
	for (int x = 0; x < lattice_nx; x++){
		for (int y = 0; y < lattice_ny; y++){
			for (int z = 0; z < lattice_nz; z++){

				//If not the end of the x row, draw a line to the next point on your x row
				if (x != lattice_nx - 1){
					glBegin(GL_LINES);
					cur = sanity[x][y][z];
					glVertex3f(cur[0], cur[1], cur[2]);
					cur = sanity[x + 1][y][z];
					glVertex3f(cur[0], cur[1], cur[2]);
					glEnd();
				}

				//If not the end of the y row, draw line to next point along y axis
				if (y != lattice_ny - 1){
					glBegin(GL_LINES);
					cur = sanity[x][y][z];
					glVertex3f(cur[0], cur[1], cur[2]);
					cur = sanity[x][y + 1][z];
					glVertex3f(cur[0], cur[1], cur[2]);
					glEnd();
				}

				//we know the drill
				if (z != lattice_nz - 1){
					glBegin(GL_LINES);
					cur = sanity[x][y][z];
					glVertex3f(cur[0], cur[1], cur[2]);
					cur = sanity[x][y][z + 1];
					glVertex3f(cur[0], cur[1], cur[2]);
					glEnd();
				}
			}//z
		}//y
	}//x
		
}

inline void DisplayModel(model& M, bool randcolor=false)
{
	//draw
	if (randcolor){
		if (model_colors.find(&M) == model_colors.end())
			model_colors[&M] = Vector3d(drand48() + 0.5, drand48() + 0.5, drand48(), +0.5).normalize() + Vector3d(0.25, 0.25, 0.25);
		glColor3dv(model_colors[&M].get());
	}
	
	//Draw facets
	glEnable( GL_POLYGON_OFFSET_FILL );
	glPolygonOffset( 0.5f, 0.5f );
	//for(list<polygon>::iterator i=M.polys.begin();i!=M.polys.end();i++)
	glBegin(GL_TRIANGLES);
	for(unsigned int i=0;i<M.t_size;i++)
	{
        const triangle & t=M.tris[i];
        glNormal3dv(M.tris[i].n.get());
        for(int k=0;k<3;k++)
		{
            const Point3d& pt=M.vertices[t.v[k]].p;

			//send pt to OpenGL
            glVertex3d(pt[0],pt[1],pt[2]);
        }
	}
	glEnd();
	glDisable( GL_POLYGON_OFFSET_FILL );
}

inline void DisplayModelWireFrame(model& M, bool randcolor=false)
{
    //Draw Edges
    if(showWire)
	{
		glBegin(GL_LINES);
        for(uint i=0;i<M.e_size;i++){
            glColor3f(0,0,0);
            const edge & e=M.edges[i];
            if(e.fid.size()==2){//normal case, check if e is sharp
                triangle& f1=M.tris[e.fid.front()];
                triangle& f2=M.tris[e.fid.back()];
                if(fabs(1-f1.n*f2.n)<1e-2){
                    if(showSharpEdgeOnly) continue; //not sharp
                    else
                        glColor3f(0.7f,0.7f,0.7f);
                }
            }

            Point3d& p1=M.vertices[e.vid[0]].p;
            Point3d& p2=M.vertices[e.vid[1]].p;
            glVertex3d(p1[0],p1[1],p1[2]);
            glVertex3d(p2[0],p2[1],p2[2]);
        }
        glEnd();
    }
}


//copied from meshlab
void DisplayBackground(void)
{
	float topcolor[]={1,1,1};
	float bottomcolor[]={1,1,0.5};
	
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(-1,1,-1,1,-1,1);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glPushAttrib(GL_ENABLE_BIT);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glBegin(GL_TRIANGLE_STRIP);
		glColor3fv(topcolor);  	glVertex2f(-1, 1);
		glColor3fv(bottomcolor);	glVertex2f(-1,-1);
		glColor3fv(topcolor);	glVertex2f( 1, 1);
		glColor3fv(bottomcolor);	glVertex2f( 1,-1);
	glEnd();
	
	glPopAttrib();
	glPopMatrix(); // restore modelview
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
}

void drawAll()
{
    glEnable(GL_LIGHTING);

    //show the inputs
    glColor3f(1,1,1);
    if(light) glEnable(GL_LIGHTING);
    else glDisable(GL_LIGHTING);

    for(list<model>::iterator i=models.begin();i!=models.end();i++){
        DisplayModel(*i,randomColor);
    }
    for(list<model>::iterator i=models.begin();i!=models.end();i++){
        DisplayModelWireFrame(*i);
    }

	//draw lattice
	DisplayLattice();
}

//-----------------------------------------------------------------------------
void Display( void )
{
    //Init Draw
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    if(background) DisplayBackground();
    
    
    glPushMatrix();
    glLoadIdentity();
    static GLfloat light_position1[] = {  100, 100, 100.0f, 1.0f };
    glLightfv(GL_LIGHT0, GL_POSITION, light_position1);
    static GLfloat light_position2[] = { -100, -100, 50.0f, 1.0f };
    glLightfv(GL_LIGHT1, GL_POSITION, light_position2);
    glPopMatrix();

    drawAll();

    glDisable(GL_LIGHTING);


}


//-----------------------------------------------------------------------------
// regular openGL callback functions
bool InitGL()
{
    // transparent
    glShadeModel(GL_SMOOTH);
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

	glEnable( GL_LINE_SMOOTH );
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );

    // others
    glEnable( GL_DEPTH_TEST);
	glEnable(GL_NORMALIZE);
    glClearColor(1.0, 1.0, 1.0, 0.0);

    //Let's have light!
    GLfloat Diffuse[] =  { 0.9f, 0.9f, 0.9f, 1.0f };
    glMaterialfv(GL_FRONT, GL_DIFFUSE, Diffuse);
    glColorMaterial(GL_FRONT, GL_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);

    GLfloat WhiteLight[] =  { 0.75f, 0.75f, 0.75f, 1.0f };
    glLightfv(GL_LIGHT0,GL_DIFFUSE,WhiteLight);
    glLightfv(GL_LIGHT1,GL_DIFFUSE,WhiteLight);

    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);

    return true;
}

void Reshape( int w, int h)
{

    glViewport( 0, 0, (GLsizei)w, (GLsizei)h );

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();

	//Perspective view
    gluPerspective( 60, 1.0f*w/h, R/100, R*100 );

	//Othogonal view
	//glOrtho(-R * 1.5, R * 1.5, -R * 1.5, R * 1.5, -R * 100, R * 100);


    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

}

void Mouse(int button, int state, int x, int y)
{
	// This handles mouse click 
	// 
	// TODO: check if the user clicks on one of the lattice nodes *done
	//

	sel = false;
	if (button == 0 && state == 0){		//Only do this if the left mouse button is down pressed
		int ctr = 0;
		int best = -1;
		float score = 999;
		float total = 0;
		//for each point in the lattice, translate world coordinates to screen coordinates using gluProject()
		//If it seems to be being clicked on, assign selected pointer to this point
		for (Point3d point : FFD_lattice){

			GLdouble obx = point[0];
			GLdouble oby = point[1];
			GLdouble obz = point[2];

			glGetDoublev(GL_MODELVIEW_MATRIX, modelM);
			glGetDoublev(GL_PROJECTION_MATRIX, proj);
			glGetIntegerv(GL_VIEWPORT, view);

			gluProject(obx, oby, obz, modelM, proj, view, winx, winy, winz);

			if (abs(x - (float)*winx) < 5 && abs((view[3] - y - 10) - (float)*winy) < 5){
				total = (abs(x - (float)*winx)) + abs((view[3] - y - 10) - (float)*winy);

				if (total < score)
				{
					best = ctr;
					projectedZ = *winz;
					
				}
			}
			ctr++;
		}
		if (best != -1){
			selected = &FFD_lattice[best];
			sel = true;
			oldPoint = FFD_lattice[best];
			selIndex = best;
		}
	}
}

bool Motion(int x, int y)
{
	//
	// This handles mouse motion when a button is pressed 
	// 
	// TODO: check if the user is moving a selected lattice node  *done
	//       if so move the node to a new location
	//

	//Only do this if a point is selected (useris pressing down left mouse button on valid point)
	//Use gluUnProject to translate to selected node to current mouse coordinates
	if (sel){

		GLdouble wnx = x;
		GLdouble wny = view[3] - y;
		
		gluUnProject(wnx, wny, projectedZ, modelM, proj, view, objx, objy, objz);

		selected[0][0] = *objx;
		selected[0][1] = *objy;
		selected[0][2] = *objz;

		// TODO: recompute the position of every vertex in the model 
		//       i.e., using FFD_parameterization and FFD_lattice 
		//       Note, you only need to do this to the first model in "models"
		//

		model& m = models.front();
		double xthing, ything, zthing;
		for (int i = 0; i < m.v_size; ++i){  //DO this for every point in the model
			xthing = ything = zthing = 0;
			int inverseInd = FFD_lattice.size() - selIndex - 1;
			xthing = m.vertices[i].p[0] - (oldPoint[0] * FFD_parameterization[i][inverseInd]) + (selected[0][0] * FFD_parameterization[i][inverseInd]);
			ything = m.vertices[i].p[1] - (oldPoint[1] * FFD_parameterization[i][inverseInd]) + (selected[0][1] * FFD_parameterization[i][inverseInd]);
			zthing = m.vertices[i].p[2] - (oldPoint[2] * FFD_parameterization[i][inverseInd]) + (selected[0][2] * FFD_parameterization[i][inverseInd]);
			/*vector<double> weights = FFD_parameterization[i];   //this is a vector of lattice weights for this model vertex
			int ctr = FFD_lattice.size()-1;
			for (Point3d point : FFD_lattice){					//for each point in the lattice, add up the product of the point times the weight
				xthing += (point[0] * weights[ctr]);
				ything += (point[1] * weights[ctr]);
				zthing += (point[2] * weights[ctr]);
				ctr--;
			}*/
			m.vertices[i].p[0] = xthing;
			m.vertices[i].p[1] = ything;
			m.vertices[i].p[2] = zthing;
		}
		oldPoint = FFD_lattice[selIndex];
		glutPostRedisplay();
		return true;   //if a point is being moved, don't rotate the model!!!!!!!!!
	}
	return false;
}

void PassiveMotion(int x, int y)
{
	// This handles mouse motion when mouse button is NOT pressed
	// does nothing now...

	glutPostRedisplay();
}

//Used for simulation/anitmation. 
void TimerCallback(int value)
{
    //in simuation state
    glutPostRedisplay();
    glutTimerFunc(30, TimerCallback,value);
}

//Handle keyboard events
void Keyboard( unsigned char key, int x, int y )
{
    // find closest colorPt3D if ctrl is pressed...
    switch( key ){
        case 27: exit(0);
        case 'w' : showWire=!showWire; break;
        case 'r' : randomColor=!randomColor; break;
		case 'R' : model_colors.clear(); break;
		case 'L' : light=!light; break;
		case 'b' : background=!background; break;
		case 'S' : showSharpEdgeOnly=!showSharpEdgeOnly;
		           for(map<model*,int>::iterator i=model_wire_gids.begin();i!=model_wire_gids.end();i++) glDeleteLists(i->second,1);
		           model_wire_gids.clear();
		           break;
    }
    glutPostRedisplay();
}



