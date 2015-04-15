#include <stdlib.h>
#include <stdio.h>
#include <GLUT/glut.h>
#include "solver.c"
/* macros */

#define IX(i,j,k) ((i)+(N+2)*(j) + (N+2)*(N+2) *(k))

/* external definitions (from solver.c) */

extern void dens_step ( int N, float * x, float * x0, float * u, float * v, float * w, float diff, float dt );
extern void temp_step ( int N, float * u, float * v, float * w, float diff, float dt );
extern void vel_step ( int N, float * u, float * v, float * w, float * u0, float * v0, float* w0, float visc, float dt );


/* global variables */

static int N;
static float dt, diff, visc;
static float force, source, temp_source;
static int dvel, dAxis, addsource;

static float * u, * v, * u_prev, * v_prev, * w, * w_prev;
static float * dens, * dens_prev;
static float * temp, * temp_prev;
static float * buo;

static int win_id;
static int win_x, win_y;
static int mouse_down[3];
static int omx, omy, mx, my;
static float source_alpha =  0.05; //for displaying density

GLfloat trans[3];
GLfloat rot[2];

/*
 ----------------------------------------------------------------------
 free/clear/allocate simulation data
 ----------------------------------------------------------------------
 */

float clamp(float x) {
    return x > 360.0f ? x-360.0f : x < -360.0f ? x+=360.0f : x;
}

static void free_data ( void )
{
    if ( u ) free ( u );
    if ( v ) free ( v );
    if ( w ) free ( w );
    if ( u_prev ) free ( u_prev );
    if ( v_prev ) free ( v_prev );
    if ( w_prev ) free ( w_prev );
    if ( dens ) free ( dens );
    if ( dens_prev ) free ( dens_prev );
    if ( temp ) free (temp);
    if ( temp_prev ) free ( temp_prev );
    if ( buo ) free ( buo );
}

static void clear_data ( void )
{
    int i, size=(N+2)*(N+2)*(N+2);
    
    for ( i=0 ; i<size ; i++ ) {
        u[i] = v[i] = w[i] = u_prev[i] = v_prev[i] = w_prev[i] = dens[i] = dens_prev[i] = buo[i] = 0.0f;
        temp[i] = temp_prev[i] = 0;
    }
}

static int allocate_data ( void )
{
    int size = (N+2)*(N+2)*(N+2);
    
    u			= (float *) malloc ( size*sizeof(float) );
    v			= (float *) malloc ( size*sizeof(float) );
    w			= (float *) malloc ( size*sizeof(float) );
    u_prev		= (float *) malloc ( size*sizeof(float) );
    v_prev		= (float *) malloc ( size*sizeof(float) );
    w_prev		= (float *) malloc ( size*sizeof(float) );
    dens		= (float *) malloc ( size*sizeof(float) );
    dens_prev	= (float *) malloc ( size*sizeof(float) );
    temp        = (float *) malloc ( size*sizeof(float) );
    temp_prev   = (float *) malloc ( size*sizeof(float) );
    buo         = (float *) malloc ( size*sizeof(float) );
    
    if ( !u || !v || !w || !u_prev || !v_prev || !w_prev || !dens || !dens_prev || !temp || !temp_prev || !buo ) {
        fprintf ( stderr, "cannot allocate data\n" );
        return ( 0 );
    }
    
    return ( 1 );
}


/*
 ----------------------------------------------------------------------
 OpenGL specific drawing routines
 ----------------------------------------------------------------------
 */

static float* getColour(

static void draw_axis( void ) {
    glLineWidth ( 1.0f );
    glBegin (GL_LINES);
    
    glColor3f(1.0f, 1.0f, 1.0f);
    glVertex3f (0.f, 0.f, 0.f);
    glVertex3f (1.0f, 0.f, 0.f);
    
    glVertex3f (0.f, 0.f, 0.f);
    glVertex3f (0.f, 1.0f, 0.f);

    glVertex3f (0.f, 0.f, 0.f);
    glVertex3f (0.f, 0.f, 1.0f);
    
    glVertex3f (1.f, 1.f, 1.f);
    glVertex3f (0.f, 1.f, 1.0f);
    
    glVertex3f (1.f, 1.f, 1.f);
    glVertex3f (1.f, 1.f, 0.0f);
    
    glVertex3f (1.f, 1.f, 1.f);
    glVertex3f (1.f, 0.f, 1.0f);
    
    glVertex3f (1.f, 0.f, 1.f);
    glVertex3f (0.f, 0.f, 1.0f);
    
    glVertex3f (1.f, 0.f, 1.f);
    glVertex3f (1.f, 0.f, 0.0f);
    
    glVertex3f (0.f, 1.f, 1.f);
    glVertex3f (0.f, 0.f, 1.f);

    glVertex3f (0.f, 1.f, 1.f);
    glVertex3f (0.f, 1.f, 0.f);

    glVertex3f (1.f, 1.f, 0.f);
    glVertex3f (1.f, 0.f, 0.f);
    
    glVertex3f (1.f, 1.f, 0.f);
    glVertex3f (0.f, 1.f, 0.f);
    glEnd();
}

static void draw_velocity ( void )
{
    int i, j, k;
    float x, y, z, h;
    
    h = 1.0f/N;
    
    glColor3f ( 1.0f, 1.0f, 1.0f );
    glLineWidth ( 1.0f );
    
    glBegin ( GL_LINES );
    
    for ( i=1 ; i<=N ; i++ ) {
        x = (i-0.5f)*h;
        for ( j=1 ; j<=N ; j++ ) {
            y = (j-0.5f)*h;
            for (k=1; k <= N ; k++) {
                z =(k-0.5)*h;
                glVertex3f ( x, y, z );
                glVertex3f ( x+u[IX(i,j,k)], y+v[IX(i,j,k)], z+v[IX(i,j,k)] );
            }
        }
    }
    
    glEnd ();
}

static void draw_density ( void )
{
    int i, j, k;
    float x, y, z, h, d000, d010, d100, d110, d001, d011, d101, d111;;
    
    h = 1.0f/N;
    
    glBegin ( GL_QUADS );
    
    for ( i=0; i<=N; i++ ) {
        x = (i-0.5f)*h;
        for ( j=0; j<=N; j++ ) {
            y = (j-0.5f)*h;
            for ( k=0; k<=N; k++ ) {
                z = (k-0.5f)*h;
                
                d000 = dens[IX(i,j,k)];
                d010 = dens[IX(i,j+1,k)];
                d100 = dens[IX(i+1,j,k)];
                d110 = dens[IX(i+1,j+1,k)];
                
                d001 = dens[IX(i,j,k+1)];
                d011 = dens[IX(i,j+1,k+1)];
                d101 = dens[IX(i+1,j,k+1)];
                d111 = dens[IX(i+1,j+1,k+1)];
                
                // draw density as a cube of quads
                
                glColor4f ( d111, d111, d111, source_alpha ); glVertex3f ( x+h,y+h,z+h );
                glColor4f ( d011, d011, d011, source_alpha ); glVertex3f ( x, y+h, z+h);
                glColor4f ( d001, d001, d001, source_alpha ); glVertex3f ( x, y, z+h );
                glColor4f ( d101, d101, d101, source_alpha ); glVertex3f ( x+h, y, z+h );
                
                glColor4f ( d110, d110, d110, source_alpha ); glVertex3f ( x+h, y+h, z );
                glColor4f ( d111, d111, d111, source_alpha ); glVertex3f ( x+h,y+h,z+h );
                glColor4f ( d101, d101, d101, source_alpha ); glVertex3f ( x+h, y, z+h );
                glColor4f ( d100, d100, d100, source_alpha ); glVertex3f ( x+h, y, z );
                
                glColor4f ( d010, d010, d010, source_alpha ); glVertex3f ( x, y+h, z );
                glColor4f ( d110, d110, d110, source_alpha ); glVertex3f ( x+h, y+h, z );
                glColor4f ( d100, d100, d100, source_alpha ); glVertex3f ( x+h, y, z );
                glColor4f ( d000, d000, d000, source_alpha ); glVertex3f ( x, y, z );
                
                glColor4f ( d011, d011, d011, source_alpha ); glVertex3f ( x, y+h, z+h);
                glColor4f ( d010, d010, d010, source_alpha ); glVertex3f ( x, y+h, z );
                glColor4f ( d000, d000, d000, source_alpha ); glVertex3f ( x, y, z );
                glColor4f ( d001, d001, d001, source_alpha ); glVertex3f ( x, y, z+h );
                
                glColor4f ( d100, d100, d100, source_alpha ); glVertex3f ( x+h, y, z );
                glColor4f ( d000, d000, d000, source_alpha ); glVertex3f ( x, y, z );
                glColor4f ( d001, d001, d001, source_alpha ); glVertex3f ( x, y, z+h );
                glColor4f ( d101, d101, d101, source_alpha ); glVertex3f ( x+h, y, z+h );
                
                glColor4f ( d110, d110, d110, source_alpha ); glVertex3f ( x+h, y+h, z );
                glColor4f ( d010, d010, d010, source_alpha ); glVertex3f ( x, y+h, z );
                glColor4f ( d011, d011, d011, source_alpha ); glVertex3f ( x, y+h, z+h);
                glColor4f ( d111, d111, d111, source_alpha ); glVertex3f ( x+h, y+h, z+h );
            }
        }
    }
    glEnd ();
}

static void draw_temperature (void)
{
    int i, j, k;
    float x, y, z, h, d000, d010, d100, d110, d001, d011, d101, d111;;
    
    h = 1.0f/N;
    
    glBegin ( GL_QUADS );
    
    for ( i=0; i<=N; i++ ) {
        x = (i-0.5f)*h;
        for ( j=0; j<=N; j++ ) {
            y = (j-0.5f)*h;
            for ( k=0; k<=N; k++ ) {
                z = (k-0.5f)*h;
                
                d000 = temp[IX(i,j,k)];
                d010 = temp[IX(i,j+1,k)];
                d100 = temp[IX(i+1,j,k)];
                d110 = temp[IX(i+1,j+1,k)];
                
                d001 = temp[IX(i,j,k+1)];
                d011 = temp[IX(i,j+1,k+1)];
                d101 = temp[IX(i+1,j,k+1)];
                d111 = temp[IX(i+1,j+1,k+1)];
                
                
                // draw density as a cube of quads (6 faces)
                
                glColor4f ( d111, 0, 22-d111, source_alpha ); glVertex3f ( x+h,y+h,z+h );
                glColor4f ( d011, 0, 22-d011, source_alpha ); glVertex3f ( x, y+h, z+h);
                glColor4f ( d001, 0, 22-d001, source_alpha ); glVertex3f ( x, y, z+h );
                glColor4f ( d101, 0, 22-d101, source_alpha ); glVertex3f ( x+h, y, z+h );
                
                glColor4f ( d110, 0, 22-d110, source_alpha ); glVertex3f ( x+h, y+h, z );
                glColor4f ( d111, 0, 22-d111, source_alpha ); glVertex3f ( x+h,y+h,z+h );
                glColor4f ( d101, 0, 22-d101, source_alpha ); glVertex3f ( x+h, y, z+h );
                glColor4f ( d100, 0, 22-d100, source_alpha ); glVertex3f ( x+h, y, z );
                
                glColor4f ( d010, 0, 22-d010, source_alpha ); glVertex3f ( x, y+h, z );
                glColor4f ( d110, 0, 22-d110, source_alpha ); glVertex3f ( x+h, y+h, z );
                glColor4f ( d100, 0, 22-d100, source_alpha ); glVertex3f ( x+h, y, z );
                glColor4f ( d000, 0, 22-d000, source_alpha ); glVertex3f ( x, y, z );
                
                glColor4f ( d011, 0, 22-d011, source_alpha ); glVertex3f ( x, y+h, z+h);
                glColor4f ( d010, 0, 22-d010, source_alpha ); glVertex3f ( x, y+h, z );
                glColor4f ( d000, 0, 22-d000, source_alpha ); glVertex3f ( x, y, z );
                glColor4f ( d001, 0, 22-d001, source_alpha ); glVertex3f ( x, y, z+h );
            
                glColor4f ( d100, 0, 22-d100, source_alpha ); glVertex3f ( x+h, y, z );
                glColor4f ( d000, 0, 22-d000, source_alpha ); glVertex3f ( x, y, z );
                glColor4f ( d001, 0, 22-d001, source_alpha ); glVertex3f ( x, y, z+h );
                glColor4f ( d101, 0, 22-d101, source_alpha ); glVertex3f ( x+h, y, z+h );
                
                glColor4f ( d110, 0, 22-d110, source_alpha ); glVertex3f ( x+h, y+h, z );
                glColor4f ( d010, 0, 22-d010, source_alpha ); glVertex3f ( x, y+h, z );
                glColor4f ( d011, 0, 22-d011, source_alpha ); glVertex3f ( x, y+h, z+h);
                glColor4f ( d111, 0, 22-d111, source_alpha ); glVertex3f ( x+h, y+h, z+h );
            }
        }
    }
    
    glEnd ();
}

/*
 ----------------------------------------------------------------------
 relates mouse movements to forces sources
 ----------------------------------------------------------------------
 */

static void get_from_UI ( float * d, float * u, float * v, float *t )
{
    int i, j, k, size = (N+2)*(N+2)*(N+2);
    
    for ( i=0 ; i<size ; i++ ) {
        u[i] = v[i] = d[i] = t[i] = 0.0f;
    }
    
    if ( !mouse_down[0] && !mouse_down[2] ) return;
    
    i = (int)((       mx /(float)win_x)*N+1);
    j = (int)(((win_y-my)/(float)win_y)*N+1);
    k = (int) (N/2);
    
    if ( i<1 || i>N || j<1 || j>N ) return;
    
    if ( mouse_down[0] ) {
        u[IX(i,j,k)] = force * (mx-omx);
        v[IX(i,j,k)] = force * (omy-my);
        w[IX(i,j,k)] = force * (omy-my);
    }
    
    if ( mouse_down[2] ) {
        d[IX(i,j,k)] = source;
        t[IX(i,j,k)] = temp_source;
    }
    
    omx = mx;
    omy = my;
    
    return;
}

/*
 ----------------------------------------------------------------------
 GLUT callback routines
 ----------------------------------------------------------------------
 */

static void key_func ( unsigned char key, int x, int y )
{
    switch ( key )
    {
        case 'c':
        case 'C':
            clear_data ();
            break;
            
        case 'q':
        case 'Q':
            free_data ();
            exit ( 0 );
            break;
            
        case 'v':
        case 'V':
            dvel = !dvel;
            break;
            
        case 't':
        case 'T':
            temp_source = -temp_source;
            break;
            
        case 'a':
        case 'A':
            dAxis = !dAxis;
            break;
            
        case 'x':       // 'X' key - add source at centre
            addsource = 1;
            break;
            
    }
}

static void mouse_func ( int button, int state, int x, int y )
{
    omx = mx = x;
    omx = my = y;
    
    mouse_down[button] = state == GLUT_DOWN;
}

static void motion_func ( int x, int y )
{
    mx = x;
    my = y;
}

static void reshape_func ( int width, int height )
{
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0, (float)width/height, 0.001, 100.0);
}

static void idle_func ( void )
{
        get_from_UI ( dens_prev, u_prev, v_prev, temp_prev );
    
    
    vel_step ( N, u, v, w, temp, u_prev, v_prev, w_prev, temp_prev, visc, dt, buo);
    dens_step ( N, dens, dens_prev, u, v, w, diff, dt );
    temp_step(N, temp, temp_prev, u, v, w, diff/100, dt);
    
    glutSetWindow ( win_id );
    glutPostRedisplay ();
}

static void pre_display ( void )
{
    glViewport ( 0, 0, win_x, win_y );
    glMatrixMode ( GL_PROJECTION );
    glLoadIdentity ();
    gluOrtho2D ( 0.0, 1.0, 0.0, 1.0 );
    glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
    glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
}

static void post_display ( void )
{
    glutSwapBuffers ();
}

static void display_func ( void )
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glPushMatrix();
    
    glTranslatef(0.0f, 0.0f, -2.5f);
    glRotatef(rot[0], 1.0f, 0.0f, 0.0f);
    glRotatef(rot[1], 0.0f, 1.0f, 0.0f);
    glTranslatef(-0.5, -0.5, -0.5);
    glPushMatrix();
    if ( dvel ) draw_temperature ();
    else		draw_density ();
    if (dAxis) draw_axis();
    glPopMatrix();
    
    glEnd();
    glPopMatrix();
    
    glEnd();
    glFlush();
    
    post_display();
}

static void special_func(int key, int x, int y) {
    
    switch(key)
    {
        case GLUT_KEY_UP:
            rot[0] += 5.0f;
            rot[0] = clamp(rot[0]);
            break;
        case GLUT_KEY_DOWN:
            rot[0] -= 5.0f;
            rot[0] = clamp(rot[0]);
            break;
        case GLUT_KEY_LEFT:
            rot[1] -= 5.0f;
            rot[1] = clamp(rot[1]);
            break;
        case GLUT_KEY_RIGHT:
            rot[1] += 5.0f;
            rot[1] = clamp(rot[1]);
            break;

    }
}

/*
 ----------------------------------------------------------------------
 open_glut_window --- open a glut compatible window and set callbacks
 ----------------------------------------------------------------------
 */

static void open_glut_window ( void )
{
    glutInitDisplayMode ( GLUT_RGBA | GLUT_DOUBLE );
    
    glutInitWindowPosition ( 0, 0 );
    glutInitWindowSize ( win_x, win_y );
    win_id = glutCreateWindow ( "3D Convection Fluids Simulation" );
    
    glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
    glClear ( GL_COLOR_BUFFER_BIT );
    glutSwapBuffers ();
    glClear ( GL_COLOR_BUFFER_BIT );
    glutSwapBuffers ();
    
    glutKeyboardFunc ( key_func );
    glutMouseFunc ( mouse_func );
    glutMotionFunc ( motion_func );
    glutReshapeFunc ( reshape_func );
    glutIdleFunc ( idle_func );
    glutDisplayFunc ( display_func );
    glutSpecialFunc(special_func);
}


/*
 ----------------------------------------------------------------------
 main --- main routine
 ----------------------------------------------------------------------
 */

int main ( int argc, char ** argv )
{
    glutInit ( &argc, argv );
    
    if ( argc != 1 && argc != 6 ) {
        fprintf ( stderr, "usage : %s N dt diff visc force source\n", argv[0] );
        fprintf ( stderr, "where:\n" );\
        fprintf ( stderr, "\t N      : grid resolution\n" );
        fprintf ( stderr, "\t dt     : time step\n" );
        fprintf ( stderr, "\t diff   : diffusion rate of the density\n" );
        fprintf ( stderr, "\t visc   : viscosity of the fluid\n" );
        fprintf ( stderr, "\t force  : scales the mouse movement that generate a force\n" );
        fprintf ( stderr, "\t source : amount of density that will be deposited\n" );
        exit ( 1 );
    }
    
    if ( argc == 1 ) {
        N = 32;
        dt = 0.4f;
        diff = 0.0001f;
        visc = 0.0f;
        force = 0.0f;
        source = 100.0f;
        temp_source = 35.0f;
        fprintf ( stderr, "Using defaults : N=%d dt=%g diff=%g visc=%g force = %g source=%g\n",
                 N, dt, diff, visc, force, source );
    } else {
        N = atoi(argv[1]);
        dt = atof(argv[2]);
        diff = atof(argv[3]);
        visc = atof(argv[4]);
        force = atof(argv[5]);
        source = atof(argv[6]);
        temp_source = atof(argv[7]);
    }
    
    printf ( "\n\nHow to use this demo:\n\n" );
    printf ( "\t Add densities with the right mouse button\n" );
    printf ( "\t Add velocities with the left mouse button and dragging the mouse\n" );
    printf ( "\t Toggle density/velocity display with the 'v' key\n" );
    printf ( "\t Clear the simulation by pressing the 'c' key\n" );
    printf ( "\t Quit by pressing the 'q' key\n" );
    
    dvel = 0;
    rot[0] = 30;
    rot[1] = -45;
    
    if ( !allocate_data () ) exit ( 1 );
    clear_data ();    
    
    win_x = 768;
    win_y = 768;
    open_glut_window ();
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    
    glEnable(GL_ALPHA_TEST);
    glAlphaFunc(GL_GREATER, 0);
    
    glutMainLoop ();
    
    exit ( 0 );
}