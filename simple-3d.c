#include <math.h>
#include <unistd.h>
#include <SDL.h>

#define SDL_MAIN_HANDLED

typedef struct{
  double node[3];
} node;

typedef struct{
  size_t connection[2];
} conn;

typedef struct{
  node* nodes;
  conn* links;
  size_t n_nodes;
  size_t n_links;
} structure;

structure read_structre_from_files(FILE *nodes_f, FILE* links_f) {

  size_t i, n_nodes, n_links;

  structure s;
  node* n;
  conn* l;

  n = (node*)malloc(sizeof(node)*1);

  n_nodes = 0;
  while(3 == fscanf(nodes_f,
		    "%lf %lf %lf",
		    n[n_nodes].node,
		    n[n_nodes].node+1,
		    n[n_nodes].node+2)) {
    n_nodes++;
    n = (node*)realloc(n,sizeof(node)*(n_nodes+1));
  }

  n = (node*)realloc(n,sizeof(node)*n_nodes);

  l = (conn*)malloc(sizeof(conn)*1);

  n_links = 0;
  while(2 == fscanf(links_f,
		    "%lu %lu",
		    l[n_links].connection,
		    l[n_links].connection+1)) {
    if(l[n_links].connection[0] < 0 || l[n_links].connection[0] > n_nodes ||
       l[n_links].connection[1] < 0 || l[n_links].connection[1] > n_nodes) {
      
      fprintf(stderr,"Links file corrupt! Exiting!\n");
      _exit(1);
    }
    n_links++;
    l = (conn*)realloc(l,sizeof(conn)*(n_links+1));
  }

  s.links = l;
  s.nodes = n;
  s.n_nodes = n_nodes;
  s.n_links = n_links;

  return(s);
}

void rescale_structure(size_t imagesize, structure* s) {

  size_t i,node;
  double center[3];
  double max_dist=0.;
  
  double d;

  double scalefactor;

  double imagesize_d = (double)imagesize;
  double imagesize_d_half = imagesize_d / 2.;
  
  memset(center,0,sizeof(double)*3);
  
  /* find center of the structure */
  for(node=0;node<s->n_nodes;node++) {
    for(i=0;i<3;i++) {
      center[i]+=s->nodes[node].node[i];
    }
  }
  for(i=0;i<3;i++) {
    center[i]/=s->n_nodes;
  }

  /* find maximum distance from center in one direction */
  for(node=0;node<s->n_nodes;node++) {
    d =0;
    for(i=0;i<3;i++) {
      d +=
      (center[i]-s->nodes[node].node[i])*(center[i]-s->nodes[node].node[i]);
    }
    d = sqrt(d);
    if (d > max_dist) {
      max_dist = d;

    }
  }

  scalefactor = (imagesize_d_half-imagesize_d_half*0.3)/max_dist;

  /* center nodes to 0 and scale*/
  for(node=0;node<s->n_nodes;node++) {
    for(i=0;i<3;i++) {
      s->nodes[node].node[i] = s->nodes[node].node[i]-center[i];
      s->nodes[node].node[i] = s->nodes[node].node[i]*scalefactor;
    }
  }
}

void rotate(double* c, double* a) {

  size_t i;
  double coord[3];

  double c_a;
  double s_a;

  c_a = cos(a[0]);
  s_a = sin(a[0]);
  
  coord[0] = c[0];
  coord[1] = c_a*c[1]-s_a*c[2];
  coord[2] = s_a*c[1]+c_a*c[2];

  for(i=0;i<3;i++) {
    c[i] = coord[i];
  }

  c_a = cos(a[1]);
  s_a = sin(a[1]);

  coord[0] = c_a*c[0]+s_a*c[2];
  coord[1] = c[1];
  coord[2] = -s_a*c[0]+c_a*c[2];

  for(i=0;i<3;i++) {
    c[i] = coord[i];
  }

  c_a = cos(a[2]);
  s_a = sin(a[2]);
  
  coord[0] = c_a*c[0]-s_a*c[1];
  coord[1] = s_a*c[0]+c_a*c[1];
  coord[2] = c[2];

  for(i=0;i<3;i++) {
    c[i] = coord[i];
  }
}

void uv_from_coordinate(size_t* uv, double* c,size_t imagesize) {

  double imagesize_d = (double) imagesize;
  double imagesize_d_half = imagesize_d/2.;
  
  double coord[3];

  long int u_l;
  long int v_l;

  double z_retraction = 0.00130*(c[2]+imagesize_d);
  
  u_l = c[0]/z_retraction+imagesize_d_half;
  v_l = c[1]/z_retraction+imagesize_d_half;

  //  u_l /= z_retraction;
  //  v_l /= z_retraction;

  uv[0] = (size_t) u_l;
  uv[1] = (size_t) v_l;
}

void drawpixel(size_t x, size_t y, unsigned int* image, size_t imagesize) {

  unsigned int* pixel_i;
  unsigned char* pixel_c;

  pixel_i = image+x*imagesize+y;

  pixel_c = (char*) pixel_i;

  pixel_c[0] = 255;
  pixel_c[1] = 255;
  pixel_c[2] = 255;
  pixel_c[3] = 0;
  
}

void bresenham_line(unsigned int* image, int imagesize,
		   size_t* uv_start,
		   size_t* uv_end) {


  int x,y;
  int prev_x, prev_y;
  
  int dx, sx;
  int dy, sy;

  int err;
  int err_times_two;

  x = uv_end[0];
  y = uv_end[1];
  prev_x = uv_start[0];
  prev_y = uv_start[1];
  
  dx = abs(x - prev_x);
  dy = -abs(y - prev_y);

  sx = x < prev_x ? 1 : -1;
  sy = y < prev_y ? 1 : -1;

  err = dx+dy;
  err_times_two = 2*err;

  while( x != prev_x || y != prev_y) {
    drawpixel(x,y,image,imagesize);
    err_times_two = 2*err;
    if(err_times_two > dy) {err += dy; x+=sx;}
    if(err_times_two < dx) {err += dx; y+=sy;}
  }
}
 
void draw_wireframe_at_rotation(structure s,
				size_t imagesize, unsigned int* image,
				double* angle) {

  size_t i;
  size_t j;
  size_t conn;

  double coord_a[3];
  double coord_b[3];

  size_t uv_a[2];
  size_t uv_b[2];
  
  double coord_transformed[3];
  
  for(conn = 0;conn<s.n_links;conn++) {
    for(i = 0; i< 3; i++) {
      coord_a[i] = s.nodes[s.links[conn].connection[0]].node[i];
      coord_b[i] = s.nodes[s.links[conn].connection[1]].node[i];
    }
    rotate(coord_a,angle);
    rotate(coord_b,angle);
    
    uv_from_coordinate(uv_a,coord_a,imagesize);
    uv_from_coordinate(uv_b,coord_b,imagesize);
    bresenham_line(image,imagesize,uv_a,uv_b);
  }
}

int main(int argc, char** argv) {

  SDL_Window * my_window = NULL;
  SDL_Renderer * my_renderer = NULL;

  SDL_Texture * my_texture;
  SDL_Event my_event;

  int imagesize = 600;
  int pitch = imagesize*4;

  double angles[3] = {0.,0.,0.};

  unsigned int *image =
    (unsigned int*)malloc(sizeof(unsigned int)*imagesize*imagesize);

  structure s;

  FILE* node_f = fopen(argv[1],"r");
  FILE* links_f = fopen(argv[2],"r");

  int update=1;
  
  s = read_structre_from_files(node_f, links_f);

  rescale_structure(imagesize,&s);
  
  SDL_SetMainReady();

  SDL_Init(SDL_INIT_VIDEO | SDL_INIT_EVENTS);

  my_window = SDL_CreateWindow("3D Wireframe Viewer", 100,100,
			       imagesize, imagesize,0);

  my_renderer = SDL_CreateRenderer(my_window,-1,SDL_RENDERER_ACCELERATED);

  my_texture = SDL_CreateTexture(my_renderer, SDL_PIXELFORMAT_ARGB8888,
				 SDL_TEXTUREACCESS_STREAMING,
				 imagesize,imagesize);

  while(1) {

    if(update) {

      SDL_LockTexture(my_texture,NULL,(void**)&image,&pitch);
      image = memset(image,0,sizeof(unsigned int)*imagesize*imagesize);
      draw_wireframe_at_rotation(s,
				 (size_t)imagesize, image,
				 angles);
      SDL_UnlockTexture(my_texture);
      
      SDL_RenderCopy(my_renderer,my_texture,NULL,NULL);
      SDL_RenderPresent(my_renderer);

      update = 0;
      
    }
    
    
    while(SDL_PollEvent(&my_event)) {
      switch(my_event.type) {
      case SDL_QUIT:
	goto finish;
	break;
      case SDL_WINDOWEVENT:
	if(my_event.window.event == SDL_WINDOWEVENT_CLOSE) {
	  SDL_Quit();
	  goto finish;
	}
	break;
      case SDL_KEYDOWN:
	if(my_event.key.keysym.sym == SDLK_q) {
	  angles[0]+=0.05;
	  if(angles[0] > 6.283) {
	    angles[0] = 0;
	  }
	  update = 1;
	}
	if(my_event.key.keysym.sym == SDLK_e) {
	 angles[0]-=0.05;
	  if(angles[0] < 0) {
	    angles[0] = 6.283;
	  }
	  update = 1;
	} 
	if(my_event.key.keysym.sym == SDLK_a) {
	  angles[1]+=0.05;
	  if(angles[1] > 6.283) {
	    angles[1] = 0.;
	  }
	  update = 1;
	}
	if(my_event.key.keysym.sym == SDLK_d) {
	  angles[1]-=0.05;
	 if(angles[1] < 0) {
	   angles[1] = 6.283;
	 }
	 update = 1;
	}
	if(my_event.key.keysym.sym == SDLK_z) {
	  angles[2]+=0.05;
	  if(angles[2] > 6.283) {
	    angles[2] = 0.;
	  }
	  update = 1;
	}
	if(my_event.key.keysym.sym == SDLK_c) {
	  angles[2]-=0.05;
	 if(angles[2] < 0) {
	   angles[2] = 6.283;
	 }
	 update = 1;
	}
	
      }
    }
    SDL_Delay(5);
  }

 finish:
  return(0);
}
  
	
