#include <stdlib.h>
#include <stdio.h>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <pthread.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "./2d_grid_drawer/stb_image_write.h"
#include <stdlib.h>

void handle_close(GLFWwindow * window);
void window_size_callback(GLFWwindow * window, int xsize, int ysize);
double*generate_test_array();
unsigned char* normalize_array(double * array, int array_size);
void diffuse(double * array, int width, int height, int passes);

void handle_close(GLFWwindow * window) {
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS){
		glfwSetWindowShouldClose(window, 1);
	}
}

void window_size_callback(GLFWwindow * window, int xsize, int ysize){
	glViewport(0, 0, xsize, ysize);
}

double rand_from_range(double min, double max){
	double range = max-min;
	double div = RAND_MAX/range;
	return min+(rand()/div);
}

double*generate_test_array(){
	int size = 1920*1080;
	double *test_array = malloc(sizeof(double)*size);
	srand(time(NULL));
	
	#pragma omp parallel for simd
	for (int i = 0; i < size; i++){
		test_array[i] = pow(rand_from_range(-1.0, 1.0),2);
	}
	
	return test_array;
}

unsigned char* normalize_array(double * array, int array_size){
	int draw_size = array_size*3;
	unsigned char* drawable_array = malloc(sizeof(unsigned char)*draw_size);
	
	double min = array[0];
	double max= array[0];
	
	#pragma omp parallel
	{
		double local_min = array[0];
		double local_max= array[0];
		for (int i = 0; i < array_size; i++){
			if (array[i] > local_max) {
				local_max = array[i];
			} else if (array[i] < local_min) {
				local_min= array[i];
			}
			
		}
		#pragma omp critial 
		{
			if (local_min < min) {
				min = local_min;
			} else if (local_max > max) {
				max = local_max;
			}
		}
	}

	double range = (max - min);
	#pragma omp parallel for simd
	for (int i = 0; i < array_size; i ++) {
		unsigned char blue = (unsigned char)255*(array[i]- min)/range;
		drawable_array[3*i] = blue;	
		drawable_array[3*i+2] = 255-blue;	
	}
	
	return drawable_array;
}

void diffuse(double * array, int width, int height, int passes) {

	double * temp = malloc(sizeof(double)*width*height);

	for (int p = 0; p < passes; p++){
		
		memcpy(temp, array, sizeof(double)*width*height);

		#pragma omp parallel for
		for (int i = 1; i < height-1; i ++) {
			#pragma omp parallel for simd
			for (int j = 1; j < width-1; j++){
				array[j+width*i] = (0.6*temp[j+width*i] + 0.1*temp[j+width*i+1] + 0.1*temp[j+width*i-1] + 0.1*temp[j+(width)*(i-1)] + 0.1*temp[j+(width)*(i+1)]);
			}
		}
	}

	free(temp);

}

typedef struct drawstruct{
	const char * print_name;
	int width;
	int height;
	unsigned int shader_program;
	unsigned int texture; 
	unsigned int VAO;
	GLFWwindow* window;
	unsigned char * data;
} drawstruct;

void * initialize_glfw(void * arg) {

	drawstruct * drawable = (drawstruct *)(arg);

	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);	

	drawable->window = glfwCreateWindow(drawable->width,drawable->height, "grid", NULL, NULL);

	glfwMakeContextCurrent(drawable->window);
	glfwSetFramebufferSizeCallback(drawable->window, window_size_callback);
	glfwPollEvents();

	glewExperimental = GL_TRUE;
	glewInit();

	glViewport(0,0,1920, 1080);


	float verts [] = {
		-0.90, -0.90,0,0,
		0.90,-0.90,0.90,0,
		-0.90,0.90,0,0.90,
		0.90,-0.90,0.90,0,
		0.90,0.90,0.90,0.90,
		-0.90,0.90,0,0.90,
	};

	
	glGenVertexArrays(1, &drawable->VAO);
	glBindVertexArray(drawable->VAO);

	unsigned int VBO;
	glGenBuffers(1,&VBO);
	glBindBuffer(GL_ARRAY_BUFFER,VBO);
	glBufferData(GL_ARRAY_BUFFER,sizeof(verts), verts, GL_STATIC_DRAW);

	int success;
	const char * vertex_shader_source = "#version 330 core\n"
	"layout (location = 0) in vec2 aPos;\n"
	"layout (location = 1) in vec2 aTexCoord;\n"
	"out vec2 TexCoord;\n"
	"void main()\n"
	"{\n"
	"gl_Position = vec4(aPos.x, aPos.y, 0.0, 1.0);\n"
	"TexCoord = aTexCoord;"
	"}\0";

	const char * fragment_shader_source = "#version 330 core\n"
	"out vec4 FragColor;\n"
	"in vec2 TexCoord;\n"
	"uniform sampler2D ourTexture;\n"
	"void main()\n"
	"{\n"
	"FragColor = texture(ourTexture,TexCoord);\n"
	"}\0";

	unsigned int vertex_shader = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vertex_shader, 1, &vertex_shader_source, NULL);
	glCompileShader(vertex_shader);
	glGetShaderiv(vertex_shader, GL_COMPILE_STATUS, &success);
	assert(success);

	unsigned int fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fragment_shader, 1, &fragment_shader_source, NULL);
	glCompileShader(fragment_shader);
	glGetShaderiv(fragment_shader, GL_COMPILE_STATUS, &success);

	drawable->shader_program= glCreateProgram();
	glAttachShader(drawable->shader_program, fragment_shader);
	glAttachShader(drawable->shader_program, vertex_shader);
	glLinkProgram(drawable->shader_program);
	glGetProgramiv(drawable->shader_program, GL_LINK_STATUS, &success);
	assert(success);

	glUseProgram(drawable->shader_program);
	glDeleteShader(vertex_shader);
	glDeleteShader(fragment_shader);

	glVertexAttribPointer(0,2, GL_FLOAT, GL_FALSE, 4*sizeof(float), (void *) 0);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4*sizeof(float), (void *)(2*sizeof(float)));
	glEnableVertexAttribArray(0);		
	glEnableVertexAttribArray(1);

	glGenTextures(1, &drawable->texture);
	glBindTexture(GL_TEXTURE_2D, drawable->texture);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	return NULL;
}

void * display_image(void * arg){

	drawstruct * drawable = arg;
	int width = drawable->width;
	int height = drawable -> height;
	unsigned int shader_program = drawable->shader_program;
	unsigned int texture = drawable->texture; 
	unsigned int VAO = drawable->VAO;
	GLFWwindow * window = drawable ->window;
	unsigned char* data = drawable->data;
		
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB ,GL_UNSIGNED_BYTE, data);
	glGenerateMipmap(GL_TEXTURE_2D);

	while(!glfwWindowShouldClose(window)) {
		handle_close(window);
		glClearColor(1.0, 1.0, 1.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, texture);

		glBindVertexArray(VAO);
		glUseProgram(shader_program);
		glDrawArrays(GL_TRIANGLES, 0, 6);

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	glfwTerminate();

	return NULL;
}

void *print_image(void * arg){
	drawstruct * drawable = (drawstruct *)arg;
	int width = drawable->width;
	int height = drawable -> height;
	unsigned char *data = drawable->data;
	stbi_write_jpg(drawable->print_name, width, height,3,data, 90);
	return NULL;
} 

typedef struct normalization_data{
	int width;
	int height; 
	int diffusion_steps;
	double * data;
	unsigned char * normalized_data;
	int free_data;
}normalization_data;

void * work_thread_work(void * arg){
	normalization_data * local_data = arg;
	if (local_data->free_data== 0) {
		double * temp_data = malloc(sizeof(double)*local_data->width*local_data->height);
		memcpy(temp_data, local_data->data, sizeof(double)*local_data->width*local_data->height);
		local_data->data= temp_data;
	}
	 
	//quesion: can I get this into a thread and run it imediately?
	diffuse(local_data->data,local_data->width,local_data->height,local_data->diffusion_steps);
	local_data-> normalized_data = normalize_array(local_data->data,local_data->width*local_data->height);
	return NULL;
}

void draw_heatmap(double * data, int width, int height, int free_data, int diffusion_steps, int show, int print, char * printname){

	pthread_t work_thread;
	pthread_t print_thread;

	drawstruct *drawable = malloc(sizeof(drawstruct));
	drawable->width = width;
	drawable->height = height;
	drawable->print_name = printname;

	normalization_data *normalization_data_struct= malloc(sizeof(normalization_data));
	normalization_data_struct->width = width;
	normalization_data_struct->height = height;
	normalization_data_struct->diffusion_steps= diffusion_steps;
	normalization_data_struct->free_data=free_data;
	normalization_data_struct->data = data;

	pthread_create(&work_thread, NULL, work_thread_work,normalization_data_struct);

	if (show == 1){
		initialize_glfw(drawable);
	}

	pthread_join(work_thread,NULL);
	
	if (print == 1) {
		drawable->print_name = printname;
		pthread_create(&print_thread, NULL, print_image, drawable);
	}	

	if (show == 1){
		drawable->data = normalization_data_struct->normalized_data;
		display_image(drawable);
	}

	if (print == 1) {
		pthread_join(print_thread,NULL);
	}
	
	free(drawable);
	free(normalization_data_struct->normalized_data);
	free(normalization_data_struct->data);
	free(normalization_data_struct);
}   

//-----------------------------------------------------------------------------------------------------------------------

typedef struct grid_point {
	double h;
	double hu;
	double hv;
} grid_point;

typedef struct grid{
	int x_size;
	int y_size;
	grid_point * points_current;
	grid_point * points_next;
	grid_point * points_temp;
	int nr_timesteps;
	double * times; 
} grid;

void set_grid_point(grid_point * point, double h, double hu, double hv){
	point->h = h;
	point->hu = hu;
	point->hv = hv;
}

grid * create_grid(int x_size, int y_size){
	grid * new_grid = calloc(1, sizeof(grid));
	new_grid -> x_size = x_size;
	new_grid->y_size = y_size;
	new_grid->points_current = calloc(x_size*y_size, sizeof(grid_point));
	new_grid->points_next= calloc(x_size*y_size, sizeof(grid_point));
	new_grid->points_temp = 0;
	new_grid->nr_timesteps = 0;
	new_grid->times = NULL;
	return new_grid;
}

double * printable_grid_generator(grid * grid){
	double * heights = calloc(grid->x_size*grid->y_size,sizeof(double));

	#pragma omp parallel for
	for (int h = 0; h<grid->y_size; h++){
		#pragma omp parallel for simd
		for (int w = 0; w<grid->x_size; w++){
			heights[w+grid->x_size*h] = grid->points_current[w+grid->x_size*h].h;
		}
	}
	
	return heights;

}

void grid_points_copy(int x_size, int y_size,grid_point * source, grid_point * destination){
	#pragma omp parallel for
	for (int h = 0; h<y_size; h++){
		#pragma omp parallel for simd
		for (int w = 0; w<x_size; w++){
			destination[w+x_size*h] = source[w+x_size*h];
		}
	}
}

void set_initial_conditions(grid * grid){
	//goal: no momentum but increased height in the middle
	int x_size = grid->x_size;
	int y_size = grid->y_size;

	#pragma omp parallel 
	{	
		#pragma omp for
		for (int h = 0; h<y_size; h++){
			#pragma omp parallel for simd
			for (int w = 0; w<x_size; w++){
				grid_point setter = {1.0, 0.0,0.0};
				grid->points_current[w+x_size*h] = setter;
				grid->points_next[w+x_size*h] = setter; 
			}
		}
	
		#pragma omp barrier
		
		#pragma omp for
		for (int h = y_size/4; h<3*y_size/4; h++){
			#pragma omp parallel for simd
			for (int w = x_size/4; w<3*x_size/4; w++){
				grid_point setter = {2.0, 0.0,0.0};
				grid->points_current[w+x_size*h] = setter;
				grid->points_next[w+x_size*h] = setter; 
			}
		}
	}
}

void diffusion_step(grid * grid){
	//this should be doable
	int x_size = grid->x_size;
	int y_size = grid->y_size;
	grid_point * current = grid->points_current;
	grid_point * next = grid->points_next;
		
	#pragma omp for
	for (int h = 1; h<y_size-1; h++){
		#pragma omp parallel for simd
		for (int w = 1; w<x_size-1; w++){
			next[w+x_size*h].h = 0.8*current[w+x_size*h].h + 0.05*current[(w+1)+x_size*h].h+ 0.05*current[(w-1)+x_size*h].h  + 0.05*current[w+x_size*(h-1)].h+ 0.05*current[w+x_size*(h+1)].h;
			next[w+x_size*h].hu = 0.8*current[w+x_size*h].hu + 0.05*current[(w+1)+x_size*h].hu+ 0.05*current[(w-1)+x_size*h].hu  + 0.05*current[w+x_size*(h-1)].hu+ 0.05*current[w+x_size*(h+1)].hu;
			next[w+x_size*h].hv = 0.8*current[w+x_size*h].hv + 0.05*current[(w+1)+x_size*h].hv+ 0.05*current[(w-1)+x_size*h].hv  + 0.05*current[w+x_size*(h-1)].hv+ 0.05*current[w+x_size*(h+1)].hv;
		}
	}
		
}

void swap_next_and_current(grid * grid) {
	grid->points_temp = grid->points_next;
	grid->points_next = grid->points_current;
	grid->points_current = grid->points_temp;
}

void multi_diffusion(grid * grid, int steps){
	for (int i = 0; i< steps; i++){
		diffusion_step(grid);
		swap_next_and_current(grid);
	}
}

double f_2(grid_point current) {
	return ((current.hu*current.hu)/current.h) + 0.5*10*(current.h*current.h);
}

double f_3(grid_point current) {
	return ((current.hu*current.hv)/current.h);
}

void convection_step_right(grid * grid){
	//test that everything else is working before implementing this!
	int x_size = grid->x_size;
	int y_size = grid->y_size;
	grid_point * current = grid->points_current;
	grid_point * next = grid->points_next;
		
	#pragma omp for
	for (int h = 1; h<y_size-1; h++){
		#pragma omp parallel for simd
		for (int w = 1; w<x_size-1; w++){


			next[w+x_size*h].h = 0.5*(current[(w+1)+x_size*h].h + current[(w-1)+x_size*h].h) - (0.1/(4*1))*(current[(w-1)+x_size*h].hu - current[(w-1)+x_size*h].hu);
			next[w+x_size*h].hu = 0.5*(current[(w+1)+x_size*h].hu + current[(w-1)+x_size*h].hu) - (0.1/(4*1))*(f_2(current[(w+1)+x_size*h]) - f_2(current[(w-1)+x_size*h]));
			next[w+x_size*h].hv = 0.5*(current[(w+1)+x_size*h].hv + current[(w-1)+x_size*h].hv) - (0.1/(4*1))*(f_3(current[(w+1)+x_size*h]) - f_3(current[(w-1)+x_size*h]));


		}
	}
}

double f_4(grid_point current) {
	return ((current.hv*current.hv)/current.h) + 0.5*10*(current.h*current.h);
}

void convection_step_north(grid * grid){
	//test that everything else is working before implementing this!
	int x_size = grid->x_size;
	int y_size = grid->y_size;
	grid_point * current = grid->points_current;
	grid_point * next = grid->points_next;
		
	#pragma omp for
	for (int h = 1; h<y_size-1; h++){
		#pragma omp parallel for simd
		for (int w = 1; w<x_size-1; w++){


			next[w+x_size*h].h = 0.5*(current[w+x_size*(h-1)].h+ current[w+x_size*(h+1)].h) - (0.1/(4))*(current[w+x_size*(h-1)].hv - current[w+x_size*(h+1)].hv);
			next[w+x_size*h].hu = 0.5*(current[w+x_size*(h-1)].hu+ current[w+x_size*(h+1)].hu)- (0.1/(4))* (f_3(current[w+x_size*(h-1)]) - f_3(current[w+x_size*(h+1)]));
			next[w+x_size*h].hv = 0.5*(current[w+x_size*(h-1)].hv+ current[w+x_size*(h+1)].hv) - (0.1/(4))* (f_4(current[w+x_size*(h-1)]) - f_4(current[w+x_size*(h+1)]));


		}
	}
}

void adjust_boundary(grid * grid) {
	#pragma omp parallel for simd
	for (int i = 0; i < grid->x_size; i++){
		grid->points_current[i]=grid->points_current[i+grid->x_size];
		grid->points_current[grid->x_size*grid->y_size-i]=grid->points_current[grid->x_size*(grid->y_size-1)-i];
	}
	
	#pragma omp parallel for simd
	for (int i = 0; i < grid->y_size; i++){
		grid->points_current[i*grid->x_size]=grid->points_current[i*grid->x_size+1];
		grid->points_current[(i+1)*grid->x_size-1]=grid->points_current[(i+1)*grid->x_size-2];
	}
}

void multi_convection(grid * grid, int steps){
	for (int i = 0; i< steps; i++){
		convection_step_right(grid);
		swap_next_and_current(grid);
		adjust_boundary(grid);
		convection_step_north(grid);
		swap_next_and_current(grid);
		adjust_boundary(grid);
		diffusion_step(grid);
		swap_next_and_current(grid);
		adjust_boundary(grid);
	}
}

void simulate(){
	grid * simulation_grid = create_grid(1000,1000);
	set_initial_conditions(simulation_grid);
	multi_convection(simulation_grid,100000);
	double * printable = printable_grid_generator(simulation_grid);
	draw_heatmap(printable, 1000, 1000, 1, 0, 1, 0, NULL);
}

int main(){
	simulate();	
}
