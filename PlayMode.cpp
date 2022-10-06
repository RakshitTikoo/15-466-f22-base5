#include "PlayMode.hpp"

#include "LitColorTextureProgram.hpp"

#include "DrawLines.hpp"
#include "Mesh.hpp"
#include "Load.hpp"
#include "gl_errors.hpp"
#include "data_path.hpp"

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/quaternion.hpp>

#include <random>

GLuint Vribe_meshes_for_lit_color_texture_program = 0;
Load< MeshBuffer > Vribe_meshes(LoadTagDefault, []() -> MeshBuffer const * {
	MeshBuffer const *ret = new MeshBuffer(data_path("Vribe.pnct"));
	Vribe_meshes_for_lit_color_texture_program = ret->make_vao_for_program(lit_color_texture_program->program);
	return ret;
});

Load< Scene > Vribe_scene(LoadTagDefault, []() -> Scene const * {
	return new Scene(data_path("Vribe.scene"), [&](Scene &scene, Scene::Transform *transform, std::string const &mesh_name){
		Mesh const &mesh = Vribe_meshes->lookup(mesh_name);

		scene.drawables.emplace_back(transform);
		Scene::Drawable &drawable = scene.drawables.back();

		drawable.pipeline = lit_color_texture_program_pipeline;

		drawable.pipeline.vao = Vribe_meshes_for_lit_color_texture_program;
		drawable.pipeline.type = mesh.type;
		drawable.pipeline.start = mesh.start;
		drawable.pipeline.count = mesh.count;

	});
});

// Add Walkmesh code

// *******************************************************************************************************************************************************************
// Lookup fails. Tried to print out meshes, but list is empty. The export-walkmeshes clearly creates the output file. Walkmesh Code Implemented. Game not implemented.
// *******************************************************************************************************************************************************************
WalkMesh const *walkmesh = nullptr;
Load< WalkMeshes > Vribe_walkmeshes(LoadTagDefault, []() -> WalkMeshes const * {
	WalkMeshes *ret = new WalkMeshes(data_path("Vribe1.w"));
	//walkmesh = &ret->lookup("PlatForm"); // Lookup not working
	walkmesh = &ret->first();
	return ret;
});


Load< Sound::Sample > Main_Song_1_sample(LoadTagDefault, []() -> Sound::Sample const * {
	return new Sound::Sample(data_path("Main5.opus"));
});
//
//Load< Sound::Sample > Main_Song_2_sample(LoadTagDefault, []() -> Sound::Sample const * {
//	return new Sound::Sample(data_path("Main2.opus"));
//});
//
//Load< Sound::Sample > Main_Song_3_sample(LoadTagDefault, []() -> Sound::Sample const * {
//	return new Sound::Sample(data_path("Main3.opus"));
//});
//
//Load< Sound::Sample > Main_Song_4_sample(LoadTagDefault, []() -> Sound::Sample const * {
//	return new Sound::Sample(data_path("Main4.opus"));
//});
//
//Load< Sound::Sample > Main_Song_5_sample(LoadTagDefault, []() -> Sound::Sample const * {
//	return new Sound::Sample(data_path("Main5.opus"));
//});
//
//
Load< Sound::Sample > Rain_sample(LoadTagDefault, []() -> Sound::Sample const * {
	return new Sound::Sample(data_path("Rain.opus"));
});

Load< Sound::Sample > Car_sample(LoadTagDefault, []() -> Sound::Sample const * {
	return new Sound::Sample(data_path("Car.opus"));
});



void PlayMode::get_transforms() {

	for(uint32_t i = 0; i < 42; i = i + 1) {
		obstacle[i] = nullptr;
	}

	//get pointers to leg for convenience:
	for (auto &transform : scene.transforms) {

		if (transform.name == "Player_Car") player.transform = &transform;
		else {
			std::string obstacle_str;
			for(uint32_t i = 0; i < 42; i = i + 1) {
				if(i < 10)
					obstacle_str = std::string("Obstacle.00") + std::to_string(i); 
				else 
					obstacle_str = std::string("Obstacle.0") + std::to_string(i); 

				if (transform.name == obstacle_str) obstacle[i] = &transform;
				
			}
		}		
	}




	if (player.transform == nullptr) throw std::runtime_error("player not found.");
	for(uint32_t i = 0; i < 42; i = i + 1) {
		if (obstacle[i] == nullptr) throw std::runtime_error("obstacle not found.");
	}

	// Update player walkpoint 
	player.at = walkmesh->nearest_walk_point(player.transform->position);

	//get pointer to camera for convenience:
	if (scene.cameras.size() != 1) throw std::runtime_error("Expecting scene to have exactly one camera, but it has " + std::to_string(scene.cameras.size()));
	camera = &scene.cameras.front();

	if (camera == nullptr) printf("No Camera found.");

	//Rain1_pos = Rain1->position;
	//Rain2_pos = Rain2->position;
	//Rain3_pos = Rain3->position;

	
}


void PlayMode::init(int i) {
	if(i == 0)
	 {
		get_transforms();

		//camera->transform->position.x = -0.0f; camera->transform->position.y = 144.520065f; camera->transform->position.z = 89.564598f;
		camera->transform->position.x = player.transform->position.x; camera->transform->position.y =  player.transform->position.y + 10.0f; camera->transform->position.z =  player.transform->position.z + 20.0f;
		camera->transform->rotation.w = 0.0f; camera->transform->rotation.x = 0.0f; camera->transform->rotation.y = -0.422618f; camera->transform->rotation.z = -0.906308f;

	 }

	rain_vol = 1.5f;
	car_vol = 0.0f;
	rain_rate = 20.0f;
	car_speed = 0.0f;
	max_speed = 60.0f;
	main1_vol = 0.0f;

}

PlayMode::PlayMode() : scene(*Vribe_scene) {
	
	init(0);
	
	//start playing music:
	RainLoop = Sound::loop(*Rain_sample,rain_vol, 0.0f);
	CarLoop = Sound::loop(*Car_sample,car_vol, 0.0f);
	Main1Loop  = Sound::loop(*Main_Song_1_sample,main1_vol, 0.0f);

	//player.transform = Player_Car;
	

}

PlayMode::~PlayMode() {
}

bool PlayMode::handle_event(SDL_Event const &evt, glm::uvec2 const &window_size) {

	if (evt.type == SDL_KEYDOWN) {
		if (evt.key.keysym.sym == SDLK_a) {
			left.downs += 1;
			left.pressed = true;
			return true;
		} else if (evt.key.keysym.sym == SDLK_d) {
			right.downs += 1;
			right.pressed = true;
			return true;
		} else if (evt.key.keysym.sym == SDLK_w) {
			up.downs += 1;
			up.pressed = true;
			return true;
		} else if (evt.key.keysym.sym == SDLK_s) {
			down.downs += 1;
			down.pressed = true;
			return true;
		}
	} else if (evt.type == SDL_KEYUP) {
		if (evt.key.keysym.sym == SDLK_a) {
			left.pressed = false;
			return true;
		} else if (evt.key.keysym.sym == SDLK_d) {
			right.pressed = false;
			return true;
		} else if (evt.key.keysym.sym == SDLK_w) {
			up.pressed = false;
			return true;
		} else if (evt.key.keysym.sym == SDLK_s) {
			down.pressed = false;
			return true;
		}
	} 

	return false;
}

void PlayMode::update(float elapsed) {

	// =======================
	// Rain Animation
	// =======================
	//Rain1->position.x -= (rain_rate*elapsed);
	//Rain1->position.z -= (rain_rate*elapsed);
	//Rain2->position.x -= (rain_rate*elapsed);
	//Rain2->position.z -= (rain_rate*elapsed);
	//Rain3->position.x -= (rain_rate*elapsed);
	//Rain3->position.z -= (rain_rate*elapsed);

	//if(Rain1->position.z < Side_Ground->position.z) Rain1->position = Rain1_pos;
	//if(Rain2->position.z < Side_Ground->position.z) Rain2->position = Rain2_pos;
	//if(Rain3->position.z < Side_Ground->position.z) Rain3->position = Rain3_pos;
	

	// ========================
	// Player Movement
	// ========================
	glm::vec2 move = glm::vec2(0.0f);
	if (left.pressed && !right.pressed) move.x =-1.0f;
	if (!left.pressed && right.pressed) move.x = 1.0f;
	move.y = 1;
	//if (down.pressed && !up.pressed) move.y =-1.0f;
	//if (!down.pressed && up.pressed) move.y = 1.0f;

	if(up.pressed == 1) {
		
		car_speed += 2.0f;
		car_vol += 0.02f;

		if(car_speed > max_speed) car_speed = max_speed;
		if(car_vol > 1.0f) car_vol = 1.0f;

		if(left.pressed == 1) {}
		if(right.pressed == 1) {}
	}
	if (down.pressed == 1) {
		car_speed -= 3.0f;
		car_vol -= 0.03f;

		if(car_speed < -1.0f*max_speed) car_speed = -1.0f*max_speed;
		if(car_vol < 0.0f) car_vol = 0.0f;
	}

	if(up.pressed == 0 && down.pressed == 0) {
		car_speed -= 2.0f;
		car_vol -= 0.02f;

		if(car_speed < 0.0f) car_speed = 0.0f;
		if(car_vol < 0.0f) car_vol = 0.0f;
	}

	if(car_speed > 0.0f) main1_vol += 0.2f;
	else main1_vol -= 0.2f;

	if(main1_vol > 10.0f) main1_vol = 10.0f;
	if(main1_vol < 0.0f) main1_vol = 0.0f;

	move.x = 0.1f*max_speed*move.x*elapsed;
	move.y = car_speed*move.y*elapsed;

	//glm::vec3 remain = player.transform->make_local_to_world() * glm::vec4(move.x, move.y, 0.0f, 0.0f);


	// =====================
	// Enemy Collision
	// =====================



	// =======================
	// Volume Control
	// =======================

	CarLoop->volume = car_vol;
	Main1Loop->volume = main1_vol;




		//move camera:
	//{
//
	//	//combine inputs into a move:
	//	constexpr float PlayerSpeed = 30.0f;
	//	glm::vec2 move = glm::vec2(0.0f);
	//	if (left.pressed && !right.pressed) move.x =-1.0f;
	//	if (!left.pressed && right.pressed) move.x = 1.0f;
	//	if (down.pressed && !up.pressed) move.y =-1.0f;
	//	if (!down.pressed && up.pressed) move.y = 1.0f;
//
	//	//make it so that moving diagonally doesn't go faster:
	//	if (move != glm::vec2(0.0f)) move = glm::normalize(move) * PlayerSpeed * elapsed;
//
	//	glm::mat4x3 frame = camera->transform->make_local_to_parent();
	//	glm::vec3 frame_right = frame[0];
	//	//glm::vec3 up = frame[1];
	//	glm::vec3 frame_forward = -frame[2];
//
	//	camera->transform->position += move.x * frame_right + move.y * frame_forward;
	//}
	
	
	//player walking:
	//{
		//combine inputs into a move:
		//constexpr float PlayerSpeed = 3.0f;
		//glm::vec2 move = glm::vec2(0.0f);
		//if (left.pressed && !right.pressed) move.x =-1.0f;
		//if (!left.pressed && right.pressed) move.x = 1.0f;
		//if (down.pressed && !up.pressed) move.y =-1.0f;
		//if (!down.pressed && up.pressed) move.y = 1.0f;
//
		//////make it so that moving diagonally doesn't go faster:
		if (move != glm::vec2(0.0f)) move = glm::normalize(move);

		//get move in world coordinate system:
		glm::vec3 remain = player.transform->make_local_to_world() * glm::vec4(move.x, move.y, 0.0f, 0.0f);

		//using a for() instead of a while() here so that if walkpoint gets stuck in
		// some awkward case, code will not infinite loop:
		for (uint32_t iter = 0; iter < 10; ++iter) {
			if (remain == glm::vec3(0.0f)) break;
			WalkPoint end;
			float time;
			walkmesh->walk_in_triangle(player.at, remain, &end, &time);
			player.at = end;
			if (time == 1.0f) {
				//finished within triangle:
				remain = glm::vec3(0.0f);
				break;
			}
			//some step remains:
			remain *= (1.0f - time);
			//try to step over edge:
			glm::quat rotation;
			if (walkmesh->cross_edge(player.at, &end, &rotation)) {
				//stepped to a new triangle:
				player.at = end;
				//rotate step to follow surface:
				remain = rotation * remain;
			} else {
				//ran into a wall, bounce / slide along it:
				glm::vec3 const &a = walkmesh->vertices[player.at.indices.x];
				glm::vec3 const &b = walkmesh->vertices[player.at.indices.y];
				glm::vec3 const &c = walkmesh->vertices[player.at.indices.z];
				glm::vec3 along = glm::normalize(b-a);
				glm::vec3 normal = glm::normalize(glm::cross(b-a, c-a));
				glm::vec3 in = glm::cross(normal, along);
//
				//check how much 'remain' is pointing out of the triangle:
				float d = glm::dot(remain, in);
				if (d < 0.0f) {
					//bounce off of the wall:
					remain += (-1.25f * d) * in;
				} else {
					//if it's just pointing along the edge, bend slightly away from wall:
					remain += 0.01f * d * in;
				}
			}
		}
//
		if (remain != glm::vec3(0.0f)) {
			std::cout << "NOTE: code used full iteration budget for walking." << std::endl;
		}
//
		//update player's position to respect walking:
		player.transform->position = walkmesh->to_world_point(player.at);
//
		{ //update player's rotation to respect local (smooth) up-vector:
			
			glm::quat adjust = glm::rotation(
				player.transform->rotation * glm::vec3(0.0f, 0.0f, 1.0f), //current up vector
				walkmesh->to_world_smooth_normal(player.at) //smoothed up vector at walk location
			);
			player.transform->rotation = glm::normalize(adjust * player.transform->rotation);
		}
//
		
		//glm::mat4x3 frame = camera->transform->make_local_to_parent();
		//glm::vec3 right = frame[0];
		////glm::vec3 up = frame[1];
		//glm::vec3 forward = -frame[2];

		//camera->transform->position += move.x * right + move.y * forward;
		
	//}

	//reset button press counters:
	left.downs = 0;
	right.downs = 0;
	up.downs = 0;
	down.downs = 0;
}

void PlayMode::draw(glm::uvec2 const &drawable_size) {
	//update camera aspect ratio for drawable:
	camera->aspect = float(drawable_size.x) / float(drawable_size.y);

	//set up light type and position for lit_color_texture_program:
	// TODO: consider using the Light(s) in the scene to do this
	glUseProgram(lit_color_texture_program->program);
	glUniform1i(lit_color_texture_program->LIGHT_TYPE_int, 1);
	glUniform3fv(lit_color_texture_program->LIGHT_DIRECTION_vec3, 1, glm::value_ptr(glm::vec3(0.0f, 0.0f,-1.0f)));
	glUniform3fv(lit_color_texture_program->LIGHT_ENERGY_vec3, 1, glm::value_ptr(glm::vec3(1.0f, 1.0f, 0.95f)));
	glUseProgram(0);

	glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
	glClearDepth(1.0f); //1.0 is actually the default value to clear the depth buffer to, but FYI you can change it.
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS); //this is the default depth comparison function, but FYI you can change it.

	scene.draw(*camera);

	/* In case you are wondering if your walkmesh is lining up with your scene, try:
	{
		glDisable(GL_DEPTH_TEST);
		DrawLines lines(player.camera->make_projection() * glm::mat4(player.camera->transform->make_world_to_local()));
		for (auto const &tri : walkmesh->triangles) {
			lines.draw(walkmesh->vertices[tri.x], walkmesh->vertices[tri.y], glm::u8vec4(0x88, 0x00, 0xff, 0xff));
			lines.draw(walkmesh->vertices[tri.y], walkmesh->vertices[tri.z], glm::u8vec4(0x88, 0x00, 0xff, 0xff));
			lines.draw(walkmesh->vertices[tri.z], walkmesh->vertices[tri.x], glm::u8vec4(0x88, 0x00, 0xff, 0xff));
		}
	}
	*/

	{ //use DrawLines to overlay some text:
		glDisable(GL_DEPTH_TEST);
		float aspect = float(drawable_size.x) / float(drawable_size.y);
		DrawLines lines(glm::mat4(
			1.0f / aspect, 0.0f, 0.0f, 0.0f,
			0.0f, 1.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 1.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f
		));

		constexpr float H = 0.09f;
		lines.draw_text("Mouse motion looks; WASD moves; escape ungrabs mouse",
			glm::vec3(-aspect + 0.1f * H, -1.0 + 0.1f * H, 0.0),
			glm::vec3(H, 0.0f, 0.0f), glm::vec3(0.0f, H, 0.0f),
			glm::u8vec4(0x00, 0x00, 0x00, 0x00));
		float ofs = 2.0f / drawable_size.y;
		lines.draw_text("Mouse motion looks; WASD moves; escape ungrabs mouse",
			glm::vec3(-aspect + 0.1f * H + ofs, -1.0 + + 0.1f * H + ofs, 0.0),
			glm::vec3(H, 0.0f, 0.0f), glm::vec3(0.0f, H, 0.0f),
			glm::u8vec4(0xff, 0xff, 0xff, 0x00));
	}
	GL_ERRORS();
}
