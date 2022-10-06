#include "Mode.hpp"

#include "Scene.hpp"
#include "WalkMesh.hpp"
#include "Sound.hpp"
#include <glm/glm.hpp>

#include <vector>
#include <deque>

struct PlayMode : Mode {
	PlayMode();
	virtual ~PlayMode();

	//functions called by main loop:
	virtual bool handle_event(SDL_Event const &, glm::uvec2 const &window_size) override;
	virtual void update(float elapsed) override;
	virtual void draw(glm::uvec2 const &drawable_size) override;

	// helper functions
	void get_transforms();
	void init(int i);
	//----- game state -----

	//input tracking:
	struct Button {
		uint8_t downs = 0;
		uint8_t pressed = 0;
	} left, right, down, up;

	//local copy of the game scene (so code can change it during gameplay):
	Scene scene;


	// Get all Transforms
	//Scene::Transform *Barrel = nullptr;
	//Scene::Transform *Enemy_Car1 = nullptr;
	//Scene::Transform *Enemy_Car2 = nullptr;
	//Scene::Transform *Enemy_Car3 = nullptr;
	//Scene::Transform *Enemy_Car4 = nullptr;
	//
	//Scene::Transform *Rain1 = nullptr;
	//Scene::Transform *Rain2 = nullptr;
	//Scene::Transform *Rain3 = nullptr;
//
	//Scene::Transform *Side_Ground = nullptr;
	//Scene::Transform *Trees1 = nullptr;
	//Scene::Transform *Trees2 = nullptr;

	Scene::Transform* obstacle[42];

	//Scene::Transform *Reference = nullptr;

		
	//Scene::Transform *WalkMesh = nullptr;

	Scene::Camera *camera = nullptr;

	// Car Logic
	float car_speed;
	float max_speed;

	// Rain Loc
	glm::vec3 Rain1_pos, Rain2_pos, Rain3_pos;
	float rain_rate;

	// Sound Value
	float rain_vol, car_vol, main1_vol;

	std::shared_ptr< Sound::PlayingSample > RainLoop;
	std::shared_ptr< Sound::PlayingSample > CarLoop;
	std::shared_ptr< Sound::PlayingSample > Main1Loop;


	//player info:
	struct Player {
		WalkPoint at;
		//transform is at player's feet and will be yawed by mouse left/right motion:
		Scene::Transform *transform = nullptr;
		//camera is at player's head and will be pitched by mouse up/down motion:
		//
	} player;
};
