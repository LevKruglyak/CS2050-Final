#include "hello_imgui/hello_imgui_include_opengl.h"
#include "imgui.h"
#include "imgui_stacklayout.h"
#include "implot/implot.h"
#include "implot/implot_internal.h"
#include "simulation.h"

#include <mpi.h>
#include <algorithm>
#include <memory>
#include "immapp/runner.h"

const int TAG_PROGRESS = 3;
const int TAG_DONE = 2;

enum class Command : int {
  None,
  CreateSim,
  DeleteSim,
  Step,
  GatherRho,
  GatherFF,
  Shutdown,
};

void broadcast_command(Command cmd) {
  int icmd = static_cast<int>(cmd);
  MPI_Bcast(&icmd, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

inline glm::vec2 cxlog(const glm::vec2& z) {
  const float r = glm::length(z);
  float theta = std::atan2(z.y, z.x);
  return glm::vec2(std::log(r), theta);
}
inline glm::vec3 hsv2rgb(const glm::vec3& hsv) {
  const glm::vec4 K(1.0f, 2.0f / 3.0f, 1.0f / 3.0f, 3.0f);
  const glm::vec3 p =
      glm::abs(glm::fract(glm::vec3(hsv.x) + glm::vec3(K.x, K.y, K.z)) * 6.0f - glm::vec3(K.w));
  return hsv.z * glm::mix(glm::vec3(K.x), glm::clamp(p - glm::vec3(K.x), 0.0f, 1.0f), hsv.y);
}

inline float hdrTone(float r, float exposure = 1.0f) {
  return 1.0f - std::exp(-exposure * r);
}

inline glm::vec3 complexColour(const glm::vec2& z) {
  const float hue = cxlog(z).y / (2.0f * M_PI);
  const float val = hdrTone(glm::length(z), 1.0);
  const glm::vec3 hsv(hue, 1.0f, val);
  return hsv2rgb(hsv);
}

class SimulationCached : public Simulation {
  float density_hdr(float input) {
    input *= (0.3 * params.RADIUS * params.RADIUS / params.MASS);
    return 1.0 - exp(-input * 3.0);
  }

  void load_textures() {
    broadcast_command(Command::GatherRho);
    auto globalDensity = gather_rho();
    densityTextureData.reserve(params.RESOLUTION * params.RESOLUTION);
#pragma omp parallel for collapse(2)
    for (int i = 0; i < params.RESOLUTION; i++) {
      for (int j = 0; j < params.RESOLUTION; j++) {
        densityTextureData[i * params.RESOLUTION + j] =
            density_hdr(globalDensity[i * params.RESOLUTION + j]);
      }
    }

    broadcast_command(Command::GatherFF);
    auto globalFF = gather_ff();
    ffTextureData.reserve(params.RESOLUTION * params.RESOLUTION * 3);
#pragma omp parallel for collapse(2)
    for (int i = 0; i < params.RESOLUTION; i++) {
      for (int j = 0; j < params.RESOLUTION; j++) {
        glm::vec3 color = complexColour(globalFF[j * params.RESOLUTION + i]);

        ffTextureData[(i * params.RESOLUTION + j) * 3 + 0] = color.x;
        ffTextureData[(i * params.RESOLUTION + j) * 3 + 1] = color.y;
        ffTextureData[(i * params.RESOLUTION + j) * 3 + 2] = color.z;
      }
    }
  }

  void init_textures() {
    load_textures();
    glGenTextures(1, &densityTexture);
    glBindTexture(GL_TEXTURE_2D, densityTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, params.RESOLUTION, params.RESOLUTION, 0, GL_RED,
                 GL_FLOAT, densityTextureData.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_R, GL_RED);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_G, GL_RED);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_B, GL_RED);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_A, GL_ONE);

    glGenTextures(1, &ffTexture);
    glBindTexture(GL_TEXTURE_2D, ffTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, params.RESOLUTION, params.RESOLUTION, 0, GL_RGB,
                 GL_FLOAT, ffTextureData.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  }

  void sync_textures() {
    load_textures();

    glBindTexture(GL_TEXTURE_2D, densityTexture);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, Nx, Ny, GL_RED, GL_FLOAT, densityTextureData.data());

    glBindTexture(GL_TEXTURE_2D, ffTexture);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, params.RESOLUTION, params.RESOLUTION, GL_RGB, GL_FLOAT,
                    ffTextureData.data());
  }

 public:
  GLuint densityTexture;
  std::vector<float> densityTextureData;

  GLuint densitykTexture;
  std::vector<float> densitykTextureData;

  GLuint phiTexture;
  std::vector<float> phiTextureData;

  GLuint ffTexture;
  std::vector<float> ffTextureData;

  SimulationCached(Params params) : Simulation(params) { init_textures(); }
  void sync() { sync_textures(); }
};

class App {
  bool busy = false;
  float particle_radius = 0.5;
  ImVec4 particle_color = ImVec4(1.0, 1.0, 1.0, 1.0);

  std::unique_ptr<SimulationCached> simulation = nullptr;
  Simulation::Params params;

  void viewport_gui() {
    ImGui::Begin("Viewport");
    if (ImPlot::BeginPlot("Viewport", ImVec2(-1.0, -1.0),
                          ImPlotFlags_Equal | ImPlotFlags_NoTitle)) {
      ImPlot::SetupAxes("", "");
      float hr = params.RADIUS / 2.0;

      if (simulation != nullptr) {
        ImPlot::PlotImage("ff", simulation->ffTexture, ImPlotPoint(-hr, -hr), ImPlotPoint(hr, hr));
        ImPlot::PlotImage("density", simulation->densityTexture, ImPlotPoint(-hr, -hr),
                          ImPlotPoint(hr, hr));
      }

      ImPlot::PushPlotClipRect();
      // Draw the bounds
      ImDrawList* draw_list = ImPlot::GetPlotDrawList();
      ImVec2 p_min = ImPlot::PlotToPixels(ImPlotPoint(-hr, -hr));
      ImVec2 p_max = ImPlot::PlotToPixels(ImPlotPoint(hr, hr));
      ImU32 col = ImGui::ColorConvertFloat4ToU32(ImGui::GetStyleColorVec4(ImGuiCol_ButtonActive));

      draw_list->AddRect(p_min, p_max, col, 0.0f, 0, 2.0f);
      ImPlot::PopPlotClipRect();

      ImPlot::EndPlot();
    }
    ImGui::End();
  }

  void simulation_params_gui() {
    ImGui::Begin("Settings");
    if (ImGui::CollapsingHeader("Viewport Controls")) {
      ImGui::ColorEdit4("Particle Color", (float*)&particle_color);
      ImGui::SliderFloat("Particle Radius", &particle_radius, 0.01, 10.0);
    }

    ImGui::BeginDisabled(simulation != nullptr);
    ImGui::SeparatorText("Phyiscal Paramters");
    ImGui::InputFloat("Gravitational Constant", &params.GRAVITY, 0.0f, 0.0f, "%.3f",
                      ImGuiInputTextFlags_CharsScientific);
    ImGui::InputFloat("Universe Radius", &params.RADIUS, 0.0f, 0.0f, "%.3f",
                      ImGuiInputTextFlags_CharsScientific);
    ImGui::InputFloat("Universe Mass", &params.MASS, 0.0f, 0.0f, "%.3f",
                      ImGuiInputTextFlags_CharsScientific);
    ImGui::Checkbox("Evolve Scale Factor", &params.USE_SCALE_FACTOR);

    ImGui::SeparatorText("Numerical Accuracy");
    ImGui::InputFloat("Integration Timestep", &params.TIMESTEP, 0.0f, 0.0f, "%.3f",
                      ImGuiInputTextFlags_CharsScientific);
    ImGui::InputFloat("Gravitational Softening", &params.SOFTENING, 0.0f, 0.0f, "%.3f",
                      ImGuiInputTextFlags_CharsScientific);
    int resolution = params.RESOLUTION;
    ImGui::InputInt("Density Texture Resolution", &resolution, 0, 0);
    params.RESOLUTION = resolution;
    ImGui::InputFloat("Particles Per Cell", &params.PARTICLES_PER_CELL, 0.0f, 0.0f, "%.3f");

    ImGui::SeparatorText("Initial Conditions");
    ImGui::InputFloat("Perlin Noise Scale", &params.PERLIN_NOISE_SCALE, 0.0f, 0.0f, "%.3f",
                      ImGuiInputTextFlags_CharsScientific);
    ImGui::InputFloat("Perlin Noise Perturbation", &params.PERLIN_NOISE_PERTURBATION, 0.0f, 0.0f,
                      "%.3f", ImGuiInputTextFlags_CharsScientific);
    ImGui::InputInt("Perlin Noise Octaves", &params.PERLIN_NOISE_OCTAVES, 0, 0);
    params.PERLIN_NOISE_OCTAVES = std::clamp(params.PERLIN_NOISE_OCTAVES, 1, 20);
    ImGui::EndDisabled();

    ImGui::SeparatorText("Actions");

    ImGui::BeginHorizontal("actions");

    ImGui::BeginDisabled(simulation != nullptr);
    if (ImGui::Button("Create Simulation")) {
      broadcast_command(Command::CreateSim);
      MPI_Bcast(&params, sizeof(Simulation::Params), MPI_BYTE, 0, MPI_COMM_WORLD);
      simulation = std::make_unique<SimulationCached>(params);
    }
    ImGui::EndDisabled();

    ImGui::BeginDisabled(simulation == nullptr);
    if (ImGui::Button("Delete Simulation")) {
      simulation = nullptr;
      broadcast_command(Command::DeleteSim);
    }

    static int num_iterations = 10;
    ImGui::BeginDisabled(busy);
    if (ImGui::Button("Advance")) {
      broadcast_command(Command::Step);
      for (int i = 0; i < num_iterations; ++i) {
        simulation->timestep();
      }
      MPI_Bcast(&num_iterations, 1, MPI_INT, 0, MPI_COMM_WORLD);
      busy = true;
    }

    static float done = 0.0;
    MPI_Status st;
    int flag;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &st);

    if (st.MPI_TAG == TAG_PROGRESS) {
      MPI_Recv(&done, 1, MPI_FLOAT, st.MPI_SOURCE, TAG_PROGRESS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else if (st.MPI_TAG == TAG_DONE) {
      char dummy;
      MPI_Recv(&dummy, 1, MPI_BYTE, st.MPI_SOURCE, TAG_DONE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      busy = false;

      simulation->sync();
    }

    ImGui::PushItemWidth(80);
    ImGui::InputInt("Timesteps", &num_iterations, 0, 0);
    ImGui::PopItemWidth();
    ImGui::EndDisabled();
    ImGui::EndHorizontal();
    ImGui::EndDisabled();

    if (busy) {
      ImGui::ProgressBar(done, ImVec2(385, 0));
    }

    ImGui::SeparatorText("Information");
    ImGui::LabelText("Number of Worker Nodes", "%d", wsize - 1);
    if (simulation) {
      ImGui::LabelText("Number of Particles", "%.2e", (double)simulation->N);
      if (params.USE_SCALE_FACTOR) {
        ImGui::LabelText("Elapsed Time", "%e", simulation->t);
        ImGui::LabelText("Scale Factor", "%e", simulation->a + simulation->t * simulation->adot);
        ImGui::LabelText("Hubble Factor", "%e",
                         simulation->adot / (simulation->a + simulation->t * simulation->adot));
      }
    }

    ImGui::End();
  }

 public:
  int wsize;
  App() {}

  void init() { MPI_Comm_size(MPI_COMM_WORLD, &wsize); }

  void gui() {
    ImGui::DockSpaceOverViewport(ImGui::GetMainViewport()->ID);

    viewport_gui();
    simulation_params_gui();
  }
};

void worker_loop() {
  std::unique_ptr<Simulation> sim = nullptr;
  Simulation::Params params;

  while (true) {
    int icmd;
    MPI_Bcast(&icmd, 1, MPI_INT, 0, MPI_COMM_WORLD);
    Command cmd = static_cast<Command>(icmd);

    switch (cmd) {
      case Command::CreateSim:
        MPI_Bcast(&params, sizeof(Simulation::Params), MPI_BYTE, 0, MPI_COMM_WORLD);
        sim = std::make_unique<Simulation>(params);
        break;
      case Command::DeleteSim:
        sim = nullptr;
        break;
      case Command::Step:
        int num_iterations;
        MPI_Bcast(&num_iterations, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (sim) {
          for (int i = 0; i < num_iterations; i++) {
            sim->timestep();

            if (sim->rank == 0) {
              float done = float(i) / num_iterations;
              MPI_Request req;
              MPI_Isend(&done, 1, MPI_FLOAT, 0, TAG_PROGRESS, MPI_COMM_WORLD, &req);
            }
          }

          {
            char dummy = 0;
            MPI_Send(&dummy, 1, MPI_BYTE, 0, TAG_DONE, MPI_COMM_WORLD);
          }
        }
        break;
      case Command::GatherRho:
        sim->gather_rho();
        break;
      case Command::GatherFF:
        sim->gather_ff();
        break;
      case Command::Shutdown:
        return;
      default:
        break;
    }
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    std::unique_ptr<App> app = std::make_unique<App>();
    HelloImGui::RunnerParams params;
    params.appWindowParams.windowTitle = "CS205 Final Project";
    params.imGuiWindowParams.menuAppTitle = "Main";
    params.appWindowParams.windowGeometry.size = {1920, 1080};
    params.imGuiWindowParams.showMenuBar = true;
    params.imGuiWindowParams.showStatusBar = false;
    params.imGuiWindowParams.defaultImGuiWindowType =
        HelloImGui::DefaultImGuiWindowType::ProvideFullScreenDockSpace;
    params.fpsIdling.enableIdling = false;
    params.callbacks.PostInit = [&app]() {
      app->init();
    };
    params.callbacks.ShowGui = [&app]() {
      app->gui();
    };
    params.callbacks.BeforeExit = [&app]() {
      broadcast_command(Command::Shutdown);
      app = nullptr;
    };

    ImmApp::AddOnsParams addons;
    addons.withImplot = true;
    ImmApp::Run(params, addons);
  } else {
    worker_loop();
  }

  MPI_Finalize();
  return 0;
}
