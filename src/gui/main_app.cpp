#include "hello_imgui/hello_imgui_include_opengl.h"
#include "imgui.h"
#include "implot/implot.h"
#include "implot/implot_internal.h"
#include "simulation_v2.h"

#include <mpi.h>
#include <memory>
#include "immapp/runner.h"
#include "visualization.h"

enum class Command : int {
  None,
  CreateSim,
  DeleteSim,
  Step,
  GatherRho,
  GatherFF,
  GatherRhok,
  GatherPhi,
  GatherParticles,
  Shutdown,
};

void broadcast_command(Command cmd) {
  int icmd = static_cast<int>(cmd);
  MPI_Bcast(&icmd, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

class SimulationCached : public Simulation {
  void sync_particles() {
    broadcast_command(Command::GatherParticles);
    cached_particles = gather_particles();
    x_data.reserve(cached_particles.size());
    y_data.reserve(cached_particles.size());
    for (uint i = 0; i < cached_particles.size(); i++) {
      x_data[i] = cached_particles[i].p.x - params.RADIUS / 2;
      y_data[i] = cached_particles[i].p.y - params.RADIUS / 2;
    }
  }

  float density_hdr(float input) {
    input *= (0.3 * params.RADIUS * params.RADIUS / params.MASS);
    // now mostly < 1
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
    //
    //     broadcast_command(Command::GatherRhok);
    //     auto globalDensityk = gather_rhok();
    //     densitykTextureData.reserve(params.RESOLUTION * params.RESOLUTION * 3);
    // #pragma omp parallel for collapse(2)
    //     for (int i = 0; i < params.RESOLUTION; i++) {
    //       for (int j = 0; j < params.RESOLUTION; j++) {
    //         glm::vec3 color = complexColour(globalDensityk[j * params.RESOLUTION + i]);
    //
    //         densitykTextureData[(i * params.RESOLUTION + j) * 3 + 0] = color.x;
    //         densitykTextureData[(i * params.RESOLUTION + j) * 3 + 1] = color.y;
    //         densitykTextureData[(i * params.RESOLUTION + j) * 3 + 2] = color.z;
    //       }
    //     }
    //
    //     broadcast_command(Command::GatherPhi);
    //     auto globalPhiGrad = gather_phi();
    //     phiTextureData.reserve(params.RESOLUTION * params.RESOLUTION);
    // #pragma omp parallel for collapse(2)
    //     for (int i = 0; i < params.RESOLUTION; i++) {
    //       for (int j = 0; j < params.RESOLUTION; j++) {
    //         phiTextureData[i * params.RESOLUTION + j] = globalPhiGrad[j * params.RESOLUTION + i];
    //       }
    //     }
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
    //
    //   glGenTextures(1, &densitykTexture);
    //   glBindTexture(GL_TEXTURE_2D, densitykTexture);
    //   glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, params.RESOLUTION, params.RESOLUTION, 0, GL_RGB, GL_FLOAT,
    //                densitykTextureData.data());
    //   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    //   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    //
    //   glGenTextures(1, &phiTexture);
    //   glBindTexture(GL_TEXTURE_2D, phiTexture);
    //   glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, params.RESOLUTION, params.RESOLUTION, 0, GL_RED, GL_FLOAT,
    //                phiTextureData.data());
    //   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    //   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    //   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_R, GL_RED);
    //   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_G, GL_RED);
    //   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_B, GL_RED);
    //   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_A, GL_ONE);
  }

  void sync_textures() {
    load_textures();

    glBindTexture(GL_TEXTURE_2D, densityTexture);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, Nx, Ny, GL_RED, GL_FLOAT, densityTextureData.data());

    glBindTexture(GL_TEXTURE_2D, ffTexture);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, params.RESOLUTION, params.RESOLUTION, GL_RGB, GL_FLOAT,
                    ffTextureData.data());
    //
    //   glBindTexture(GL_TEXTURE_2D, densityTexture);
    //   glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, params.RESOLUTION, params.RESOLUTION, GL_RGB, GL_FLOAT,
    //                   densitykTextureData.data());
    //
    //   glBindTexture(GL_TEXTURE_2D, phiTexture);
    //   glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, params.RESOLUTION, params.RESOLUTION, GL_RED, GL_FLOAT,
    //                   phiTextureData.data());
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

  std::vector<Particle> cached_particles;
  std::vector<double> x_data = {};
  std::vector<double> y_data = {};

  SimulationCached(Params params) : Simulation(params) {
    sync_particles();
    init_textures();
  }

  void sync() {
    sync_particles();
    sync_textures();
  }
};

class App {
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

        ImPlot::PlotImage("density", simulation->densityTexture, ImPlotPoint(-hr, -hr),
                          ImPlotPoint(hr, hr));
        ImPlot::PlotImage("ff", simulation->ffTexture, ImPlotPoint(-hr, -hr), ImPlotPoint(hr, hr));

        // // Draw the particles
        // ImPlot::GetStyle().MarkerSize = particle_radius;
        // ImPlot::SetNextLineStyle(particle_color);
        // ImPlot::PlotScatter("particles", simulation->x_data.data(), simulation->y_data.data(),
        //                     std::min(simulation->cached_particles.size(), (size_t)100000));
      }

      ImPlot::PushPlotClipRect();
      // Draw the bounds
      ImDrawList* draw_list = ImPlot::GetPlotDrawList();
      ImVec2 p_min = ImPlot::PlotToPixels(ImPlotPoint(-hr, -hr));
      ImVec2 p_max = ImPlot::PlotToPixels(ImPlotPoint(hr, hr));
      ImU32 col = IM_COL32(80, 150, 80, 255);
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

    if (ImGui::CollapsingHeader("Simulation Controls")) {
      ImGui::BeginDisabled(simulation != nullptr);
      ImGui::InputFloat("Gravitational Constant", &params.GRAVITY, 0.0f, 0.0f, "%.6e",
                        ImGuiInputTextFlags_CharsScientific);
      ImGui::InputFloat("Universe Radius", &params.RADIUS, 0.0f, 0.0f, "%.6e",
                        ImGuiInputTextFlags_CharsScientific);
      ImGui::InputFloat("Universe Mass", &params.MASS, 0.0f, 0.0f, "%.6e",
                        ImGuiInputTextFlags_CharsScientific);
      ImGui::InputFloat("Integration Timestep", &params.TIMESTEP, 0.0f, 0.0f, "%e",
                        ImGuiInputTextFlags_CharsScientific);
      ImGui::InputFloat("Gravitational Softening", &params.SOFTENING, 0.0f, 0.0f, "%e",
                        ImGuiInputTextFlags_CharsScientific);
      ImGui::Checkbox("Evolve Scale Factor", &params.USE_SCALE_FACTOR);
      int resolution = params.RESOLUTION;
      ImGui::SliderInt("Density Texture Resolution", &resolution, 32, 4096);
      params.RESOLUTION = resolution;
      ImGui::SliderFloat("Particles Per Cell", &params.PARTICLES_PER_CELL, 0.01f, 100, "%f",
                         ImGuiSliderFlags_Logarithmic);
      ImGui::EndDisabled();

      if (!simulation && ImGui::Button("Create Simulation")) {
        broadcast_command(Command::CreateSim);
        MPI_Bcast(&params, sizeof(Simulation::Params), MPI_BYTE, 0, MPI_COMM_WORLD);
        simulation = std::make_unique<SimulationCached>(params);
      }

      if (simulation && ImGui::Button("Delete Simulation")) {
        simulation = nullptr;
        broadcast_command(Command::DeleteSim);
      }

      if (simulation) {
        static int num_iterations = 1;
        ImGui::SliderInt("Number of Iterations", &num_iterations, 1, 10000);

        if (ImGui::Button("Update")) {

          broadcast_command(Command::Step);
          for (int i = 0; i < num_iterations; ++i) {
            simulation->timestep();
          }
          MPI_Bcast(&num_iterations, 1, MPI_INT, 0, MPI_COMM_WORLD);
          MPI_Barrier(MPI_COMM_WORLD);
          simulation->sync();
        }
      }
    }

    if (simulation) {
      if (ImGui::CollapsingHeader("Simulation Observables")) {
        ImGui::LabelText("number of particles", "%.2e",
                         (double)simulation->cached_particles.size());
        if (params.USE_SCALE_FACTOR) {
          ImGui::LabelText("elapsed time", "%e", simulation->t);
          ImGui::LabelText("scale factor", "%e", simulation->a + simulation->t * simulation->adot);
          ImGui::LabelText("Hubble factor", "%e",
                           simulation->adot / (simulation->a + simulation->t * simulation->adot));
        }
      }
    }

    ImGui::End();
  }

 public:
  App() {}

  void init() {}

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
          for (int i = 0; i < num_iterations; i++)
            sim->timestep();
        }
        MPI_Barrier(MPI_COMM_WORLD);

        break;
      case Command::GatherRho:
        sim->gather_rho();
        break;
      case Command::GatherFF:
        sim->gather_ff();
        break;
      case Command::GatherRhok:
        // sim->gather_rhok();
        break;
      case Command::GatherPhi:
        // sim->gather_phi();
        break;
      case Command::GatherParticles:
        sim->gather_particles();
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
