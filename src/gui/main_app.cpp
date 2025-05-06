#include "hello_imgui/hello_imgui_include_opengl.h"
#include "imgui.h"
#include "implot/implot.h"
#include "implot/implot_internal.h"
#include "simulation.h"

#include "immapp/runner.h"
#include "simulation.h"
#include <memory>
#include <mpi.h>

enum class Command : int {
  None,
  CreateSim,
  DeleteSim,
  Step,
  GatherRho,
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
    particles = gather_particles();
    x_data.reserve(particles.size());
    y_data.reserve(particles.size());
    for (uint i = 0; i < particles.size(); i++) {
      x_data[i] = particles[i].p.x;
      y_data[i] = particles[i].p.y;
    }
  }

  void load_density() {
    broadcast_command(Command::GatherRho);
    auto globalDensity = gather_phi();
    densityTextureData.reserve(params.RESOLUTION * params.RESOLUTION);
#pragma omp parallel for collapse(2)
    for (int i = 0; i < params.RESOLUTION; i++) {
      for (int j = 0; j < params.RESOLUTION; j++) {
        densityTextureData[i * params.RESOLUTION + j] =
            globalDensity[j * params.RESOLUTION + i];
      }
    }
  }

  void init_density() {
    load_density();
    glGenTextures(1, &densityTexture);
    glBindTexture(GL_TEXTURE_2D, densityTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, params.RESOLUTION,
                 params.RESOLUTION, 0, GL_RED, GL_FLOAT,
                 densityTextureData.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_R, GL_RED);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_G, GL_RED);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_B, GL_RED);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_A, GL_ONE);
  }

  void sync_density() {
    load_density();
    glBindTexture(GL_TEXTURE_2D, densityTexture);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, params.RESOLUTION,
                    params.RESOLUTION, GL_RED, GL_FLOAT,
                    densityTextureData.data());
  }

public:
  GLuint densityTexture;
  std::vector<float> densityTextureData;

  std::vector<Particle> particles;
  std::vector<double> x_data = {};
  std::vector<double> y_data = {};

  SimulationCached(Params params) : Simulation(params) {
    sync_particles();
    init_density();
  }

  void sync() {
    sync_particles();
    sync_density();
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

      if (simulation != nullptr) {
        float hr = simulation->params.RADIUS / 2.0;

        // Draw the density
        ImPlot::PlotImage("density", simulation->densityTexture,
                          ImPlotPoint(-hr, hr), ImPlotPoint(hr, -hr));

        // Draw the bounds
        if (ImPlot::BeginItem("bounds", 0)) {
          ImDrawList *draw_list = ImPlot::GetPlotDrawList();
          ImVec2 p_min = ImPlot::PlotToPixels(ImPlotPoint(-hr, -hr));
          ImVec2 p_max = ImPlot::PlotToPixels(ImPlotPoint(hr, hr));
          draw_list->AddRect(p_min, p_max, ImPlot::GetCurrentItem()->Color,
                             0.0f, 0, 2.0f);
          ImPlot::EndItem();
        }

        // if (ImPlot::BeginItem("MPI bounds", 0)) {
        //   ImDrawList *draw_list = ImPlot::GetPlotDrawList();
        //   float dx = simulation->dx;
        //   int resolution = simulation->params.RESOLUTION;
        //
        //   int n_ranks = simulation->size;
        //
        //   for (int r = 0; r < n_ranks; ++r) {
        //     float x0 =
        //         -hr + r * ceil(simulation->params.RESOLUTION / n_ranks) * dx;
        //     float x1 = -hr + (r + 1) *
        //                          ceil(simulation->params.RESOLUTION /
        //                          n_ranks) * dx;
        //
        //     ImVec2 p_min = ImPlot::PlotToPixels(ImPlotPoint(x0, -hr));
        //     ImVec2 p_max = ImPlot::PlotToPixels(ImPlotPoint(x1, hr));
        //     draw_list->AddRect(p_min, p_max, ImPlot::GetCurrentItem()->Color,
        //                        0.0f, 0, 1.5f);
        //   }
        //
        //   ImPlot::EndItem();
        // }

        // Draw the particles
        ImPlot::GetStyle().MarkerSize = particle_radius;
        ImPlot::SetNextLineStyle(particle_color);
        ImPlot::PlotScatter("particles", simulation->x_data.data(),
                            simulation->y_data.data(),
                            simulation->particles.size());
      }

      ImPlot::EndPlot();
    }
    ImGui::End();
  }

  void simulation_params_gui() {
    ImGui::Begin("Settings");
    if (ImGui::CollapsingHeader("Viewport Controls")) {
      ImGui::ColorEdit4("Particle Color", (float *)&particle_color);
      ImGui::SliderFloat("Particle Radius", &particle_radius, 0.01, 10.0);
    }

    if (ImGui::CollapsingHeader("Simulation Controls")) {
      ImGui::BeginDisabled(simulation != nullptr);
      ImGui::SliderFloat("Gravitational Constant", &params.GRAVITY, 1e-4, 1e4,
                         "%e", ImGuiSliderFlags_Logarithmic);
      ImGui::SliderFloat("Universe Radius", &params.RADIUS, 1, 1000, "%e");
      ImGui::SliderFloat("Universe Mass", &params.MASS, 1e-4, 1e4, "%e",
                         ImGuiSliderFlags_Logarithmic);
      ImGui::SliderFloat("Integration Timestep", &params.TIMESTEP, 1e-4, 1e4,
                         "%e", ImGuiSliderFlags_Logarithmic);
      ImGui::SliderFloat("Gravitational Softening", &params.SOFTENING, 1e-4,
                         1e4, "%e", ImGuiSliderFlags_Logarithmic);
      int resolution = params.RESOLUTION;
      ImGui::SliderInt("Density Texture Resolution", &resolution, 32, 4096);
      params.RESOLUTION = resolution;
      ImGui::SliderFloat("Particles Per Cell", &params.PARTICLES_PER_CELL,
                         0.01f, 100, "%f", ImGuiSliderFlags_Logarithmic);
      ImGui::EndDisabled();

      if (!simulation && ImGui::Button("Create Simulation")) {
        broadcast_command(Command::CreateSim);
        MPI_Bcast(&params, sizeof(Simulation::Params), MPI_BYTE, 0,
                  MPI_COMM_WORLD);
        simulation = std::make_unique<SimulationCached>(params);
      }

      if (simulation && ImGui::Button("Delete Simulation")) {
        simulation = nullptr;
        broadcast_command(Command::DeleteSim);
      }

      if (simulation && ImGui::Button("Update")) {
        broadcast_command(Command::Step);
        MPI_Barrier(MPI_COMM_WORLD);
        simulation->sync();
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
      MPI_Bcast(&params, sizeof(Simulation::Params), MPI_BYTE, 0,
                MPI_COMM_WORLD);
      sim = std::make_unique<Simulation>(params);
      break;
    case Command::DeleteSim:
      sim = nullptr;
      break;
    case Command::Step:
      if (sim) {
        sim->timestep();
        MPI_Barrier(MPI_COMM_WORLD);
      }

      break;
    case Command::GatherRho:
      sim->gather_phi();
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

int main(int argc, char **argv) {
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
    params.callbacks.PostInit = [&app]() { app->init(); };
    params.callbacks.ShowGui = [&app]() { app->gui(); };
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
