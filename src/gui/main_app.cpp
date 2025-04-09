#include "imgui.h"
#include "imgui_internal.h"
#include "implot/implot.h"
#include "implot/implot_internal.h"
#include "simulation.h"

#include "immapp/runner.h"
#include <iostream>
#include <random>

#include "hello_imgui/hello_imgui_include_opengl.h"
#include "imgui.h"
#include "simulation.h"

class GravitySimulator {
  std::unique_ptr<Simulation> sim;

  std::vector<double> x_data = {};
  std::vector<double> y_data = {};
  bool updated = false;

  std::vector<float> density = {};
  GLuint densityTexture;

  ImVec4 particle_color = ImVec4(1.0, 1.0, 1.0, 1.0);
  float particle_radius = 0.5;
  ImVec4 bounds_color = ImVec4(1.0, 1.0, 1.0, 1.0);

public:
  GravitySimulator() {
    std::random_device rd;
    std::mt19937 rng(rd());

    std::vector<Particle> particles;
    const double total_mass = 10.0;
    const uint num_particles = 100000;

    std::uniform_real_distribution<double> sample(-5.0, 5.0);

    for (uint i = 0; i < num_particles; i++) {
      vec2 x = vec2(sample(rng), sample(rng));

      particles.push_back(
          Particle{.pos = x, .mass = total_mass / num_particles});
    }

    sim = std::make_unique<Simulation>(particles);
    density.resize(sim->resolution * sim->resolution, 0.0);
  }

  bool running = false;

  void init() {
    glGenTextures(1, &densityTexture);
    glBindTexture(GL_TEXTURE_2D, densityTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, sim->resolution, sim->resolution, 0,
                 GL_RED, GL_FLOAT, density.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_R, GL_RED);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_G, GL_RED);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_B, GL_RED);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_A, GL_ONE);
  }

  void update() {
    auto io = ImGui::GetIO();

    ImGui::Columns(2, "viewport", false);
    ImGui::SetColumnWidth(0, ImGui::GetIO().DisplaySize.x * 0.7);

    if (!updated) {
      x_data.reserve(sim->N);
      y_data.reserve(sim->N);
      for (uint i = 0; i < sim->N; i++) {
        x_data[i] = sim->particles[i].pos.x;
        y_data[i] = sim->particles[i].pos.y;
      }

      for (uint i = 0; i < sim->resolution * sim->resolution; i++) {
        density[i] = (sim->density[i] / sim->total_mass) *
                     (sim->resolution * sim->resolution) * 0.2;
      }

      glBindTexture(GL_TEXTURE_2D, densityTexture);
      glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, sim->resolution, sim->resolution,
                      GL_RED, GL_FLOAT, density.data());
      updated = true;
    }

    if (ImPlot::BeginPlot("Viewport", ImVec2(-1.0, -1.0), ImPlotFlags_Equal)) {
      ImPlot::SetupAxes("x [Mpc]", "y [Mpc]");
      ImPlot::GetStyle().MarkerSize = particle_radius;

      float hr = sim->radius / 2;
      ImPlot::PlotImage("density", densityTexture, ImPlotPoint(-hr, hr),
                        ImPlotPoint(hr, -hr));
      ImPlot::SetNextLineStyle(particle_color);
      ImPlot::PlotScatter("particles", x_data.data(), y_data.data(), sim->N);

      ImPlot::SetNextFillStyle(bounds_color);
      if (ImPlot::BeginItem("bounds", 0)) {
        ImDrawList *draw_list = ImPlot::GetPlotDrawList();
        ImVec2 p_min = ImPlot::PlotToPixels(ImPlotPoint(-hr, -hr));
        ImVec2 p_max = ImPlot::PlotToPixels(ImPlotPoint(hr, hr));
        draw_list->AddRect(p_min, p_max, ImPlot::GetCurrentItem()->Color, 0.0f,
                           0, 2.0f);
        ImPlot::EndItem();
      }

      ImPlot::EndPlot();
      ImGui::NextColumn();

      ImGui::Separator();
      if (ImGui::CollapsingHeader("Viewport Settings")) {
        ImGui::ColorEdit4("Particle Color", (float *)&particle_color);
        ImGui::SliderFloat("Particle Radius", &particle_radius, 0.01, 10.0);
        ImGui::ColorEdit4("Bounds Color", (float *)&bounds_color);
      }

      ImGui::Separator();
      ImGui::Text("Number of Particles: %d", (int)sim->N);
      ImGui::Text("softening: %f", sim->softening);
      ImGui::Text("dt: %f", sim->dt);
      ImGui::Text("G: %f", sim->grav_constant);
      ImGui::Text("R: %f", sim->radius);
      ImGui::Text("total mass: %f", sim->total_mass);
      ImGui::Text("resolution: %d", sim->resolution);

      ImGui::Separator();

      ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
                  1000.0f / io.Framerate, io.Framerate);
      ImGui::Checkbox("Running", &running);

      ImGui::EndColumns();
    }

    if (running) {
      sim->update();
      updated = false;
    }
  }
};

void showWindow(std::string name) {
  auto windowPtr =
      HelloImGui::GetRunnerParams()->dockingParams.dockableWindowOfName(name);
  if (windowPtr) {
    windowPtr->isVisible = true;
  }
}

int main(int, char **) {
  GravitySimulator simulator;

  HelloImGui::RunnerParams runnerParams;
  runnerParams.appWindowParams.windowTitle =
      "2D Cosmological Structure Formation Simulator";
  runnerParams.imGuiWindowParams.menuAppTitle = "Main";
  runnerParams.appWindowParams.windowGeometry.size = {1920, 1080};
  runnerParams.imGuiWindowParams.showMenuBar = true;
  runnerParams.imGuiWindowParams.showStatusBar = false;
  runnerParams.fpsIdling.enableIdling = false;
  runnerParams.callbacks.PostInit = [&simulator]() { simulator.init(); };
  runnerParams.callbacks.ShowGui = [&simulator]() { simulator.update(); };

  ImmApp::AddOnsParams addOnsParams;
  addOnsParams.withImplot = true;

  ImmApp::Run(runnerParams, addOnsParams);

  return 0;
}
