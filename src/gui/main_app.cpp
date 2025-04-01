#include "glm/gtc/constants.hpp"
#include "imgui.h"
#include "imgui_app.h"
#include "implot.h"
#include "implot_internal.h"
#include "simulation.h"
#include <memory>
#include <random>

class GravitySimulator : public ImGuiApp {
  std::unique_ptr<Simulation> simulation;
  bool running = false;

public:
  GravitySimulator() {
    auto specs = AppSpecs();
    specs.name = "Newtonian Gravity Simulator";
    specs.width = 1000;
    specs.height = 1000;

    create(specs);

    std::random_device rd;
    std::mt19937 rng(rd());

    // std::uniform_real_distribution<double> angle_sample(0.0,
    //                                                     glm::tau<double>());
    // std::uniform_real_distribution<double> radius_sample(4.0, 6.0);
    //
    // double black_hole_mass = 10.0;
    //
    // std::vector<particle> particles;
    // const uint num_particles = 5000;
    //
    // particles.push_back(particle{.mass = black_hole_mass});
    //
    // for (uint i = 0; i < num_particles; i++) {
    //   double radius = radius_sample(rng);
    //   double angle = angle_sample(rng);
    //
    //   double x = cos(angle);
    //   double y = sin(angle);
    //   double speed = sqrt(GRAV * black_hole_mass / radius);
    //
    //   particles.push_back(particle{.pos = radius * vec3(x, y, 0.0),
    //                                .vel = speed * vec3(-y, x, 0.0),
    //                                .mass = 0.00000000001});
    // }
    SimulationParams params = SimulationParams{
        .grav_constant = 0.01,
        .softening = 0.001,
        .dt = 0.01,
        .radius = 10.0,
    };

    std::vector<Particle> particles;
    const double total_mass = 1.0;
    const uint num_particles = 5000;

    std::uniform_real_distribution<double> sample(-5.0, 5.0);

    for (uint i = 0; i < num_particles; i++) {
      while (true) {
        vec3 x = vec3(sample(rng), sample(rng), 0.0);

        // if (dot(x, x) < 1.0) {
        particles.push_back(
            Particle{.pos = x, .mass = total_mass / num_particles});
        break;
        // }
      }
    }

    simulation = std::make_unique<Simulation>(particles, params);
  }

  ~GravitySimulator() { destroy(); }

  void update() override {
    auto io = ImGui::GetIO();

    std::vector<double> x_data = {};
    std::vector<double> y_data = {};

    for (const auto &p : simulation->particles) {
      x_data.push_back(p.pos.x);
      y_data.push_back(p.pos.y);
    }

    {
      ImGui::SetNextWindowPos(ImVec2(0, 0));
      ImGui::SetNextWindowSize(ImGui::GetIO().DisplaySize);
      ImGui::Begin("BackPanel", nullptr,
                   ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                       ImGuiWindowFlags_NoMove |
                       ImGuiWindowFlags_NoBringToFrontOnFocus |
                       ImGuiWindowFlags_NoNav);

      ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
                  1000.0f / io.Framerate, io.Framerate);
      ImGui::Checkbox("Running", &running);

      ImGui::Columns(2, "simulation", false);
      ImGui::SetColumnWidth(0, ImGui::GetIO().DisplaySize.x * 0.5);
      if (ImPlot::BeginPlot("Viewport", ImVec2(-1.0, -1.0),
                            ImPlotFlags_Equal)) {
        ImPlot::SetupAxes("x [Mpc]", "y [Mpc]");
        ImPlot::GetStyle().MarkerSize = 1;
        ImPlot::GetStyle().FillAlpha = 0.8;
        ImPlot::PlotScatter("particles", x_data.data(), y_data.data(),
                            simulation->particles.size());

        if (ImPlot::BeginItem("bounds")) {
          float hr = simulation->params.radius / 2;
          ImVec2 p_min = ImPlot::PlotToPixels(ImPlotPoint(-hr, -hr));
          ImVec2 p_max = ImPlot::PlotToPixels(ImPlotPoint(hr, hr));
          ImDrawList *draw_list = ImPlot::GetPlotDrawList();
          draw_list->AddRect(p_min, p_max, ImPlot::GetCurrentItem()->Color,
                             0.0f, 0, 2.0f);
          ImPlot::EndItem();
        }

        ImPlot::EndPlot();
      }
      ImGui::NextColumn();

      if (ImPlot::BeginPlot("Energy", ImVec2(-1.0, 0.0))) {
        ImPlot::SetupAxes("Time [s]", "Energy", ImPlotAxisFlags_AutoFit,
                          ImPlotAxisFlags_AutoFit);
        ImPlot::SetupLegend(ImPlotLocation_NorthEast);
        ImPlot::PlotLine("total", simulation->es.data(), simulation->es.size(),
                         simulation->params.dt);
        ImPlot::PlotLine("kinetic", simulation->kes.data(),
                         simulation->kes.size(), simulation->params.dt);
        ImPlot::PlotLine("potential", simulation->pes.data(),
                         simulation->pes.size(), simulation->params.dt);
        ImPlot::EndPlot();
      }
      ImGui::EndColumns();

      ImGui::End();
    }

    if (running) {
      simulation->update();
    }
  }
};

// Main code
int main(int, char **) {
  GravitySimulator simulator;
  simulator.run();
}
