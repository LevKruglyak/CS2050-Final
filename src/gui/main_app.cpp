#include "imgui.h"
#include "imgui_app.h"
#include "implot.h"

class GravitySimulator : public ImGuiApp {
public:
  GravitySimulator() {
    auto specs = AppSpecs();
    specs.name = "Newtonian Gravity Simulator";
    specs.width = 1000;
    specs.height = 1000;

    create(specs);
  }
  ~GravitySimulator() { destroy(); }

  void update() override {
    auto io = ImGui::GetIO();

    float x_data[1000] = {1.0, 0.0};
    float y_data[1000] = {1.0, 0.0};

    {
      ImGui::Begin("Window");
      ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
                  1000.0f / io.Framerate, io.Framerate);

      if (ImPlot::BeginPlot("My Plot", ImVec2(-1.0, -1.0))) {
        ImPlot::PlotScatter("Scatter Plot", x_data, y_data, 1000);
        ImPlot::EndPlot();
      }

      ImGui::End();
    }
  }
};

// Main code
int main(int, char **) {
  GravitySimulator simulator;
  simulator.run();
}
