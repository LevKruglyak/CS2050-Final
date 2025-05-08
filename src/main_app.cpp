#include "hello_imgui/hello_imgui_include_opengl.h"
#include "imgui.h"
#include "imgui_stacklayout.h"
#include "implot/implot.h"
#include "implot/implot_internal.h"
#include "simulation.h"

#include <fftw3.h>
#include <mpi.h>
#include <algorithm>
#include <memory>
#include "immapp/runner.h"

namespace ImGui {
struct LkxProfilerTask {
  double duration;
  const char* name;
  uint32_t color;

  LkxProfilerTask(double d, const char* n, uint32_t c) : duration(d), name(n), color(c) {}
};

class LkxProfiler {
 public:
  int frameWidth = 3;
  int frameSpacing = 1;
  bool useColoredLegendText = false;
  float maxFrameTime = 1.f / 30.f;

  explicit LkxProfiler(size_t framesCount = 300) : _frames(framesCount), _head(0) {}

  void LoadFrameData(const LkxProfilerTask* tasks, size_t count, double total) {
    Frame& dst = _frames[_head];
    dst.tasks.assign(tasks, tasks + count);
    dst.total = total;
    _head = (_head + 1) % _frames.size();
  }

  void Draw(int legendWidth, int height, int frameIndexOffset = 0) {
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const ImVec2 org = ImGui::GetCursorScreenPos();
    const float fullW = ImGui::CalcItemWidth();
    const float graphW = fullW - legendWidth;
    const ImVec2 graphTL(org.x, org.y);
    const ImVec2 graphBR(org.x + graphW, org.y + height);
    const ImVec2 legendTL(org.x + graphW, org.y);
    const ImVec2 legendBR(org.x + fullW, org.y + height);

    float maxTotal = 0.0f;
    for (auto& frame : _frames) {
      float sum = 0.0f;
      for (auto& t : frame.tasks)
        sum += float(t.duration);
      maxTotal = std::max(maxTotal, sum);
    }
    if (maxTotal <= 0.0f)
      maxTotal = maxFrameTime;
    maxTotal *= 1.2;

    dl->AddRectFilled(graphTL, graphBR, IM_COL32(40, 40, 40, 255));

    ImGui::PushClipRect(graphTL, graphBR, true);
    for (size_t f = 0; f < _frames.size(); ++f) {
      size_t idx = (_head + _frames.size() - 1 - frameIndexOffset - f) % _frames.size();
      float x1 = graphBR.x - (f + 1) * (frameWidth + frameSpacing);
      float x0 = x1 - frameWidth;
      if (x1 < graphTL.x)
        break;

      float yCurr = graphBR.y;
      for (auto& t : _frames[idx].tasks) {
        float h = std::max((float(t.duration) / maxTotal) * height, 1.0f);
        float y0 = yCurr;
        float y1 = yCurr - h;
        dl->AddRectFilled(ImVec2(x0, y1), ImVec2(x1, y0), t.color);
        yCurr = y1 - frameSpacing;
      }
    }
    ImGui::PopClipRect();
    dl->AddRect(graphTL, graphBR, IM_COL32(180, 180, 180, 160));

    const auto& latest = _frames[(_head + _frames.size() - 1) % _frames.size()];
    float yBaseL = legendBR.y;
    float yBase = legendBR.y;
    const float textH = ImGui::GetTextLineHeight() + 2.f;

    if (latest.tasks.size()) {
      ImVec2 R0(legendTL.x + 3.f + 5.f + 30.f, yBase - 3.f);
      ImVec2 R1(R0.x + 10.f, R0.y - 10.f);
      char buf[64];
      float ms = float(latest.total) * 1000.f;
      snprintf(buf, 64, "[%.2fms] Total", ms);
      dl->AddText(ImVec2(R1.x + 5.f, R1.y - 3.f), ImGui::GetColorU32(ImGuiCol_Text), buf);
      yBase -= (textH);
    }

    for (size_t i = 0; i < latest.tasks.size(); ++i) {
      const auto& t = latest.tasks[i];
      if (yBase - textH < legendTL.y)
        break;

      float endH = std::max((float(t.duration) / maxTotal) * height, 1.0f);
      ImVec2 L0(legendTL.x + 3.f, yBaseL);
      ImVec2 L1(L0.x + 5.f, yBaseL - endH);
      dl->AddRectFilled(L0, L1, t.color);
      yBaseL -= endH + frameSpacing;

      ImVec2 R0(legendTL.x + 3.f + 5.f + 30.f, yBase - 3.f);
      ImVec2 R1(R0.x + 10.f, R0.y - 10.f);
      dl->AddRectFilled(R0, R1, t.color);
      yBase -= (textH);

      ImVec2 pts[4] = {{L1.x, L0.y}, {L1.x, L1.y}, {R0.x, R1.y}, {R0.x, R0.y}};
      dl->AddConvexPolyFilled(pts, 4, t.color);

      char buf[64];
      float ms = float(t.duration) * 1000.f;
      snprintf(buf, 64, "[%.2fms] %s", ms, t.name);
      dl->AddText(ImVec2(R1.x + 5.f, R1.y - 3.f), t.color, buf);
    }

    ImGui::Dummy(ImVec2(fullW, height));
  }

  bool empty() { return _frames.empty(); }

 private:
  struct Frame {
    std::vector<LkxProfilerTask> tasks;
    double total;
  };
  std::vector<Frame> _frames;
  size_t _head = 0;
};
}  // namespace ImGui

#define RGBA_LE(col)                                                   \
  (((col & 0xff000000) >> (3 * 8)) + ((col & 0x00ff0000) >> (1 * 8)) + \
   ((col & 0x0000ff00) << (1 * 8)) + ((col & 0x000000ff) << (3 * 8)))

const int TAG_PROGRESS = 3;
const int TAG_DONE = 2;

enum class Command : int {
  None,
  CreateSim,
  DeleteSim,
  Step,
  GatherRho,
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

constexpr glm::vec3 catmullRom(const glm::vec3& p0, const glm::vec3& p1, const glm::vec3& p2,
                               const glm::vec3& p3, float t) {
  const float t2 = t * t;
  const float t3 = t2 * t;
  return 0.5f * ((2.0f * p1) + (-p0 + p2) * t + (2.0f * p0 - 5.0f * p1 + 4.0f * p2 - p3) * t2 +
                 (-p0 + 3.0f * p1 - 3.0f * p2 + p3) * t3);
}

inline constexpr std::array<glm::vec3, 256> makeCosmicLUT() {
  struct Stop {
    float pos;
    glm::vec3 c;
  };

  constexpr Stop stops[] = {{0.00f, {0 / 255.f, 0 / 255.f, 0 / 255.f}},
                            {0.20f, {45 / 255.f, 20 / 255.f, 71 / 255.f}},
                            {0.40f, {92 / 255.f, 43 / 255.f, 111 / 255.f}},
                            {0.60f, {125 / 255.f, 55 / 255.f, 116 / 255.f}},
                            {0.8f, {245 / 255.f, 179 / 255.f, 50 / 255.f}},
                            {0.95, {248 / 255.f, 239 / 255.f, 159 / 255.f}},
                            {1.00f, {255 / 255.f, 255 / 255.f, 255 / 255.f}}};
  constexpr std::size_t N = std::size(stops);

  std::array<glm::vec3, 256> lut{};

  for (std::size_t i = 0; i < lut.size(); ++i) {
    const float t = static_cast<float>(i) / 255.0f;

    std::size_t s = 0;
    while (s + 1 < N && t > stops[s + 1].pos)
      ++s;
    const Stop& P0 = stops[(s == 0) ? s : s - 1];
    const Stop& P1 = stops[s];
    const Stop& P2 = stops[s + 1];
    const Stop& P3 = stops[(s + 2 < N) ? s + 2 : s + 1];

    const float segT = std::clamp((t - P1.pos) / (P2.pos - P1.pos), 0.0f, 1.0f);

    lut[i] = catmullRom(P0.c, P1.c, P2.c, P3.c, segT);
  }

  return lut;
}

inline glm::vec3 mapHDRtoColor(float v, float maxV, const std::array<glm::vec3, 256>& lut) {
  if (v <= 0.0f)
    return lut[0];
  float t = std::log(1.0f + v) / std::log(1.0f + maxV);
  t = std::clamp(t, 0.0f, 1.0f);
  int idx = int(t * 255.0f + 0.5f);
  return lut[idx];
}

class SimulationCached : public Simulation {
  glm::fvec3 density_hdr(float input) {
    input *= params.RADIUS * params.RADIUS / params.MASS;
    static constexpr std::array<glm::vec3, 256> lut = makeCosmicLUT();
    return mapHDRtoColor(1.0 - exp(-0.4 * input), 1.0, lut);
  }

  void load_textures() {
    broadcast_command(Command::GatherRho);
    auto globalDensity = gather_rho();
    densityTextureData.reserve(params.RESOLUTION * params.RESOLUTION * 3);
#pragma omp parallel for collapse(2)
    for (int i = 0; i < params.RESOLUTION; i++) {
      for (int j = 0; j < params.RESOLUTION; j++) {
        glm::fvec3 col = density_hdr(globalDensity[i * params.RESOLUTION + j]);
        densityTextureData[(i * params.RESOLUTION + j) * 3 + 0] = col.x;
        densityTextureData[(i * params.RESOLUTION + j) * 3 + 1] = col.y;
        densityTextureData[(i * params.RESOLUTION + j) * 3 + 2] = col.z;
      }
    }
  }

  void init_textures() {
    load_textures();

    glGenTextures(1, &densityTexture);
    glBindTexture(GL_TEXTURE_2D, densityTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, params.RESOLUTION, params.RESOLUTION, 0, GL_RGB,
                 GL_FLOAT, densityTextureData.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  }

  void sync_textures() {
    load_textures();

    glBindTexture(GL_TEXTURE_2D, densityTexture);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, Nx, Ny, GL_RGB, GL_FLOAT, densityTextureData.data());
  }

 public:
  GLuint densityTexture;
  std::vector<float> densityTextureData;

  ImGui::LkxProfiler profiler{300};

  SimulationCached(Params params) : Simulation(params) { init_textures(); }
  void sync() { sync_textures(); }
};

class App {
  bool busy = false;
  std::unique_ptr<SimulationCached> simulation = nullptr;
  Simulation::Params params;

  void viewport_gui() {
    ImGui::Begin("Viewport");
    if (ImPlot::BeginPlot("Viewport", ImVec2(-1.0, -1.0),
                          ImPlotFlags_Equal | ImPlotFlags_NoTitle | ImPlotFlags_NoLegend)) {
      ImPlot::SetupAxes("", "");
      float r = params.RADIUS;
      if (simulation != nullptr) {
        ImPlot::PlotImage("density", simulation->densityTexture, ImPlotPoint(0, 0),
                          ImPlotPoint(r, r));
      }

      ImPlot::PushPlotClipRect();
      // Draw the bounds
      ImDrawList* draw_list = ImPlot::GetPlotDrawList();
      ImVec2 p_min = ImPlot::PlotToPixels(ImPlotPoint(0, 0));
      ImVec2 p_max = ImPlot::PlotToPixels(ImPlotPoint(r, r));
      ImU32 col = ImGui::ColorConvertFloat4ToU32(ImGui::GetStyleColorVec4(ImGuiCol_ButtonActive));

      draw_list->AddRect(p_min, p_max, col, 0.0f, 0, 2.0f);
      ImPlot::PopPlotClipRect();

      ImPlot::EndPlot();
    }
    ImGui::End();
  }

  void simulation_params_gui() {
    ImGui::Begin("Settings");
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
      MPI_Bcast(&num_iterations, 1, MPI_INT, 0, MPI_COMM_WORLD);
      busy = true;
    }

    static float progress = 0.0;
    MPI_Status st;
    int flag;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &st);

    if (st.MPI_TAG == TAG_PROGRESS) {
      MPI_Recv(&progress, 1, MPI_FLOAT, st.MPI_SOURCE, TAG_PROGRESS, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      {
        auto profile = simulation->timestep();

        const static uint32_t color1 = RGBA_LE(0x1abc9cffu);
        const static uint32_t color2 = RGBA_LE(0x16a085ffu);
        const static uint32_t color3 = RGBA_LE(0x2ecc71ffu);
        const static uint32_t color4 = RGBA_LE(0x27ae60ffu);
        const static uint32_t color5 = RGBA_LE(0x3498dbffu);
        const static uint32_t color6 = RGBA_LE(0x2980b9ffu);
        const static uint32_t color7 = RGBA_LE(0x9b59b6ffu);
        const static uint32_t color8 = RGBA_LE(0x8e44adffu);
        const static uint32_t color9 = RGBA_LE(0xf1c40fffu);
        const static uint32_t color10 = RGBA_LE(0xf39c12ffu);
        const static uint32_t color11 = RGBA_LE(0xe67e22ffu);
        const static uint32_t color12 = RGBA_LE(0xd35400ffu);
        const static uint32_t color13 = RGBA_LE(0xe74c3cffu);
        const static uint32_t color14 = RGBA_LE(0xc0392bffu);

        std::vector<ImGui::LkxProfilerTask> tasks;
        // durations only:
        tasks.emplace_back(profile.update_positions, "Update Positions", color1);
        tasks.emplace_back(profile.mass_local_accum, "Accumulate Local Masses", color2);
        tasks.emplace_back(profile.mass_halo_exchange, "Exchange Mass Halos", color3);
        tasks.emplace_back(profile.fft_forward, "Forward FFT", color4);
        tasks.emplace_back(profile.spectral_solve, "Spectral Gradient", color5);
        tasks.emplace_back(profile.fft_backward, "Backward FFT", color6);
        tasks.emplace_back(profile.force_halo_exchange, "Exchange Force Halos", color1);
        tasks.emplace_back(profile.reassign_particles, "Reassign Particles", color2);

        simulation->profiler.LoadFrameData(tasks.data(), tasks.size(), profile.total_time);
      }
    } else if (st.MPI_TAG == TAG_DONE) {
      char dummy;
      MPI_Recv(&dummy, 1, MPI_BYTE, st.MPI_SOURCE, TAG_DONE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      busy = false;
      progress = 0.0;

      simulation->sync();
    }

    ImGui::PushItemWidth(80);
    ImGui::InputInt("Timesteps", &num_iterations, 0, 0);
    ImGui::PopItemWidth();
    ImGui::EndDisabled();
    ImGui::EndHorizontal();
    ImGui::EndDisabled();

    if (simulation) {
      ImGui::ProgressBar(progress, ImVec2(385, 0));
      if (!simulation->profiler.empty())
        simulation->profiler.Draw(100, 200);
    }

    ImGui::SeparatorText("Information");
    ImGui::LabelText("Number of Worker Nodes", "%d", wsize - 1);
    if (simulation) {
      ImGui::LabelText("Number of Particles", "%.2e", (double)simulation->N);
      ImGui::LabelText("Elapsed Time", "%e", simulation->t);
      if (params.USE_SCALE_FACTOR) {
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
      case Command::Shutdown:
        return;
      default:
        break;
    }
  }
}

int main(int argc, char** argv) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  fftw_init_threads();
  fftw_mpi_init();

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (size < 2) {
    printf("Needs at least two MPI processes to run!");
    exit(-1);
  }

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

  fftw_cleanup_threads();
  MPI_Finalize();
  return 0;
}
