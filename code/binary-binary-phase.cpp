#include "SpaceHub/src/spaceHub.hpp"
#include "SpaceX/SpaceHubWrapper.hpp"
using ConFile = space::multi_thread::ConcurrentFile;
using Particle = spacex::SpaceXsim::Particle;
using ParticleSys = spacex::SpaceXsim::RunArgs::ParticleSys;

USING_NAMESPACE_SPACEHUB_ALL;

double v_inf = 1_kms;

constexpr double Q_max = 100_AU;

constexpr double r_start = 1000_AU;

auto create_jupiter_system(double anomaly) {
  double a_jup = 5.2_AU;

  double m_jup = 1_Mj;

  double r_jup = 1_Rj;

  Particle sun{1_Ms, 1_Rs}, jupiter{m_jup, r_jup};

  auto jupiter_orbit = EllipOrbit{sun.mass, jupiter.mass, a_jup, 0.0, 0.0, 0.0, anomaly, 0.0};

  move_particles(jupiter_orbit, jupiter);

  move_to_COM_frame(sun, jupiter);

  return std::make_tuple(sun, jupiter, jupiter_orbit);
}

void binary_binary_Adrain_phase(size_t th_id, std::string const &dir, size_t sim_num, double a_s, double v_inf) {
  char name[100];
  sprintf(name, "%.0lf-%.2lf-%d.txt", a_s, v_inf / kms, static_cast<int>(th_id));

  std::fstream out_file{dir + "two-binaries-phase-" + name, std::fstream::out};

  for (size_t i = 0; i < sim_num; ++i) {
    double Q = random::Uniform(0, Q_max);

    double anomaly = random::Uniform(0, 2 * consts::pi);

    double b_phase = 0.5 * consts::pi * static_cast<int>(random::Uniform(0, 4));

    bool is_collided = false;

    int coll_i = -1;

    int coll_j = -1;

    Particle star1{1_Ms, 1_Rs}, star2{1_Ms, 1_Rs};

    auto [sun, jupiter, jupiter_orbit] = create_jupiter_system(anomaly);

    auto binary_orbit = EllipOrbit{star1.mass, star2.mass, a_s, 0.0, 0.0, 0.0, b_phase, 0.0};

    move_particles(binary_orbit, star2);

    auto const w = consts::pi * static_cast<int>(random::Uniform(0, 2));  // random::Uniform(0, 2 * consts::pi);

    double A = -M_tot(sun, jupiter, star1, star2) / (v_inf * v_inf);

    double E = 1 - Q / A;

    double B = -A * sqrt(E * E - 1);

    auto const in_orbit =
        orbit::HyperOrbit(M_tot(sun, jupiter), M_tot(star1, star2), v_inf, B, w, 0.0, 0.0, r_start, orbit::Hyper::in);

    move_particles(in_orbit, star1, star2);

    move_to_COM_frame(sun, jupiter, star1, star2);

    double end_time = 2 * time_to_periapsis(cluster(sun, jupiter), cluster(star1, star2));

    spacex::SpaceXsim::RunArgs args;

    args.add_stop_condition([&](auto &ptc, auto dt) -> bool {
      size_t number = ptc.number();
      for (size_t i = 0; i < number; ++i) {
        for (size_t j = i + 1; j < number; ++j) {
          if (ptc.Collide(i, j)) {
            is_collided = true;
            coll_i = static_cast<int>(i);
            coll_j = static_cast<int>(j);
            return true;
          }
        }
      }
      return false;
    });

    args.add_stop_condition(end_time);

    args.add_stop_point_operation([&](auto &ptc, auto dt) {
      space::display(out_file, i, w, Q, anomaly, is_collided, coll_i, coll_j, ptc, b_phase, "\r\n");
    });

    spacex::SpaceXsim simulator{0, sun, jupiter, star1, star2};

    simulator.run(args);
  }
}

int main(int argc, char **argv) {
  size_t sim_num;
  std::string output_dir;
  double a = 1;
  double v_inf = 1;

  tools::read_command_line(argc, argv, output_dir, sim_num, a, v_inf);

  v_inf *= kms;

  auto_indexed_multi_thread(binary_binary_Adrain_phase, output_dir, sim_num, a, v_inf);
}
