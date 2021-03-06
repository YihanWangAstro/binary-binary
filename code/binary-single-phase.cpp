#include "SpaceHub/src/spaceHub.hpp"
#include "SpaceX/SpaceHubWrapper.hpp"
using ConFile = space::multi_thread::ConcurrentFile;
using Particle = spacex::SpaceXsim::Particle;
using ParticleSys = spacex::SpaceXsim::RunArgs::ParticleSys;

USING_NAMESPACE_SPACEHUB_ALL;

double v_inf = 1_kms;

constexpr double Q_max = 40_AU;

constexpr double r_start = 500_AU;

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

double anomaly_at_peri(double finf, double t, double M, double a) {
  double T = 2 * consts::pi * sqrt(a * a * a / M);
  double df = 2 * consts::pi * (t / T - static_cast<int>(t / T));
  double f = finf + df;
  return 2 * consts::pi * (f / (2 * consts::pi) - static_cast<int>(f / (2 * consts::pi)));
}

void binary_single_Adrain_phase(size_t th_id, std::string const &dir, size_t sim_num, double v_inf) {
  char name[100];
  sprintf(name, "%.2lf-%d.txt", v_inf / kms, static_cast<int>(th_id));

  std::fstream phase_file{dir + "bs-out-" + name, std::fstream::out};

  // double A = Q / (1 - E);

  for (size_t i = 0; i < sim_num; ++i) {
    double Q = random::Uniform(-Q_max, Q_max);

    double anomaly = random::Uniform(0, 2 * consts::pi);

    bool is_collided = false;

    int coll_i = -1;

    int coll_j = -1;

    Particle star1{2_Ms, 1.68_Rs};

    auto [sun, jupiter, jupiter_orbit] = create_jupiter_system(anomaly);

    auto const w = 0.0 + static_cast<int>(Q < 0) * consts::pi;

    double A = -M_tot(sun, jupiter, star1) / (v_inf * v_inf);

    double E = 1 - fabs(Q) / A;

    // double V = sqrt(-M_tot(sun, jupiter, star1) / A);

    double B = -A * sqrt(E * E - 1);

    auto const in_orbit =
        orbit::HyperOrbit(M_tot(sun, jupiter), star1.mass, v_inf, B, w, 0.0, 0.0, r_start, orbit::Hyper::in);

    move_particles(in_orbit, star1);

    move_to_COM_frame(sun, jupiter, star1);

    double end_time = 2 * time_to_periapsis(cluster(sun, jupiter), star1);

    double anomaly_p = anomaly_at_peri(anomaly, end_time / 2, sun.mass + jupiter.mass, jupiter_orbit.a);

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
      auto out = scattering::to_hierarchical(ptc);
      space::display(phase_file, i, w, Q, anomaly, anomaly_p, is_collided, coll_i, coll_j, out, "\r\n");
    });

    spacex::SpaceXsim simulator{0, sun, jupiter, star1};

    simulator.run(args);
  }
}

int main(int argc, char **argv) {
  size_t sim_num;
  std::string output_dir;
  double v_inf = 1;
  tools::read_command_line(argc, argv, output_dir, sim_num, v_inf);
  v_inf *= kms;

  auto_indexed_multi_thread(binary_single_Adrain_phase, output_dir, sim_num, v_inf);
}
