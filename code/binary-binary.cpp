#include "SpaceHub/src/spaceHub.hpp"
#include "SpaceX/SpaceHubWrapper.hpp"
using ConFile = space::multi_thread::ConcurrentFile;
using Particle = spacex::SpaceXsim::Particle;
using ParticleSys = spacex::SpaceXsim::RunArgs::ParticleSys;

USING_NAMESPACE_SPACEHUB_ALL;

constexpr bool coll_detect = true;

constexpr double interact_factor = 0.001;

double v_inf = 10_kms;

// double num_den = 500;

auto collision = [](auto &ptc) -> bool {
  size_t number = ptc.number();
  for (size_t i = 0; i < number; ++i) {
    for (size_t j = i + 1; j < number; ++j) {
      if (ptc.Collide(i, j)) return true;
    }
  }
  return false;
};

auto post_ae = [](auto const &ptc) {
  auto const &star = ptc[0];

  auto const &jup = ptc[1];

  auto const &star1 = ptc[2];

  auto const &star2 = ptc[3];

  auto const [a_j, e_j] = calc_a_RL_vector(star.m + jup.m, jup.pos - star.pos, jup.vel - star.vel);

  auto const [a_s, e_s] = calc_a_RL_vector(star1.m + star2.m, star1.pos - star2.pos, star1.vel - star2.vel);

  return std::make_tuple(a_j, e_j, a_s, e_s);
};

auto jupiter_ejection = [](auto const &ptc) -> bool {
  auto const &star = ptc[0];

  auto const &jup = ptc[1];

  auto const [a_j, e_j] = calc_a_e(star.m + jup.m, jup.pos - star.pos, jup.vel - star.vel);

  auto R = norm(jup.pos - star.pos);

  if ((e_j < 0 || e_j > 1) && R > 1000_AU) {
    return true;
  } else {
    return false;
  }
};

auto create_jupiter_system(double inc) {
  double a_jup = 5.2_AU;

  double m_jup = 1_Mj;

  double r_jup = 1_Rj;

  Particle sun{1_Ms, 1_Rs}, jupiter{m_jup, r_jup};

  auto jupiter_orbit = EllipOrbit{sun.mass, jupiter.mass, a_jup, 0, inc, isotherm, isotherm, isotherm};

  move_particles(jupiter_orbit, jupiter);

  move_to_COM_frame(sun, jupiter);

  return std::make_tuple(sun, jupiter, jupiter_orbit);
}
/*
auto create_iso_jupiter_system() {
  double a_jup = 5.2_AU;

  double m_jup = 1_Mj;

  double r_jup = 1_Rj;

  Particle sun{1_Ms, 1_Rs}, jupiter{m_jup, r_jup};

  auto jupiter_orbit = EllipOrbit{sun.mass, jupiter.mass, a_jup, 0, isotherm, isotherm, isotherm, isotherm};

  move_particles(jupiter_orbit, jupiter);

  move_to_COM_frame(sun, jupiter);

  return std::make_tuple(sun, jupiter, jupiter_orbit);
}*/

void binary_binary_Adrain(size_t th_id, std::string const &dir, size_t sim_num, double a_s, double E, double Q) {
  char name[100];
  sprintf(name, "%.0lf-%.2lf-%d.txt", a_s, E, static_cast<int>(th_id));

  std::fstream out_file{dir + "two-binaries-" + name, std::fstream::out};

  std::fstream state_file{dir + "state-" + name, std::fstream::out};

  double const tidal_factor = 1e-7;

  double inc = 0.0;

  double A = Q / (1 - E);

  for (size_t i = 0; i < sim_num; ++i) {
    bool is_collided = false;

    Particle star1{1_Ms, 1_Rs}, star2{1_Ms, 1_Rs};

    auto [sun, jupiter, jupiter_orbit] = create_jupiter_system(inc);

    auto binary_orbit = EllipOrbit{star1.mass, star2.mass, a_s, 0, inc, isotherm, isotherm, isotherm};

    move_particles(binary_orbit, star2);

    auto const w = consts::pi * static_cast<int>(random::Uniform(0, 2));  // random::Uniform(0, 2 * consts::pi);

    auto const r_start = std::max(orbit::tidal_radius(tidal_factor, M_tot(sun, jupiter), M_tot(star1, star2),
                                                      2 * jupiter_orbit.a, 2 * binary_orbit.a),
                                  5 * Q);

    double V = sqrt(-M_tot(sun, jupiter, star1, star2) / A);

    double B = -A * sqrt(E * E - 1);

    auto const in_orbit =
        orbit::HyperOrbit(M_tot(sun, jupiter), M_tot(star1, star2), V, B, w, 0, 0, r_start, orbit::Hyper::in);

    move_particles(in_orbit, star1, star2);

    move_to_COM_frame(sun, jupiter, star1, star2);

    double end_time = 2 * time_to_periapsis(cluster(sun, jupiter), cluster(star1, star2));

    spacex::SpaceXsim::RunArgs args;

    if constexpr (coll_detect) {
      args.add_stop_condition([&](auto &ptc) -> bool {
        size_t number = ptc.number();
        for (size_t i = 0; i < number; ++i) {
          for (size_t j = i + 1; j < number; ++j) {
            if (ptc.Collide(i, j)) {
              is_collided = true;
              return true;
            }
          }
        }
        return false;
      });
    }

    args.add_stop_condition(end_time);

    args.add_stop_point_operation([&](auto &ptc) {
      auto [a_j, e_j, a_s, e_s] = post_ae(ptc);

      space::display(out_file, i, w, is_collided, a_j, e_j, a_s, e_s);
      out_file << std::endl;

      space::display(state_file, i, w, is_collided, ptc, "\r\n");
    });

    spacex::SpaceXsim simulator{0, sun, jupiter, star1, star2};

    simulator.run(args);
  }
}

void binary_binary(size_t th_id, std::string const &dir, size_t sim_num, double a_s, double b_factor) {
  char name[100];
  sprintf(name, "%.0lf-%.0lf.txt", a_s, 2 * b_factor);

  std::fstream out_file{dir + "two-binaries-" + name, std::fstream::out};

  double const tidal_factor = 1e-5;

  double inc = 0.0;

  for (size_t i = 0; i < sim_num; ++i) {
    Particle star1{1_Ms, 1_Rs}, star2{1_Ms, 1_Rs};

    auto [sun, jupiter, jupiter_orbit] = create_jupiter_system(inc);

    auto binary_orbit = EllipOrbit{star1.mass, star2.mass, a_s, 0, inc, isotherm, isotherm, isotherm};

    move_particles(binary_orbit, star2);

    auto const b = (jupiter_orbit.a + 0.5 * binary_orbit.a) * b_factor;

    auto const w = consts::pi * static_cast<int>(random::Uniform(0, 2));  // random::Uniform(0, 2 * consts::pi);

    auto const r_start = orbit::tidal_radius(tidal_factor, M_tot(sun, jupiter), M_tot(star1, star2),
                                             2 * jupiter_orbit.a, 2 * binary_orbit.a);

    auto const in_orbit =
        orbit::HyperOrbit(M_tot(sun, jupiter), M_tot(star1, star2), v_inf, b, w, 0, 0, r_start, orbit::Hyper::in);

    move_particles(in_orbit, star1, star2);

    move_to_COM_frame(sun, jupiter, star1, star2);

    double end_time = 2 * time_to_periapsis(cluster(sun, jupiter), cluster(star1, star2));

    spacex::SpaceXsim::RunArgs args;

    if constexpr (coll_detect) {
      args.add_stop_condition(collision);
    }

    args.add_stop_condition(end_time);

    args.add_stop_point_operation([&](auto &ptc) {
      auto [a_j, e_j, a_s, e_s] = post_ae(ptc);

      space::display(out_file, i, w, a_j, e_j, a_s, e_s);
      out_file << std::endl;
    });

    spacex::SpaceXsim simulator{0, sun, jupiter, star1, star2};

    simulator.run(args);
  }
}

int main(int argc, char **argv) {
  size_t sim_num;
  std::string output_dir;
  double a = 1;
  tools::read_command_line(argc, argv, output_dir, sim_num, a);

  // std::array<double, 1> incs = {0_deg};
  // std::array<double, 4> a_s = {1_AU, 5_AU, 25_AU, 125_AU};
  // std::array<double, 7> b_factors = {0.125, 0.25, 0.5, 1, 2, 4, 8};
  std::vector<std::thread> threads;

  size_t th_id = 0;
  double min_Q = 1_AU;
  double max_Q = 1000_AU;
  double E = 1.5;

  double dp = (log10(max_Q) - log10(min_Q)) / 40;

  for (double p = 0; p <= log10(max_Q); p += dp) {
    threads.emplace_back(std::thread{binary_binary_Adrain, th_id, output_dir, sim_num, a, E, pow(10.0, p)});
    th_id++;
  }

  // double b_factor_max = 20;
  // double db = 0.5;
  // for (auto a : a_s) {
  // for (auto b : b_factors) {
  /*for (double b = 0.5; b <= b_factor_max; b += db) {
    threads.emplace_back(std::thread{binary_binary, th_id, output_dir, sim_num, a, b});
    th_id++;
  }*/
  //}

  for (auto &th : threads) {
    th.join();
  }
}
