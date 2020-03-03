#include "SpaceHub/src/spaceHub.hpp"
#include "SpaceX/SpaceHubWrapper.hpp"
using ConFile = space::multi_thread::ConcurrentFile;
using Particle = spacex::SpaceXsim::Particle;
using ParticleSys = spacex::SpaceXsim::RunArgs::ParticleSys;

USING_NAMESPACE_SPACEHUB_ALL;

constexpr bool coll_detect = true;

constexpr double interact_factor = 0.001;

double v_inf = 10_kms;

double num_den = 500;

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
/*
auto create_jupiter_system(double inc) {
  double a_jup = 5.2_AU;

  double m_jup = 1_Mj;

  double r_jup = 1_Rj;

  Particle sun{1_Ms, 1_Rs}, jupiter{m_jup, r_jup};

  auto jupiter_orbit = EllipOrbit{sun.mass, jupiter.mass, a_jup, 0, inc, 0, 0, random::Uniform(0, 2 * consts::pi)};

  move_particles(jupiter_orbit, jupiter);

  move_to_COM_frame(sun, jupiter);

  return std::make_tuple(sun, jupiter, jupiter_orbit);
}*/

auto create_iso_jupiter_system() {
  double a_jup = 5.2_AU;

  double m_jup = 1_Mj;

  double r_jup = 1_Rj;

  Particle sun{1_Ms, 1_Rs}, jupiter{m_jup, r_jup};

  auto jupiter_orbit = EllipOrbit{sun.mass, jupiter.mass, a_jup, 0, isotherm, isotherm, isotherm, isotherm};

  move_particles(jupiter_orbit, jupiter);

  move_to_COM_frame(sun, jupiter);

  return std::make_tuple(sun, jupiter, jupiter_orbit);
}

void binary_binary(size_t th_id, std::string const &dir, size_t sim_num, double a_s, double b_factor) {
  char name[100];
  sprintf(name, "%.0lf-%.0lf.txt", a_s, 8 * b_factor);

  std::fstream out_file{dir + "two-binaries-" + name, std::fstream::out};

  double const tidal_factor = 1e-5;

  for (size_t i = 0; i < sim_num; ++i) {
    Particle star1{1_Ms, 1_Rs}, star2{1_Ms, 1_Rs};

    auto [sun, jupiter, jupiter_orbit] = create_iso_jupiter_system();

    auto binary_orbit = EllipOrbit{star1.mass, star2.mass, a_s, 0, isotherm, isotherm, isotherm, isotherm};

    move_particles(binary_orbit, star2);

    auto const b = (jupiter_orbit.a + 0.5 * binary_orbit.a) * b_factor;

    auto const w = random::Uniform(0, 2 * consts::pi);

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

      space::display(out_file, i, ptc, a_j, e_j, a_s, e_s);
      out_file << std::endl;
    });

    spacex::SpaceXsim simulator{0, sun, jupiter, star1, star2};

    simulator.run(args);
  }
}

int main(int argc, char **argv) {
  size_t sim_num;
  std::string output_dir;

  tools::read_command_line(argc, argv, output_dir, sim_num);

  // std::array<double, 1> incs = {0_deg};

  std::array<double, 4> a_s = {1_AU, 5_AU, 25_AU, 125_AU};
  std::array<double, 7> b_factors = {0.125, 0.25, 0.5, 1, 2, 4, 8};

  std::vector<std::thread> threads;

  size_t th_id = 0;

  for (auto a : a_s) {
    for (auto b : b_factors) {
      threads.emplace_back(std::thread{binary_binary, th_id, output_dir, sim_num, a, b});
      th_id++;
    }
  }

  for (auto &th : threads) {
    th.join();
  }
}
