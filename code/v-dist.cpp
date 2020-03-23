#include "SpaceHub/src/spaceHub.hpp"
#include "SpaceX/SpaceHubWrapper.hpp"
using ConFile = space::multi_thread::ConcurrentFile;
using Particle = spacex::SpaceXsim::Particle;
using ParticleSys = spacex::SpaceXsim::RunArgs::ParticleSys;

USING_NAMESPACE_SPACEHUB_ALL;

constexpr bool coll_detect = true;

constexpr double interact_factor = 0.001;

double v_inf = 0.1_kms;

// double num_den = 500;
/*
auto collision = [](auto &ptc, auto dt) -> bool {
  size_t number = ptc.number();
  for (size_t i = 0; i < number; ++i) {
    for (size_t j = i + 1; j < number; ++j) {
      if (ptc.Collide(i, j)) return true;
    }
  }
  return false;
};*/

auto post_ae = [](auto const &ptc) {
  auto const &star = ptc[0];

  auto const &jup = ptc[1];

  auto const &star1 = ptc[2];

  auto const &star2 = ptc[3];

  auto const [a_j, e_j] = calc_a_RL_vector(star.m + jup.m, jup.pos - star.pos, jup.vel - star.vel);

  auto const [a_s, e_s] = calc_a_RL_vector(star1.m + star2.m, star1.pos - star2.pos, star1.vel - star2.vel);

  return std::make_tuple(a_j, e_j, a_s, e_s);
};
/*
auto jupiter_ejection = [](auto const &ptc, auto dt) -> bool {
  auto const &star = ptc[0];

  auto const &jup = ptc[1];

  auto const [a_j, e_j] = calc_a_e(star.m + jup.m, jup.pos - star.pos, jup.vel - star.vel);

  auto R = norm(jup.pos - star.pos);

  if ((e_j < 0 || e_j > 1) && R > 1000_AU) {
    return true;
  } else {
    return false;
  }
};*/

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

inline double start_v(double m1, double m2, double v_inf, double r_start) {
  return sqrt(v_inf * v_inf + 2 * (m1 + m2) / r_start);
}

inline double fall_free_time(double m1, double m2, double r) {
  return 0.5 * consts::pi * sqrt(r * r * r / (2 * (m1 + m2)));
}

void binary_binary_Adrain(size_t th_id, std::string const &dir, size_t sim_num, double a_s) {
  char name[100];
  sprintf(name, "%.0lf-%d.txt", a_s, static_cast<int>(th_id));

  std::fstream out_file{dir + "fall-free-" + name, std::fstream::out};

  std::fstream state_file{dir + "fall-free-state-" + name, std::fstream::out};

  double const tidal_factor = 1e-5;

  double inc = 0.0;

  for (size_t i = 0; i < sim_num; ++i) {
    bool is_collided = false;

    Particle star1{1_Ms, 1_Rs}, star2{1_Ms, 1_Rs};

    auto [sun, jupiter, jupiter_orbit] = create_jupiter_system(inc);

    auto binary_orbit = EllipOrbit{star1.mass, star2.mass, a_s, 0, inc, isotherm, isotherm, isotherm};

    move_particles(binary_orbit, star2);

    auto const r_start = orbit::tidal_radius(tidal_factor, M_tot(sun, jupiter), M_tot(star1, star2),
                                             2 * jupiter_orbit.a, 2 * binary_orbit.a);

    double v_start = start_v(M_tot(sun, jupiter), M_tot(star1, star2), v_inf, r_start);

    move_particles(Vector3d(r_start, 0, 0), Vector3d(-v_start, 0, 0), star1, star2);

    move_to_COM_frame(sun, jupiter, star1, star2);

    double end_time = 20 * fall_free_time(M_tot(sun, jupiter), M_tot(star1, star2), r_start);

    spacex::SpaceXsim::RunArgs args;

    if constexpr (coll_detect) {
      args.add_stop_condition([&](auto &ptc, auto dt) -> bool {
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

    args.add_stop_point_operation([&](auto &ptc, auto dt) {
      auto [a_j, e_j, a_s, e_s] = post_ae(ptc);

      space::display(out_file, i, is_collided, a_j, e_j, a_s, e_s);
      out_file << std::endl;

      space::display(state_file, i, is_collided, ptc, "\r\n");
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

  std::vector<std::thread> threads;

  for (size_t th_id = 0; th_id < machine_thread_num; ++th_id) {
    threads.emplace_back(std::thread{binary_binary_Adrain, th_id, output_dir, sim_num, a});
  }

  for (auto &th : threads) {
    th.join();
  }
}
