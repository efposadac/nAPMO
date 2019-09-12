/*file: libint_iface.cpp
nAPMO package
Copyright (c) 2016, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co
*/

#include "libint_iface.h"
#include <libint2/util/small_vector.h>

/*
LibintInterface
*/
LibintInterface::LibintInterface(const int id)
    : max_nprim(0), nbasis(0), max_l(0), sID(id) {

  // set up thread pool
  using napmo::nthreads;
  set_nthreads();

  // initializes the Libint integrals library ... now ready to compute
  libint2::initialize();
}

void LibintInterface::init_2body_ints() {

  // compute OBS non-negligible shell-pair list

  obs_shellpair_list = compute_shellpair_list(shells);
  size_t nsp = 0;
  for (auto &sp : obs_shellpair_list) {
    nsp += sp.second.size();
  }

  // std::cout << "# of {all,non-negligible} shell-pairs = {"
  //           << shells.size() * (shells.size() + 1) / 2 << "," << nsp << "}"
  //           << std::endl;
}

void LibintInterface::add_pointcharges(const int z, const double *center) {

  // add atom information
  libint2::Atom aux;
  aux.atomic_number = z;
  aux.x = center[0];
  aux.y = center[1];
  aux.z = center[2];
  point_charges.push_back(aux);

  // std::cout<<aux.atomic_number<<std::endl;
  // std::cout<<aux.x<<" "<<aux.y<<" "<<aux.z<<" "<<std::endl;
}

void LibintInterface::add_basis(BasisSet *basis) {

  // add basis-set

  nbasis += basis->get_nbasis();
  max_nprim = std::max(max_nprim, basis->get_max_nprim());
  max_l = std::max(max_l, basis->get_max_l());

  // std::cout<<nbasis<<' '<<max_nprim<<' '<<max_l<<std::endl;

  libint2::Shell::do_enforce_unit_normalization(false);

  for (auto cont : basis->get_cont()) {

    if ((cont.get_l()[0] + cont.get_l()[1] + cont.get_l()[2]) !=
        cont.get_l()[0])
      continue;

    // std::vector<double> exponents;
    // std::vector<double> coefficients;
    boost::container::small_vector<double,6> exponents;
    boost::container::small_vector<double,6> coefficients;

    for (auto prim : cont.get_prim()) {
      exponents.push_back(prim.get_zeta());
      coefficients.push_back(prim.get_coeff());
    }

    shells.push_back(
        {{exponents},
         {
             {cont.get_l()[0], false, {coefficients}},
         },
         {{cont.get_origin()[0], cont.get_origin()[1], cont.get_origin()[2]}}});

    // Renormalize
    const auto &shell = shells.back();

    // std::cout << shell << std::endl;

    libint2::Engine engine(libint2::Operator::overlap, shell.nprim(), max_l, 0);
    const auto &buf = engine.results();
    engine.compute(shell, shell);
    Eigen::Map<const Matrix> buf_mat(buf[0], shell.size(), shell.size());
    for (unsigned int i = 0; i < shell.size(); ++i) {
      norma.push_back(1 / sqrt(buf_mat(i, i)));
    }
  }
}

std::vector<size_t> LibintInterface::map_shell_to_basis_function() {
  std::vector<size_t> result;
  result.reserve(shells.size());

  size_t n = 0;
  for (auto shell : shells) {
    result.push_back(n);
    n += shell.size();
  }

  return result;
}

Matrix LibintInterface::compute_shellblock_norm(const Matrix &A) {
  const auto nsh = shells.size();
  Matrix Ash(nsh, nsh);

  auto shell2bf = map_shell_to_basis_function();
  for (size_t s1 = 0; s1 != nsh; ++s1) {
    const auto &s1_first = shell2bf[s1];
    const auto &s1_size = shells[s1].size();
    for (size_t s2 = 0; s2 != nsh; ++s2) {
      const auto &s2_first = shell2bf[s2];
      const auto &s2_size = shells[s2].size();

      Ash(s1, s2) = A.block(s1_first, s2_first, s1_size, s2_size)
                        .lpNorm<Eigen::Infinity>();
    }
  }

  return Ash;
}

Matrix LibintInterface::compute_1body_ints(libint2::Operator obtype) {
  const auto n = nbasis;
  Matrix result(n, n);

  // construct the integrals engine
  libint2::Engine engine(obtype, max_nprim, max_l, 0);

  // nuclear attraction ints engine needs to know where the charges sit ...
  // the nuclei are charges in this case
  if (obtype == libint2::Operator::nuclear) {
    std::vector<std::pair<double, std::array<double, 3>>> q;
    for (const auto &point : point_charges) {
      q.push_back({static_cast<double>(point.atomic_number),
                   {{point.x, point.y, point.z}}});
    }
    engine.set_params(q);
  }

  auto shell2bf = map_shell_to_basis_function();
  const auto &buf = engine.results();

  // loop over unique shell pairs, {s1,s2} such that s1 >= s2
  // this is due to the permutational symmetry of the real integrals over
  // Hermitian operators: (1|2) = (2|1)
  for (unsigned int s1 = 0; s1 != shells.size(); ++s1) {
    auto bf1 = shell2bf[s1]; // first basis function in this shell
    auto n1 = shells[s1].size();

    for (unsigned int s2 = 0; s2 <= s1; ++s2) {
      auto bf2 = shell2bf[s2];
      auto n2 = shells[s2].size();

      // compute shell pair; return is the pointer to the buffer
      engine.compute(shells[s1], shells[s2]);

      // "map" buffer to a const Eigen Matrix, and copy it to the
      // corresponding
      // blocks of the result
      Eigen::Map<const Matrix> buf_mat(buf[0], n1, n2);
      result.block(bf1, bf2, n1, n2) = buf_mat;
      if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1}
                    // block, note the transpose!
        result.block(bf2, bf1, n2, n1) = buf_mat.transpose();
    }
  }

  // normalize
  for (int i = 0; i < result.cols(); ++i) {
    for (int j = 0; j < result.rows(); ++j) {
      result(i, j) *= norma[i] * norma[j];
    }
  }

  return result;
}

std::vector<QuartetBuffer> *
LibintInterface::compute_2body_ints(const Matrix &D, const Matrix &Schwartz,
                                    double precision) {

  using napmo::nthreads;
  using libint2::Engine;

  const auto nshells = shells.size();

  auto n = (nbasis * (nbasis + 1)) / 2;
  n = (n * (n + 1)) / 2;
  n = ceil(n / nthreads) + 1;

  // std::cout<<"Size per buffer: "<<n<<std::endl;

  const auto do_schwartz_screen = Schwartz.cols() != 0 && Schwartz.rows() != 0;
  Matrix D_shblk_norm; // matrix of infty-norms of shell blocks
  if (do_schwartz_screen) {
    D_shblk_norm = compute_shellblock_norm(D);
  }

  auto fock_precision = precision;

  // engine precision controls primitive truncation, assume worst-casescenario
  // (all primitive combinations add up constructively)
  auto max_nprim4 = max_nprim * max_nprim * max_nprim * max_nprim;
  auto engine_precision = std::min(fock_precision / D_shblk_norm.maxCoeff(),
                                   std::numeric_limits<double>::epsilon()) /
                          max_nprim4;

  // construct the 2-electron repulsion integrals engine pool
  std::vector<Engine> engines(nthreads);
  engines[0] = Engine(libint2::Operator::coulomb, max_nprim, max_l, 0);
  engines[0].set_precision(engine_precision); // shellset-dependentprecision
                                              // control will likely break
                                              // positive definiteness
                                              // stick with this simple recipe

  // std::cout << "compute_2body_direct:precision = " << precision << std::endl;
  // std::cout << "Engine::precision = " << engines[0].precision() <<
  // std::endl;

  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::atomic<size_t> num_ints_computed{0};

  // Setting buffers up
  std::vector<QuartetBuffer> *buffers = 0;
  buffers = new std::vector<QuartetBuffer>(nthreads);

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = map_shell_to_basis_function();

  auto lambda = [&](unsigned int thread_id) {

    auto &engine = engines[thread_id];

    buffers->at(thread_id).p.resize(n);
    buffers->at(thread_id).q.resize(n);
    buffers->at(thread_id).r.resize(n);
    buffers->at(thread_id).s.resize(n);
    buffers->at(thread_id).val.resize(n);

    const auto &buf = engines[thread_id].results();

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto &timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    auto counter = 0;
    // loop over permutationally-unique set of shells
    for (unsigned int s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = shells[s1].size();   // number of basis functions in this shell

      for (const auto &s2 : obs_shellpair_list[s1]) {
        auto bf2_first = shell2bf[s2];
        auto n2 = shells[s2].size();

        const auto Dnorm12 = do_schwartz_screen ? D_shblk_norm(s1, s2) : 0.;

        for (unsigned int s3 = 0; s3 <= s1; ++s3) {
          auto bf3_first = shell2bf[s3];
          auto n3 = shells[s3].size();

          const auto Dnorm123 =
              do_schwartz_screen
                  ? std::max(D_shblk_norm(s1, s3),
                             std::max(D_shblk_norm(s2, s3), Dnorm12))
                  : 0.;

          const auto s4_max = (s1 == s3) ? s2 : s3;
          for (const auto &s4 : obs_shellpair_list[s3]) {
            if (s4 > s4_max)
              break; // for each s3, s4 are stored in monotonically increasing
                     // order

            if ((s1234++) % nthreads != thread_id)
              continue;

            const auto Dnorm1234 =
                do_schwartz_screen
                    ? std::max(
                          D_shblk_norm(s1, s4),
                          std::max(D_shblk_norm(s2, s4),
                                   std::max(D_shblk_norm(s3, s4), Dnorm123)))
                    : 0.;

            if (do_schwartz_screen &&
                Dnorm1234 * Schwartz(s1, s2) * Schwartz(s3, s4) <
                    fock_precision)
              continue;

            auto bf4_first = shell2bf[s4];
            auto n4 = shells[s4].size();

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx,
                            0>(shells[s1], shells[s2], shells[s3], shells[s4]);

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif

            if (buf[0] == nullptr)
              continue; // all screened out

            auto intIter = IntraIntsIt(n1, n2, n3, n4, bf1_first, bf2_first,
                                       bf3_first, bf4_first);

            for (intIter.first(); intIter.is_done() == false; intIter.next()) {
              if (std::abs(buf[0][intIter.index()]) > 1.0e-10) {
                buffers->at(thread_id).p[counter] = intIter.i();
                buffers->at(thread_id).q[counter] = intIter.j();
                buffers->at(thread_id).r[counter] = intIter.k();
                buffers->at(thread_id).s[counter] = intIter.l();
                buffers->at(thread_id).val[counter] =
                    buf[0][intIter.index()] * norma[intIter.i()] *
                    norma[intIter.j()] * norma[intIter.k()] *
                    norma[intIter.l()];

                ++num_ints_computed;
                ++counter;
              }
            }
          }
        }
      }
    }

    buffers->at(thread_id).p.shrink_to_fit();
    buffers->at(thread_id).q.shrink_to_fit();
    buffers->at(thread_id).r.shrink_to_fit();
    buffers->at(thread_id).s.shrink_to_fit();
    buffers->at(thread_id).val.shrink_to_fit();

  }; // end of lambda

  napmo::parallel_do(lambda);

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto &t : timers) {
    time_for_ints += t.read(0);
  }
  std::cout << " Time for intras-pecies integrals = " << time_for_ints
            << std::endl;
  for (int t = 0; t != nthreads; ++t)
    engines[t].print_timers();
#endif

  // std::cout << " Number of unique integrals for species: " << sID << " = " << num_ints_computed << std::endl;

  return buffers;
}

void LibintInterface::compute_2body_disk(const char *filename, const Matrix &D,
                                         const Matrix &Schwartz,
                                         double precision) {
  const auto nshells = shells.size();

  using napmo::nthreads;

  const auto do_schwartz_screen = Schwartz.cols() != 0 && Schwartz.rows() != 0;
  Matrix D_shblk_norm; // matrix of infty-norms of shell blocks
  if (do_schwartz_screen) {
    D_shblk_norm = compute_shellblock_norm(D);
  }

  auto fock_precision = precision;

  // engine precision controls primitive truncation, assume worst-casescenario
  // (all primitive combinations add up constructively)
  auto max_nprim4 = max_nprim * max_nprim * max_nprim * max_nprim;
  auto engine_precision = std::min(fock_precision / D_shblk_norm.maxCoeff(),
                                   std::numeric_limits<double>::epsilon()) /
                          max_nprim4;

  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  engines[0] = Engine(libint2::Operator::coulomb, max_nprim, max_l, 0);
  engines[0].set_precision(engine_precision); // shellset-dependentprecision
                                              // control will likely break
                                              // positive definiteness
                                              // stick with this simple recipe

  // std::cout << "compute_2body_disk:precision = " << precision << std::endl;
  // std::cout << "Engine::precision = " << engines[0].precision() << std::endl;

  for (unsigned int i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::atomic<size_t> num_ints_computed{0};

  // setting buffers
  std::vector<QuartetBuffer> buffers(nthreads);

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = map_shell_to_basis_function();

  auto lambda = [&](unsigned int thread_id) {

    std::string file = std::to_string(thread_id);
    file += filename;

    std::ofstream outfile;
    outfile.open(file, std::ios::binary);

    const auto &buf = engines[thread_id].results();
    auto &engine = engines[thread_id];
    auto &buffer = buffers[thread_id];

    buffer.p.reserve(STACK_SIZE);
    buffer.q.reserve(STACK_SIZE);
    buffer.r.reserve(STACK_SIZE);
    buffer.s.reserve(STACK_SIZE);
    buffer.val.reserve(STACK_SIZE);

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto &timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    // loop over permutationally-unique set of shells
    int counter = 0;
    for (unsigned int s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = shells[s1].size();   // number of basis functions in this shell

      for (const auto &s2 : obs_shellpair_list[s1]) {
        auto bf2_first = shell2bf[s2];
        auto n2 = shells[s2].size();

        const auto Dnorm12 = do_schwartz_screen ? D_shblk_norm(s1, s2) : 0.;

        for (unsigned int s3 = 0; s3 <= s1; ++s3) {
          auto bf3_first = shell2bf[s3];
          auto n3 = shells[s3].size();

          const auto Dnorm123 =
              do_schwartz_screen
                  ? std::max(D_shblk_norm(s1, s3),
                             std::max(D_shblk_norm(s2, s3), Dnorm12))
                  : 0.;

          const auto s4_max = (s1 == s3) ? s2 : s3;
          for (const auto &s4 : obs_shellpair_list[s3]) {
            if (s4 > s4_max)
              break; // for each s3, s4 are stored in monotonically increasing
                     // order

            if ((s1234++) % nthreads != thread_id)
              continue;

            const auto Dnorm1234 =
                do_schwartz_screen
                    ? std::max(
                          D_shblk_norm(s1, s4),
                          std::max(D_shblk_norm(s2, s4),
                                   std::max(D_shblk_norm(s3, s4), Dnorm123)))
                    : 0.;

            if (do_schwartz_screen &&
                Dnorm1234 * Schwartz(s1, s2) * Schwartz(s3, s4) <
                    fock_precision)
              continue;

            auto bf4_first = shell2bf[s4];
            auto n4 = shells[s4].size();

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx,
                            0>(shells[s1], shells[s2], shells[s3], shells[s4]);

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif
            if (buf[0] == nullptr)
              continue; // all screened out

            auto intIter = IntraIntsIt(n1, n2, n3, n4, bf1_first, bf2_first,
                                       bf3_first, bf4_first);

            for (intIter.first(); intIter.is_done() == false; intIter.next()) {
              if (std::abs(buf[0][intIter.index()]) > 1.0e-10) {
                buffer.p[counter] = intIter.i() + 1;
                buffer.q[counter] = intIter.j() + 1;
                buffer.r[counter] = intIter.k() + 1;
                buffer.s[counter] = intIter.l() + 1;
                buffer.val[counter] = buf[0][intIter.index()] *
                                      norma[intIter.i()] * norma[intIter.j()] *
                                      norma[intIter.k()] * norma[intIter.l()];

                ++num_ints_computed;
                ++counter;
              }

              if (counter == STACK_SIZE) {
                write_buffer(buffer, outfile);
                counter = 0;
              }
            }
          }
        }
      }
    }

    buffer.p[counter] = -1;
    write_buffer(buffer, outfile);
    outfile.close();

  }; // end of lambda

  napmo::parallel_do(lambda);

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto &t : timers) {
    time_for_ints += t.read(0);
  }
  std::cout << " Time for intra-species integrals = " << time_for_ints
            << std::endl;
  for (int t = 0; t != nthreads; ++t)
    engines[t].print_timers();
#endif

  // std::cout << " Number of unique integrals for species: " << sID << " = "
  //           << num_ints_computed << std::endl;
}

Matrix LibintInterface::compute_2body_direct(const Matrix &D,
                                             const Matrix &Schwartz,
                                             double precision) {
  const auto n = nbasis;
  const auto nshells = shells.size();

  using napmo::nthreads;

  std::vector<Matrix> G(nthreads, Matrix::Zero(n, n));

  const auto do_schwartz_screen = Schwartz.cols() != 0 && Schwartz.rows() != 0;
  Matrix D_shblk_norm; // matrix of infty-norms of shell blocks
  if (do_schwartz_screen) {
    D_shblk_norm = compute_shellblock_norm(D);
  }

  auto fock_precision = precision;

  // engine precision controls primitive truncation, assume worst-casescenario
  // (all primitive combinations add up constructively)
  auto max_nprim4 = max_nprim * max_nprim * max_nprim * max_nprim;
  auto engine_precision = std::min(fock_precision / D_shblk_norm.maxCoeff(),
                                   std::numeric_limits<double>::epsilon()) /
                          max_nprim4;

  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  engines[0] = Engine(libint2::Operator::coulomb, max_nprim, max_l, 0);
  engines[0].set_precision(engine_precision); // shellset-dependentprecision
                                              // control will likely break
                                              // positive definiteness
                                              // stick with this simple recipe

  // std::cout << "compute_2body_direct:precision = " << precision << std::endl;
  // std::cout << "Engine::precision = " << engines[0].precision() <<
  // std::endl;

  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::atomic<size_t> num_ints_computed{0};

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = map_shell_to_basis_function();

  auto lambda = [&](unsigned int thread_id) {

    auto &engine = engines[thread_id];
    auto &g = G[thread_id];
    const auto &buf = engines[thread_id].results();

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto &timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    // loop over permutationally-unique set of shells
    for (unsigned int s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = shells[s1].size();   // number of basis functions in this shell

      for (const auto &s2 : obs_shellpair_list[s1]) {
        auto bf2_first = shell2bf[s2];
        auto n2 = shells[s2].size();

        const auto Dnorm12 = do_schwartz_screen ? D_shblk_norm(s1, s2) : 0.;

        for (unsigned int s3 = 0; s3 <= s1; ++s3) {
          auto bf3_first = shell2bf[s3];
          auto n3 = shells[s3].size();

          const auto Dnorm123 =
              do_schwartz_screen
                  ? std::max(D_shblk_norm(s1, s3),
                             std::max(D_shblk_norm(s2, s3), Dnorm12))
                  : 0.;

          const auto s4_max = (s1 == s3) ? s2 : s3;
          for (const auto &s4 : obs_shellpair_list[s3]) {
            if (s4 > s4_max)
              break; // for each s3, s4 are stored in monotonically increasing
                     // order

            if ((s1234++) % nthreads != thread_id)
              continue;

            const auto Dnorm1234 =
                do_schwartz_screen
                    ? std::max(
                          D_shblk_norm(s1, s4),
                          std::max(D_shblk_norm(s2, s4),
                                   std::max(D_shblk_norm(s3, s4), Dnorm123)))
                    : 0.;

            if (do_schwartz_screen &&
                Dnorm1234 * Schwartz(s1, s2) * Schwartz(s3, s4) <
                    fock_precision)
              continue;

            auto bf4_first = shell2bf[s4];
            auto n4 = shells[s4].size();

            num_ints_computed += n1 * n2 * n3 * n4;

            // compute the permutational degeneracy (i.e. # of equivalents) of
            // the given shell set
            auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
            auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
            auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
            auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx,
                            0>(shells[s1], shells[s2], shells[s3], shells[s4]);

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif

            if (buf[0] == nullptr)
              continue; // all screened out

            for (unsigned int f1 = 0, f1234 = 0; f1 != n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for (unsigned int f2 = 0; f2 != n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for (unsigned int f3 = 0; f3 != n3; ++f3) {
                  const auto bf3 = f3 + bf3_first;
                  for (unsigned int f4 = 0; f4 != n4; ++f4, ++f1234) {
                    const auto bf4 = f4 + bf4_first;

                    const auto value = buf[0][f1234];

                    const auto value_scal_by_deg = value * s1234_deg *
                                                   norma[bf1] * norma[bf2] *
                                                   norma[bf3] * norma[bf4];

                    g(bf1, bf2) += D(bf3, bf4) * value_scal_by_deg;
                    g(bf3, bf4) += D(bf1, bf2) * value_scal_by_deg;
                    g(bf1, bf3) -= 0.25 * D(bf2, bf4) * value_scal_by_deg;
                    g(bf2, bf4) -= 0.25 * D(bf1, bf3) * value_scal_by_deg;
                    g(bf1, bf4) -= 0.25 * D(bf2, bf3) * value_scal_by_deg;
                    g(bf2, bf3) -= 0.25 * D(bf1, bf4) * value_scal_by_deg;
                  }
                }
              }
            }
          }
        }
      }
    }

  }; // end of lambda

  napmo::parallel_do(lambda);

  // accumulate contributions from all threads
  for (size_t i = 1; i != nthreads; ++i) {
    G[0] += G[i];
  }

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto &t : timers) {
    time_for_ints += t.read(0);
  }
  std::cout << " Time for intras-pecies integrals = " << time_for_ints
            << std::endl;
  for (int t = 0; t != nthreads; ++t)
    engines[t].print_timers();
#endif

  // std::cout << " Number of unique integrals for species: " << sID << "
  // = "
  //           << num_ints_computed << std::endl;
  // symmetrize the result and return
  // symmetrize the result and return

  Matrix GG = 0.25 * (G[0] + G[0].transpose());
  return GG;
}

void LibintInterface::compute_coupling_disk(LibintInterface &other,
                                            const char *filename,
                                            double precision) {
  const auto nshells = shells.size();
  const auto oshells = other.get_shells();
  const auto onshells = oshells.size();

  using napmo::nthreads;

  auto fock_precision = precision;

  // engine precision controls primitive truncation, assume worst-casescenario
  // (all primitive combinations add up constructively)
  auto nprim_max = std::max(max_nprim, other.max_nprim);
  auto l_max = std::max(max_l, other.max_l);

  auto max_nprim4 = nprim_max * nprim_max * nprim_max * nprim_max;

  auto engine_precision =
      std::min(fock_precision, std::numeric_limits<double>::epsilon()) /
      max_nprim4;

  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  engines[0] = Engine(libint2::Operator::coulomb, nprim_max, l_max, 0);
  engines[0].set_precision(engine_precision); // shellset-dependentprecision
                                              // control will likely break
                                              // positive definiteness
                                              // stick with this simple recipe

  // std::cout << "compute_coupling_direct:precision = " << precision <<
  // std::endl;
  // std::cout << "Engine::precision = " << engines[0].precision() << std::endl;

  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::atomic<size_t> num_ints_computed{0};

  // setting buffers
  std::vector<QuartetBuffer> buffers(nthreads);

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = map_shell_to_basis_function();
  auto oshell2bf = other.map_shell_to_basis_function();

  auto lambda = [&](unsigned int thread_id) {

    std::string file = std::to_string(thread_id);
    file += filename;

    std::ofstream outfile;
    outfile.open(file, std::ios::binary);

    auto &engine = engines[thread_id];
    const auto &buf = engines[thread_id].results();
    auto &buffer = buffers[thread_id];

    buffer.p.reserve(STACK_SIZE);
    buffer.q.reserve(STACK_SIZE);
    buffer.r.reserve(STACK_SIZE);
    buffer.s.reserve(STACK_SIZE);
    buffer.val.reserve(STACK_SIZE);

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto &timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    // loop over permutationally-unique set of shells
    int counter = 0;
    for (unsigned int s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = shells[s1].size();   // number of basis functions in this shell

      for (unsigned int s2 = s1; s2 != nshells; ++s2) {
        auto bf2_first = shell2bf[s2];
        auto n2 = shells[s2].size();

        for (unsigned int os1 = 0l; os1 != onshells; ++os1) {
          auto obf1_first = oshell2bf[os1];
          auto on1 = oshells[os1].size();

          for (unsigned int os2 = os1; os2 != onshells; ++os2) {
            if ((s1234++) % nthreads != thread_id)
              continue;

            auto obf2_first = oshell2bf[os2];
            auto on2 = oshells[os2].size();

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx,
                            0>(shells[s1], shells[s2], oshells[os1],
                               oshells[os2]);

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif

            if (buf[0] == nullptr)
              continue; // all screened out

            for (unsigned int f1 = 0, f1212 = 0; f1 != n1; ++f1) {
              const auto bf1 = f1 + bf1_first;

              for (unsigned int f2 = 0; f2 != n2; ++f2) {
                const auto bf2 = f2 + bf2_first;

                for (unsigned int of1 = 0; of1 != on1; ++of1) {
                  const auto obf1 = of1 + obf1_first;

                  for (unsigned int of2 = 0; of2 != on2; ++of2, ++f1212) {
                    const auto obf2 = of2 + obf2_first;

                    if (bf1 <= bf2 && obf1 <= obf2) {
                      if (std::abs(buf[0][f1212]) < 1.0e-10)
                        continue;

                      buffer.p[counter] = bf1 + 1;
                      buffer.q[counter] = bf2 + 1;
                      buffer.r[counter] = obf1 + 1;
                      buffer.s[counter] = obf2 + 1;
                      buffer.val[counter] =
                          buf[0][f1212] * get_norma(bf1) * get_norma(bf2) *
                          other.get_norma(obf1) * other.get_norma(obf2);

                      ++num_ints_computed;
                      ++counter;

                      if (counter == STACK_SIZE) {
                        write_buffer(buffer, outfile);
                        counter = 0;
                      }
                    }
                  }
                }
              }
            } // end cartesian loop
          }
        }
      }
    } // end shells loop

    buffer.p[counter] = -1;
    write_buffer(buffer, outfile);
    outfile.close();

  }; // end of lambda

  napmo::parallel_do(lambda);

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto &t : timers) {
    time_for_ints += t.read(0);
  }
  std::cout << " Time for inter-species integrals = " << time_for_ints
            << std::endl;
  for (int t = 0; t != nthreads; ++t)
    engines[t].print_timers();
#endif

  // std::cout << " Number of unique integrals for species: " << sID << " / "
  //           << other.get_sID() << " = " << num_ints_computed << std::endl;
}

Matrix LibintInterface::compute_coupling_direct(LibintInterface &other,
                                                const Matrix &D,
                                                const bool permuted,
                                                double precision) {
  const auto n = permuted ? other.get_nbasis() : get_nbasis();

  const auto nshells = shells.size();
  const auto oshells = other.get_shells();
  const auto onshells = oshells.size();

  using napmo::nthreads;

  std::vector<Matrix> B(nthreads, Matrix::Zero(n, n));

  auto fock_precision = precision;

  // engine precision controls primitive truncation, assume worst-casescenario
  // (all primitive combinations add up constructively)
  auto nprim_max = std::max(max_nprim, other.max_nprim);
  auto l_max = std::max(max_l, other.max_l);

  auto max_nprim4 = nprim_max * nprim_max * nprim_max * nprim_max;

  auto engine_precision =
      std::min(fock_precision, std::numeric_limits<double>::epsilon()) /
      max_nprim4;

  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  engines[0] = Engine(libint2::Operator::coulomb, nprim_max, l_max, 0);
  engines[0].set_precision(engine_precision); // shellset-dependentprecision
                                              // control will likely break
                                              // positive definiteness
                                              // stick with this simple recipe

  // std::cout << "compute_coupling_direct:precision = " << precision <<
  // std::endl;
  // std::cout << "Engine::precision = " << engines[0].precision() <<
  // std::endl;

  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::atomic<size_t> num_ints_computed{0};

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = map_shell_to_basis_function();
  auto oshell2bf = other.map_shell_to_basis_function();

  auto lambda = [&](unsigned int thread_id) {

    auto &engine = engines[thread_id];
    const auto &buf = engines[thread_id].results();
    auto &b = B[thread_id];

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto &timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    // loop over permutationally-unique set of shells
    for (unsigned int s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = shells[s1].size();   // number of basis functions in this shell

      for (unsigned int s2 = s1; s2 != nshells; ++s2) {
        auto bf2_first = shell2bf[s2];
        auto n2 = shells[s2].size();

        for (unsigned int os1 = 0l; os1 != onshells; ++os1) {
          auto obf1_first = oshell2bf[os1];
          auto on1 = oshells[os1].size();

          for (unsigned int os2 = os1; os2 != onshells; ++os2) {
            if ((s1234++) % nthreads != thread_id)
              continue;

            auto obf2_first = oshell2bf[os2];
            auto on2 = oshells[os2].size();

// printf("(%d %d | %d %d) \n", s1, s2, os1, os2);

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx,
                            0>(shells[s1], shells[s2], oshells[os1],
                               oshells[os2]);

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif

            if (buf[0] == nullptr)
              continue; // all screened out

            for (unsigned int f1 = 0, f1212 = 0; f1 != n1; ++f1) {
              const auto bf1 = f1 + bf1_first;

              for (unsigned int f2 = 0; f2 != n2; ++f2) {
                const auto bf2 = f2 + bf2_first;

                for (unsigned int of1 = 0; of1 != on1; ++of1) {
                  const auto obf1 = of1 + obf1_first;

                  for (unsigned int of2 = 0; of2 != on2; ++of2, ++f1212) {
                    const auto obf2 = of2 + obf2_first;

                    if (bf1 <= bf2 && obf1 <= obf2) {
                      if (std::abs(buf[0][f1212]) < 1.0e-10)
                        continue;

                      ++num_ints_computed;

                      // printf("(%ld %ld | %ld %ld) %lf \n", bf1, bf2, obf1,
                      // obf2,
                      //        buf[0][f1212] * get_norma(bf1) * get_norma(bf2)
                      //        *
                      //            other.get_norma(obf1) *
                      //            other.get_norma(obf2));

                      if (permuted) {
                        auto coulomb = D(bf1, bf2) * buf[0][f1212] *
                                       get_norma(bf1) * get_norma(bf2) *
                                       other.get_norma(obf1) *
                                       other.get_norma(obf2);

                        b(obf2, obf1) += coulomb;
                        if (bf1 != bf2)
                          b(obf2, obf1) += coulomb;
                      } else {
                        auto coulomb = D(obf1, obf2) * buf[0][f1212] *
                                       get_norma(bf1) * get_norma(bf2) *
                                       other.get_norma(obf1) *
                                       other.get_norma(obf2);

                        // printf("%lf\n", coulomb);
                        b(bf2, bf1) += coulomb;
                        if (obf1 != obf2)
                          b(bf2, bf1) += coulomb;
                      }
                    }
                  }
                }
              }
            } // end cartesian loop
          }
        }
      }
    } // end shells loop

  }; // end of lambda

  napmo::parallel_do(lambda);

  // accumulate contributions from all threads
  for (size_t i = 1; i != nthreads; ++i) {
    B[0] += B[i];
  }

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto &t : timers) {
    time_for_ints += t.read(0);
  }
  std::cout << " Time for inter-species integrals = " << time_for_ints
            << std::endl;
  for (int t = 0; t != nthreads; ++t)
    engines[t].print_timers();
#endif

  // std::cout << " Number of unique integrals for species: " << sID << " / " << other.get_sID() << " = " << num_ints_computed << std::endl;

  return B[0];
}

/*
Auxiliary functions
*/
shellpair_list_t compute_shellpair_list(const std::vector<libint2::Shell> &bs1,
                                        const std::vector<libint2::Shell> &_bs2,
                                        const double threshold) {
  const std::vector<libint2::Shell> &bs2 = (_bs2.empty() ? bs1 : _bs2);
  const auto nsh1 = bs1.size();
  const auto nsh2 = bs2.size();
  const auto bs1_equiv_bs2 = (&bs1 == &bs2);

  using napmo::nthreads;

  // construct the 2-electron repulsion integrals engine
  using libint2::Engine;
  libint2::BasisSet obs_aux;
  std::vector<Engine> engines;
  engines.reserve(nthreads);
  engines.emplace_back(libint2::Operator::overlap,
                       std::max(obs_aux.max_nprim(bs1), obs_aux.max_nprim(bs2)),
                       std::max(obs_aux.max_l(bs1), obs_aux.max_l(bs2)), 0);

  for (size_t i = 1; i != nthreads; ++i) {
    engines.push_back(engines[0]);
  }

  // std::cout << "computing non-negligible shell-pair list ... ";

  libint2::Timers<1> timer;
  timer.set_now_overhead(25);
  timer.start(0);

  shellpair_list_t result;

  std::mutex mx;

  auto compute = [&](unsigned int thread_id) {

    auto &engine = engines[thread_id];
    const auto &buf = engines[thread_id].results();

    // loop over permutationally-unique set of shells
    for (unsigned int s1 = 0l, s12 = 0l; s1 != nsh1; ++s1) {
      mx.lock();
      if (result.find(s1) == result.end())
        result.insert(std::make_pair(s1, std::vector<size_t>()));
      mx.unlock();

      auto n1 = bs1[s1].size(); // number of basis functions in this shell

      auto s2_max = bs1_equiv_bs2 ? s1 : nsh2 - 1;
      for (unsigned int s2 = 0; s2 <= s2_max; ++s2, ++s12) {
        if (s12 % nthreads != thread_id)
          continue;

        auto on_same_center = (bs1[s1].O == bs2[s2].O);
        bool significant = on_same_center;
        if (not on_same_center) {
          auto n2 = bs2[s2].size();
          engine.compute(bs1[s1], bs2[s2]);
          Eigen::Map<const Matrix> buf_mat(buf[0], n1, n2);
          auto norm = buf_mat.norm();
          significant = (norm >= threshold);
        }

        if (significant) {
          mx.lock();
          result[s1].emplace_back(s2);
          mx.unlock();
        }
      }
    }
  }; // end of compute

  napmo::parallel_do(compute);

  // resort shell list in increasing order, i.e. result[s][s1] < result[s][s2]
  // if s1 < s2
  auto sort = [&](unsigned int thread_id) {
    for (unsigned int s1 = 0l; s1 != nsh1; ++s1) {
      if (s1 % nthreads == thread_id) {
        auto &list = result[s1];
        std::sort(list.begin(), list.end());
      }
    }
  }; // end of sort

  napmo::parallel_do(sort);

  timer.stop(0);
  // std::cout << "done (" << timer.read(0) << " s)" << std::endl;

  return result;
}

template <libint2::Operator Kernel>
Matrix compute_schwartz_ints(
    const std::vector<libint2::Shell> &bs1,
    const std::vector<libint2::Shell> &_bs2, bool use_2norm,
    typename libint2::operator_traits<Kernel>::oper_params_type params) {
  const std::vector<libint2::Shell> &bs2 = (_bs2.empty() ? bs1 : _bs2);
  const auto nsh1 = bs1.size();
  const auto nsh2 = bs2.size();
  const auto bs1_equiv_bs2 = (&bs1 == &bs2);

  Matrix K = Matrix::Zero(nsh1, nsh2);

  // construct the 2-electron repulsion integrals engine
  using libint2::Engine;
  using napmo::nthreads;
  std::vector<Engine> engines(nthreads);

  // !!! very important: cannot screen primitives in Schwartz computation !!!
  libint2::BasisSet obs_aux;
  auto epsilon = 0.;
  engines[0] = Engine(
      Kernel, std::max(obs_aux.max_nprim(bs1), obs_aux.max_nprim(bs2)),
      std::max(obs_aux.max_l(bs1), obs_aux.max_l(bs2)), 0, epsilon, params);
  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  // std::cout << "computing Schwartz bound prerequisites (kernel=" <<
  // (int)Kernel
  //           << ") ... ";

  libint2::Timers<1> timer;
  timer.set_now_overhead(25);
  timer.start(0);

  auto lambda = [&](unsigned int thread_id) {
    const auto &buf = engines[thread_id].results();

    // loop over permutationally-unique set of shells
    for (unsigned int s1 = 0l, s12 = 0l; s1 != nsh1; ++s1) {
      auto n1 = bs1[s1].size(); // number of basis functions in this shell

      auto s2_max = bs1_equiv_bs2 ? s1 : nsh2 - 1;
      for (unsigned int s2 = 0; s2 <= s2_max; ++s2, ++s12) {
        if (s12 % nthreads != thread_id)
          continue;

        auto n2 = bs2[s2].size();
        auto n12 = n1 * n2;

        engines[thread_id].compute(bs1[s1], bs2[s2], bs1[s1], bs2[s2]);

        // the diagonal elements are the Schwartz ints ... use Map.diagonal()
        Eigen::Map<const Matrix> buf_mat(buf[0], n12, n12);
        auto norm2 = use_2norm ? buf_mat.diagonal().norm()
                               : buf_mat.diagonal().lpNorm<Eigen::Infinity>();
        K(s1, s2) = std::sqrt(norm2);
        if (bs1_equiv_bs2)
          K(s2, s1) = K(s1, s2);
      }
    }
  }; // thread lambda

  napmo::parallel_do(lambda);

  timer.stop(0);
  // std::cout << "done (" << timer.read(0) << " s)" << std::endl;

  return K;
}

__inline__ void write_buffer(const QuartetBuffer &buffer,
                             std::ofstream &outfile) {
  outfile.write(reinterpret_cast<const char *>(&buffer.p[0]),
                STACK_SIZE * sizeof(int));
  outfile.write(reinterpret_cast<const char *>(&buffer.q[0]),
                STACK_SIZE * sizeof(int));
  outfile.write(reinterpret_cast<const char *>(&buffer.r[0]),
                STACK_SIZE * sizeof(int));
  outfile.write(reinterpret_cast<const char *>(&buffer.s[0]),
                STACK_SIZE * sizeof(int));
  outfile.write(reinterpret_cast<const char *>(&buffer.val[0]),
                STACK_SIZE * sizeof(double));
}

/*
Python functions
*/

LibintInterface *LibintInterface_new(int id) { return new LibintInterface(id); }

void LibintInterface_del(LibintInterface *lint) { lint->~LibintInterface(); }

void LibintInterface_add_pointcharges(LibintInterface *lint, const int z,
                                      const double *center) {
  lint->add_pointcharges(z, center);
}

void LibintInterface_add_basis(LibintInterface *lint, BasisSet *basis) {
  lint->add_basis(basis);
}

void LibintInterface_compute_1body_ints(LibintInterface *lint,
                                        int integral_kind, double *result) {

  // Default case
  libint2::Operator obtype = libint2::Operator::overlap;

  switch (integral_kind) {
  case 1: // Overlap
    obtype = libint2::Operator::overlap;
    break;
  case 2: // Kinetic
    obtype = libint2::Operator::kinetic;
    break;
  case 3: // Nuclear
    obtype = libint2::Operator::nuclear;
    break;
  }

  auto tmp = lint->compute_1body_ints(obtype);

  for (int i = 0; i < tmp.size(); ++i) {
    result[i] = tmp.array()(i);
  }
}

int LibintInterface_get_nbasis(LibintInterface *lint) {
  return lint->get_nbasis();
}

void LibintInterface_init_2body_ints(LibintInterface *lint) {
  lint->init_2body_ints();
}

std::vector<QuartetBuffer> *
LibintInterface_compute_2body_ints(LibintInterface *lint, double *dens) {

  auto n = lint->get_nbasis();
  Matrix D(n, n);

  for (int i = 0; i < D.size(); ++i) {
    D.array()(i) = dens[i];
  }

  auto shells = lint->get_shells();
  auto K = compute_schwartz_ints<>(shells);
  return lint->compute_2body_ints(D, K);
}

void LibintInterface_compute_2body_direct(LibintInterface *lint, double *dens,
                                          double *result) {
  auto n = lint->get_nbasis();
  Matrix D(n, n);

  for (int i = 0; i < D.size(); ++i) {
    D.array()(i) = dens[i];
  }

  auto shells = lint->get_shells();
  auto K = compute_schwartz_ints<>(shells);
  auto G = lint->compute_2body_direct(D, K);

  for (int i = 0; i < G.size(); ++i) {
    result[i] = G.array()(i);
  }
}

void LibintInterface_compute_2body_disk(LibintInterface *lint,
                                        const char *filename, double *dens) {
  auto n = lint->get_nbasis();
  Matrix D(n, n);

  for (int i = 0; i < D.size(); ++i) {
    D.array()(i) = dens[i];
  }
  auto shells = lint->get_shells();
  auto K = compute_schwartz_ints<>(shells);

  lint->compute_2body_disk(filename, D, K);
}

void LibintInterface_compute_coupling_direct(LibintInterface *lint,
                                             LibintInterface *olint,
                                             double *dens, double *result) {
  auto sid = lint->get_sID();
  auto osid = olint->get_sID();

  auto permuted = sid > osid;

  auto n = olint->get_nbasis();

  Matrix D(n, n);

  for (int i = 0; i < D.size(); ++i) {
    D.array()(i) = dens[i];
  }
  using namespace std;

  // cout<<D<<endl;

  auto B = permuted ? olint->compute_coupling_direct(*lint, D, permuted)
                    : lint->compute_coupling_direct(*olint, D, permuted);

  for (int i = 0; i < B.size(); ++i) {
    result[i] = B.array()(i);
  }
}

void LibintInterface_compute_coupling_disk(LibintInterface *lint,
                                           LibintInterface *olint,
                                           const char *filename) {
  lint->compute_coupling_disk(*olint, filename);
}

libint2::DIIS<Matrix> *LibintInterface_diis_new(int iter) {
  return new libint2::DIIS<Matrix>(iter);
}

void LibintInterface_diis(libint2::DIIS<Matrix> *diis, WaveFunction *psi) {
  int nbasis = psi->nbasis;

  MMap S(psi->S, nbasis, nbasis);
  MMap D(psi->D, nbasis, nbasis);
  MMap F(psi->F, nbasis, nbasis);

  // compute SCF error
  Matrix FD_comm = F * D * S - S * D * F;
  Matrix F_diis = F;

  diis->extrapolate(F_diis, FD_comm);

  for (int i = 0; i < F_diis.size(); ++i) {
    psi->F[i] = F_diis.array()(i);
  }
}
