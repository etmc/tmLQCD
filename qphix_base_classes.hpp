// Copyright © 2017 Martin Ueding <dev@martin-ueding.de>
// Licensed unter the [BSD-3-Clause](https://opensource.org/licenses/BSD-3-Clause).

/**
  \file Additions to QPhiX that are only needed for tmLQCD.

  In the original QPhiX, there are only Wilson fermions and Wilson clover
  fermions. The Dslash operators have a different call signature (the latter
  requiring a clover term), so there is no common base class. With the addition
  of Wilson twisted mass (Mario) and Wilson twisted clover (Peter), there are
  now two instances of the Dslash that have the same signature. In order to
  write a more general even-odd source preparation and solution reconstruction
  code, a common base class for non-clover and clover is desired. In order to
  leave the QPhiX code untouched (for now), this code lives here in tmLQCD.
  */

#pragma once

#include <qphix/blas_new_c.h>
#include <qphix/clover_dslash_def.h>
#include <qphix/dslash_def.h>
#include <qphix/tm_clov_dslash_def.h>
#include <qphix/tm_dslash_def.h>

#include <cassert>

namespace tmlqcd {

namespace {
size_t constexpr re = 0;
size_t constexpr im = 1;
int const n_blas_simt = 1;
}

/**
  Complex muliplication accumulate.

  Computes \f$ (r + \mathrm i i) += (a + \mathrm i b) * (c + \mathrm i d) \f$.
  */
template <typename FT>
void cplx_mul_acc(FT &r_out, FT &i_out, FT const &a, FT const &b, FT const &c, FT const &d) {
  r_out += a * c - b * d;
  i_out += a * d + b * c;
}

template <typename FT, int veclen, int soalen, bool compress12>
void clover_product(
    typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock *const out,
    typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock const *const in,
    typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::CloverBlock *local_clover,
    ::QPhiX::Geometry<FT, veclen, soalen, compress12> &geom) {
  ::QPhiX::zeroSpinor<FT, veclen, soalen, compress12>(out, geom, n_blas_simt);

  // Iterate through all the block.
  auto const num_blocks = geom.get_num_blocks();
  for (auto block = 0u; block < num_blocks; ++block) {
    // The clover term is block-diagonal in spin. Therefore we need
    // to iterate over the two blocks of spin.
    for (auto s_block : {0, 1}) {
      // Extract the diagonal and triangular parts.
      auto const &diag_in = s_block == 0 ? local_clover[block].diag1 : local_clover[block].diag2;
      auto const &off_diag_in = s_block == 0 ? local_clover[block].off_diag1 : local_clover[block].off_diag1;
      // Input two-spinor component.
      for (auto two_s_in : {0, 1}) {
        // Reconstruct four spinor index.
        auto const four_s_in = 2 * s_block + two_s_in;
        // Output two-spinor component.
        for (auto two_s_out : {0, 1}) {
          // Reconstruct four spinor index.
          auto const four_s_out = 2 * s_block + two_s_out;
          // Input color.
          for (auto c_in : {0, 1, 2}) {
            // Spin-color index (0, ..., 5).
            auto const sc_in = 3 * two_s_in + c_in;
            // Output color.
            for (auto c_out : {0, 1, 2}) {
              // Spin-color index (0, ..., 5).
              auto const sc_out = 3 * two_s_out + c_out;
              // SIMD vector.
              for (auto v = 0; v < veclen; ++v) {
                  if (sc_out == sc_in) {
                    cplx_mul_acc(out[block][c_out][four_s_out][re][v],
                                 out[block][c_out][four_s_out][im][v],
                                 diag_in[sc_in][v],
                                 0.0,
                                 in[block][c_in][four_s_in][re][v],
                                 in[block][c_in][four_s_in][im][v]);
                  }
                  else if (sc_out < sc_in) {
                    auto const idx15 = sc_in * (sc_in - 1) / 2 + sc_out;
                    cplx_mul_acc(out[block][c_out][four_s_out][re][v],
                                 out[block][c_out][four_s_out][im][v],
                                 off_diag_in[idx15][re][v],
                                 off_diag_in[idx15][im][v],
                                 in[block][c_in][four_s_in][re][v],
                                 in[block][c_in][four_s_in][im][v]);
                  }
              }
            }
          }
        }
      }
    }
  }
}

template <typename FT, int veclen, int soalen, bool compress12>
void full_clover_product(
    typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock *const out,
    typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock const *const in,
    typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FullCloverBlock *local_clover,
    ::QPhiX::Geometry<FT, veclen, soalen, compress12> &geom) {
  ::QPhiX::zeroSpinor<FT, veclen, soalen, compress12>(out, geom, n_blas_simt);

  // Iterate through all the block.
  auto const num_blocks = geom.get_num_blocks();
  for (auto block = 0u; block < num_blocks; ++block) {
    // The clover term is block-diagonal in spin. Therefore we need
    // to iterate over the two blocks of spin.
    for (auto s_block : {0, 1}) {
      // Extract the spin block as a handy alias.
      auto const &block_in = s_block == 0 ? local_clover[block].block1 : local_clover[block].block2;
      // Input two-spinor component.
      for (auto two_s_in : {0, 1}) {
        // Reconstruct four spinor index.
        auto const four_s_in = 2 * s_block + two_s_in;
        // Output two-spinor component.
        for (auto two_s_out : {0, 1}) {
          // Reconstruct four spinor index.
          auto const four_s_out = 2 * s_block + two_s_out;
          // Input color.
          for (auto c_in : {0, 1, 2}) {
            // Spin-color index (0, ..., 5).
            auto const sc_in = 3 * two_s_in + c_in;
            // Output color.
            for (auto c_out : {0, 1, 2}) {
              // Spin-color index (0, ..., 5).
              auto const sc_out = 3 * two_s_out + c_out;
              // SIMD vector.
              for (auto v = 0; v < veclen; ++v) {
                cplx_mul_acc(out[block][c_out][four_s_out][re][v],
                             out[block][c_out][four_s_out][im][v],
                             block_in[sc_out][sc_in][re][v],
                             block_in[sc_out][sc_in][im][v],
                             in[block][c_in][four_s_in][re][v],
                             in[block][c_in][four_s_in][im][v]);
              }
            }
          }
        }
      }
    }
  }
}

template <typename FT, int veclen, int soalen, bool compress12>
class FourSpinorCBWrapper {
 public:
  typedef
      typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock FourSpinorBlock;

  FourSpinorCBWrapper(::QPhiX::Geometry<FT, veclen, soalen, compress12> &geom_)
      : geom(geom_), spinor(geom.allocCBFourSpinor()) {}

  ~FourSpinorCBWrapper() { geom.free(spinor); }

  FourSpinorBlock *data() const { return spinor; }
  size_t num_blocks() const { return geom.get_num_blocks(); }

 private:
  ::QPhiX::Geometry<FT, veclen, soalen, compress12> &geom;
  FourSpinorBlock *const spinor;
};

template <typename FT, int veclen, int soalen, bool compress12>
class Dslash {
 public:
  typedef typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock Spinor;
  typedef typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::SU3MatrixBlock SU3MatrixBlock;

  Dslash(typename ::QPhiX::Geometry<FT, veclen, soalen, compress12> &geom_) : tmp(geom_) {}

  virtual ~Dslash() {}

  /**
    Computes \f$ \psi_\mathrm o = A_\mathrm{oo} \chi_\mathrm o \f$.

    The actual definition of the matrix \f$ A_\mathrm{oo} \f$ is
    implementation dependent and can be the mass factor \f$ \alpha = 4 + m
    \f$ for plain Wilson or something more complicated for twisted mass.

    \param[out] out Output spinor \f$ \psi \f$.
    \param[in] in Input spinor \f$ \chi \f$.
    */
  virtual void A_chi(Spinor *const out, Spinor const *const in, int const isign) = 0;

  /**
    Computes \f$ \psi_\mathrm e = A_\mathrm{ee}^{-1} \chi_\mathrm e \f$.

    \param[out] out Output spinor \f$ \psi \f$.
    \param[in] in Input spinor \f$ \chi \f$.
    */
  virtual void A_inv_chi(Spinor *const out, Spinor const *const in, int const isign) = 0;

  /**
    Forwarder for the `dslash`.

    \todo Make this member function `const`. For this the member function in
    QPhiX that is called internally must be marked `const` as well.
    */
  virtual void dslash(Spinor *const res,
                      const Spinor *const psi,
                      const SU3MatrixBlock *const u,
                      int const isign,
                      int const cb) = 0;

  /**
    Plain Wilson Dslash.

    In the source preparation and solution reconstruction, there is a step
    where the plain Wilson Dslash has to be used. This could have been solved
    by creating an instance of the Wilson Dslash in the solver. The solution
    chosen is that the actual Dslash is used and the result is then corrected
    by multiplying with the clover term again. This way a plain Wilson Dslash
    is provided.

    This implementation is meant for all Dslash types which are _not_ the plain
    Wilson Dslash. That particular one has to override this method with does
    _no_ correction.
    */
  virtual void wilson_dslash(Spinor *const res,
                             const Spinor *const psi,
                             const SU3MatrixBlock *const u,
                             int const isign,
                             int const cb) {
    assert(isign == 1 && "isign = -1 is not implemented");
    dslash(tmp.data(), psi, u, isign, cb);
    A_chi(res, tmp.data(), isign);
  }

  /**
    Forwarder for the `achimbdpsi`.

    \todo Make this member function `const`. For this the member function in QPhiX that is called
    internally must be marked `const` as well.
    */
  virtual void achimbdpsi(Spinor *const res,
                          const Spinor *const psi,
                          const Spinor *const chi,
                          const SU3MatrixBlock *const u,
                          double const alpha,
                          double const beta,
                          int const isign,
                          int const cb) = 0;

 protected:
  FourSpinorCBWrapper<FT, veclen, soalen, compress12> tmp;
};

template <typename FT, int veclen, int soalen, bool compress12>
class WilsonDslash : public Dslash<FT, veclen, soalen, compress12> {
 public:
  typedef typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock Spinor;
  typedef typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::SU3MatrixBlock SU3MatrixBlock;

  WilsonDslash(::QPhiX::Geometry<FT, veclen, soalen, compress12> *geom_,
               double const t_boundary_,
               double const aniso_coeff_S_,
               double const aniso_coeff_T_,
               double const mass_)
      : Dslash<FT, veclen, soalen, compress12>(*geom_),
        upstream_dslash(geom_, t_boundary_, aniso_coeff_S_, aniso_coeff_T_),
        mass_factor_alpha(4.0 + mass_),
        mass_factor_beta(1.0 / (4.0 * mass_factor_alpha)) {}

  void A_chi(Spinor *const out, Spinor const *const in, int const isign) override {
    ::QPhiX::axy(mass_factor_beta, in, out, upstream_dslash.getGeometry(), n_blas_simt);
  }

  void A_inv_chi(Spinor *const out, Spinor const *const in, int const isign) override {
    ::QPhiX::axy(1.0 / mass_factor_beta, in, out, upstream_dslash.getGeometry(), n_blas_simt);
  }

  void dslash(Spinor *const res,
              const Spinor *const psi,
              const SU3MatrixBlock *const u,
              int const isign,
              int const cb) override {
    upstream_dslash.dslash(res, psi, u, isign, cb);
  }

  void wilson_dslash(Spinor *const res,
              const Spinor *const psi,
              const SU3MatrixBlock *const u,
              int const isign,
              int const cb) override {
    dslash(res, psi, u, isign, cb);
  }

  void achimbdpsi(Spinor *const res,
                  const Spinor *const psi,
                  const Spinor *const chi,
                  const SU3MatrixBlock *const u,
                  double const alpha,
                  double const beta,
                  int const isign,
                  int const cb) override {
    upstream_dslash.dslashAChiMinusBDPsi(res, psi, chi, u, alpha, beta, isign, cb);
  }

 private:
  ::QPhiX::Dslash<FT, veclen, soalen, compress12> upstream_dslash;

  double const mass_factor_alpha;
  double const mass_factor_beta;
};

template <typename FT, int veclen, int soalen, bool compress12>
class WilsonTMDslash : public Dslash<FT, veclen, soalen, compress12> {
 public:
  typedef typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock Spinor;
  typedef typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::SU3MatrixBlock SU3MatrixBlock;

  WilsonTMDslash(::QPhiX::Geometry<FT, veclen, soalen, compress12> *geom_,
                 double const t_boundary_,
                 double const aniso_coeff_S_,
                 double const aniso_coeff_T_,
                 double const mass_,
                 double const twisted_mass_)
      : Dslash<FT, veclen, soalen, compress12>(*geom_),
        upstream_dslash(geom_, t_boundary_, aniso_coeff_S_, aniso_coeff_T_),
        mass_factor_alpha(4.0 + mass_),
        mass_factor_beta(0.25),
        derived_mu(twisted_mass_ / mass_factor_alpha),
        derived_mu_inv(mass_factor_alpha /
                       (mass_factor_alpha * mass_factor_alpha + twisted_mass_ * twisted_mass_)) {}

  void A_chi(Spinor *const out, Spinor const *const in, int const isign) override {
      helper_A_chi(out, in, derived_mu * isign, mass_factor_alpha);
  }


  void A_inv_chi(Spinor *const out, Spinor const *const in, int const isign) override {
      helper_A_chi(out, in, derived_mu * isign, derived_mu_inv);
  }

  void dslash(Spinor *const res,
              const Spinor *const psi,
              const SU3MatrixBlock *const u,
              int const isign,
              int const cb) override {
    upstream_dslash.tmdslash(res, psi, u, derived_mu, derived_mu_inv, isign, cb);
  }

  void achimbdpsi(Spinor *const res,
                  const Spinor *const psi,
                  const Spinor *const chi,
                  const SU3MatrixBlock *const u,
                  double const alpha,
                  double const beta,
                  int const isign,
                  int const cb) override {
    upstream_dslash.tmdslashAChiMinusBDPsi(
        res, psi, chi, u, derived_mu, derived_mu_inv, isign, cb);
  }

 private:
  void helper_A_chi(Spinor *const out,
                    Spinor const *const in,
                    double const factor_a,
                    double const factor_b);

  ::QPhiX::TMDslash<FT, veclen, soalen, compress12> upstream_dslash;

  double const mass_factor_alpha;
  double const mass_factor_beta;
  double const derived_mu;
  double const derived_mu_inv;
};

template <typename FT, int veclen, int soalen, bool compress12>
void WilsonTMDslash<FT, veclen, soalen, compress12>::helper_A_chi(Spinor *const out,
                                                                  Spinor const *const in,
                                                                  double const factor_a,
                                                                  double const factor_b) {
  size_t const num_blocks = upstream_dslash.getGeometry().get_num_blocks();
  for (size_t block = 0u; block < num_blocks; ++block) {
    for (int color = 0; color < 3; ++color) {
      for (int spin_block = 0; spin_block < 2; ++spin_block) {
        // Implement the $\gamma_5$ structure.
        auto const signed_factor_a = factor_a * (spin_block == 0 ? 1.0 : -1.0);

        for (int half_spin = 0; half_spin < 2; ++half_spin) {
          auto const four_spin = 2 * spin_block + half_spin;
          for (int v = 0; v < soalen; ++v) {
            auto &out_bcs = out[block][color][four_spin];
            auto const &in_bcs = in[block][color][four_spin];

            out_bcs[re][v] = factor_b * (in_bcs[re][v] + signed_factor_a * in_bcs[im][v]);
            out_bcs[im][v] = factor_b * (in_bcs[im][v] - signed_factor_a * in_bcs[re][v]);
          }
        }
      }
    }
  }
};

template <typename FT, int veclen, int soalen, bool compress12>
class WilsonClovDslash : public Dslash<FT, veclen, soalen, compress12> {
 public:
  typedef typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock Spinor;
  typedef typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::SU3MatrixBlock SU3MatrixBlock;
  typedef typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::CloverBlock CloverBlock;

  WilsonClovDslash(::QPhiX::Geometry<FT, veclen, soalen, compress12> *geom_,
                   double const t_boundary_,
                   double const aniso_coeff_S_,
                   double const aniso_coeff_T_,
                   double const mass_,
                   CloverBlock *const clover_,
                   CloverBlock *const inv_clover_)
      : Dslash<FT, veclen, soalen, compress12>(*geom_),
        upstream_dslash(geom_, t_boundary_, aniso_coeff_S_, aniso_coeff_T_),
        mass_factor_alpha(4.0 + mass_),
        mass_factor_beta(1.0 / (4.0 * mass_factor_alpha)),
        clover(clover_),
        inv_clover(inv_clover_) {}

  void A_chi(Spinor *const out, Spinor const *const in, int const isign_ignored) override {
    clover_product(out, in, clover, upstream_dslash.getGeometry());
  }

  void A_inv_chi(Spinor *const out, Spinor const *const in, int const isign_ignored) override {
    clover_product(out, in, inv_clover, upstream_dslash.getGeometry());
  }

  void dslash(Spinor *const res,
              const Spinor *const psi,
              const SU3MatrixBlock *const u,
              int const isign,
              int const cb) override {
    upstream_dslash.dslash(res, psi, u, inv_clover, isign, cb);
  }

  void achimbdpsi(Spinor *const res,
                  const Spinor *const psi,
                  const Spinor *const chi,
                  const SU3MatrixBlock *const u,
                  double const alpha,
                  double const beta,
                  int const isign,
                  int const cb) override {
    upstream_dslash.dslashAChiMinusBDPsi(res, psi, chi, u, clover, mass_factor_beta, isign, cb);
  }

 private:
  ::QPhiX::ClovDslash<FT, veclen, soalen, compress12> upstream_dslash;

  double const mass_factor_alpha;
  double const mass_factor_beta;

  CloverBlock *const clover;
  CloverBlock *const inv_clover;
};

template <typename FT, int veclen, int soalen, bool compress12>
class WilsonClovTMDslash : public Dslash<FT, veclen, soalen, compress12> {
 public:
  typedef typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock Spinor;
  typedef typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::SU3MatrixBlock SU3MatrixBlock;
  typedef
      typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FullCloverBlock FullCloverBlock;

  WilsonClovTMDslash(::QPhiX::Geometry<FT, veclen, soalen, compress12> *geom_,
                     double const t_boundary_,
                     double const aniso_coeff_S_,
                     double const aniso_coeff_T_,
                     double const mass_,
                     double const twisted_mass_,
                     FullCloverBlock *const clover_[2],
                     FullCloverBlock *const inv_clover_[2])
      : Dslash<FT, veclen, soalen, compress12>(*geom_),
        upstream_dslash(geom_, t_boundary_, aniso_coeff_S_, aniso_coeff_T_),
        mass_factor_alpha(4.0 + mass_),
        mass_factor_beta(0.25),
        derived_mu(twisted_mass_ / mass_factor_alpha),
        derived_mu_inv(mass_factor_alpha /
                       (mass_factor_alpha * mass_factor_alpha + twisted_mass_ * twisted_mass_)),
        clover(clover_),
        inv_clover(inv_clover_) {}

  void A_chi(Spinor *const out, Spinor const *const in, int const isign) override {
    full_clover_product(out, in, clover[isign], upstream_dslash.getGeometry());
  }


  void A_inv_chi(Spinor *const out, Spinor const *const in, int const isign) override {
    full_clover_product(out, in, inv_clover[isign], upstream_dslash.getGeometry());
  }

  void dslash(Spinor *const res,
              const Spinor *const psi,
              const SU3MatrixBlock *const u,
              int const isign,
              int const cb) override {
    upstream_dslash.dslash(res, psi, u, inv_clover, isign, cb);
  }

  void achimbdpsi(Spinor *const res,
                  const Spinor *const psi,
                  const Spinor *const chi,
                  const SU3MatrixBlock *const u,
                  double const alpha,
                  double const beta,
                  int const isign,
                  int const cb) override {
    upstream_dslash.dslashAChiMinusBDPsi(res, psi, chi, u, clover, mass_factor_beta, isign, cb);
  }

 private:
  ::QPhiX::TMClovDslash<FT, veclen, soalen, compress12> upstream_dslash;

  double const mass_factor_alpha;
  double const mass_factor_beta;
  double const derived_mu;
  double const derived_mu_inv;

  FullCloverBlock *const clover[2];
  FullCloverBlock *const inv_clover[2];
};
}
