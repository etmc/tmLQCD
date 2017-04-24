/***********************************************************************
 * Copyright (C) 2017 Martin Ueding <dev@martin-ueding.de>
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/

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
#include <qphix/dslash_def.h>
#include <qphix/tm_dslash_def.h>

namespace tmlqcd {

namespace {
size_t constexpr re = 0;
size_t constexpr im = 1;
}

template <typename FT, int veclen, int soalen, bool compress12>
class Dslash {
 public:
  typedef typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock Spinor;
  typedef typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::SU3MatrixBlock SU3MatrixBlock;

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

    \todo Make this member function `const`. For this the member function in QPhiX that is called
    internally must be marked `const` as well.
    */
  virtual void dslash(Spinor *const res,
                      const Spinor *const psi,
                      const SU3MatrixBlock *const u,
                      int const isign,
                      int const cb) = 0;

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
      : upstream_dslash(geom_, t_boundary_, aniso_coeff_S_, aniso_coeff_T_),
        mass_factor_alpha(4.0 + mass_),
        mass_factor_beta(1.0 / (4.0 * mass_factor_alpha)) {}

  void A_chi(Spinor *const out, Spinor const *const in, int const isign) override {
    int const n_blas_simt = 1;
    ::QPhiX::axy(mass_factor_alpha, in, out, upstream_dslash.getGeometry(), n_blas_simt);
  }

  void A_inv_chi(Spinor *const out, Spinor const *const in, int const isign) override {
    int const n_blas_simt = 1;
    ::QPhiX::axy(1.0 / mass_factor_alpha, in, out, upstream_dslash.getGeometry(), n_blas_simt);
  }

  void dslash(Spinor *const res,
              const Spinor *const psi,
              const SU3MatrixBlock *const u,
              int const isign,
              int const cb) override {
    upstream_dslash.dslash(res, psi, u, isign, cb);
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
      : upstream_dslash(geom_, t_boundary_, aniso_coeff_S_, aniso_coeff_T_),
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
}

