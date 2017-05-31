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
#include <qphix/geometry.h>
#include <qphix/tm_clov_dslash_def.h>
#include <qphix/tm_dslash_def.h>

#include <cassert>

namespace tmlqcd {

namespace {
size_t constexpr re = 0;
size_t constexpr im = 1;
int const n_blas_simt = 1;

// The even checkerboard is given by ( (x + y + z + t ) & 1 == 0 ) -> cb0 is even
int constexpr cb_even = 0;
int constexpr cb_odd = 1;
}

/**
  Complex multiplication accumulate.

  Computes \f$ (r + \mathrm i i) += (a + \mathrm i b) * (c + \mathrm i d) \f$.
  */
template <typename FT>
void cplx_mul_acc(FT &r_out, FT &i_out, FT const &a, FT const &b, FT const &c, FT const &d) {
  r_out += a * c - b * d;
  i_out += a * d + b * c;
}

/**
  Wrapper for the clover multiplication function.

  The `struct` is needed in order to allow for partial template specialization in the `Clover`
  parameter.

  \tparam Clover Type of clover block to use, must be a type from Geometry such that there exists a
  specialization for it.
  */
template <typename FT, int veclen, int soalen, bool compress12, typename Clover>
struct InnerCloverProduct {
  /**
  Multiplies the clover term for a single lattice size to a spinor.

  This function is intended to be used in a loop over all lattice sites. It is expected from the
  caller to have figured out all the correct indices. There are template specializations for the two
  different types of clover term that are used in QPhiX.

  \param[out] out Output spinor block. It is assumed to be zeroed properly, the function will just
  accumulate values into that output variable. Use \ref QPhiX::zeroSpinor for that.
  \param[in] in Input spinor block.
  \param[in] clover Single clover block that contains the lattice site of the spinor.
  \param[in] xi SIMD index for the arrays with length `soalen`, as in the spinors.
  \param[in] veclen_idx SIMD index for the arrays with length `veclen`, as in the clover term.
  */
  static void multiply(
      typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock &out,
      typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock const &in,
      Clover const &clover, int const xi, int const veclen_idx);
};

template <typename FT, int veclen, int soalen, bool compress12>
struct InnerCloverProduct<FT, veclen, soalen, compress12,
                          typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::CloverBlock> {
  static void multiply(
      typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock &spinor_out,
      typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock const &spinor_in,
      typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::CloverBlock const &clov_block,
      int const xi, int const veclen_idx) {
    // The clover term is block-diagonal in spin. Therefore we need
    // to iterate over the two blocks of spin.
    for (auto s_block : {0, 1}) {
      // Extract the diagonal and triangular parts.
      auto const &diag_in = s_block == 0 ? clov_block.diag1 : clov_block.diag2;
      auto const &off_diag_in = s_block == 0 ? clov_block.off_diag1 : clov_block.off_diag2;
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

              // See `qphix-codegen` file `dslash_common.cc`
              // function
              // `clover_term` for the index manipulations done
              // here.

              // Using separate loops over the actual indices is
              // probably
              // faster than the branching in the innermost loop.

              if (sc_out == sc_in) {
                cplx_mul_acc(spinor_out[c_out][four_s_out][re][xi],
                             spinor_out[c_out][four_s_out][im][xi], diag_in[sc_in][veclen_idx],
                             FT{0}, spinor_in[c_in][four_s_in][re][xi],
                             spinor_in[c_in][four_s_in][im][xi]);
              } else if (sc_out < sc_in) {
                auto const idx15 = sc_in * (sc_in - 1) / 2 + sc_out;
                cplx_mul_acc(
                    spinor_out[c_out][four_s_out][re][xi], spinor_out[c_out][four_s_out][im][xi],
                    off_diag_in[idx15][re][veclen_idx], -off_diag_in[idx15][im][veclen_idx],
                    spinor_in[c_in][four_s_in][re][xi], spinor_in[c_in][four_s_in][im][xi]);
              } else {
                auto const idx15 = sc_out * (sc_out - 1) / 2 + sc_in;
                cplx_mul_acc(
                    spinor_out[c_out][four_s_out][re][xi], spinor_out[c_out][four_s_out][im][xi],
                    off_diag_in[idx15][re][veclen_idx], off_diag_in[idx15][im][veclen_idx],
                    spinor_in[c_in][four_s_in][re][xi], spinor_in[c_in][four_s_in][im][xi]);
              }
            }
          }
        }
      }
    }
  }
};

template <typename FT, int veclen, int soalen, bool compress12>
struct InnerCloverProduct<
    FT, veclen, soalen, compress12,
    typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FullCloverBlock> {
  static void multiply(
      typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock &spinor_out,
      typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock const &spinor_in,
      typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FullCloverBlock const &clov_block,
      int const xi, int const veclen_idx) {
    // The clover term is block-diagonal in spin. Therefore we need
    // to iterate over the two blocks of spin.
    for (auto s_block : {0, 1}) {
      // Extract the diagonal and triangular parts.
      auto const &block_in = s_block == 0 ? clov_block.block1 : clov_block.block2;
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

              cplx_mul_acc(
                  spinor_out[c_out][four_s_out][re][xi], spinor_out[c_out][four_s_out][im][xi],
                  block_in[sc_out][sc_in][re][veclen_idx], block_in[sc_out][sc_in][im][veclen_idx],
                  spinor_in[c_in][four_s_in][re][xi], spinor_in[c_in][four_s_in][im][xi]);
            }
          }
        }
      }
    }
  }
};

/**
  Multiplies a checkerboarded QPhiX Clover term with a checkerboarded QPhiX spinor.

  Padding is taken care of. A test case for (a copy of) this function exists in QPhiX.

  If the preprocessor macro `PRINT_MAPPING` is defined, it will print out the mapping of `(x, y, z,
  t)` coordinates to block indices. Also it will check that each block is accessed the proper number
  of times, that is `soalen` for spinors and `veclen` for clover blocks.

  \param[out] out Output spinor
  \param[in] in Input spinor
  \param[in] clover Clover block
  \param[in] geom Geometry object holding the dimension of clover and spinor
  */
template <typename FT, int veclen, int soalen, bool compress12, typename Clover>
void clover_product(
    typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock *const out,
    typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock const *const in,
    Clover *clover, ::QPhiX::Geometry<FT, veclen, soalen, compress12> &geom) {
  ::QPhiX::zeroSpinor<FT, veclen, soalen, compress12>(out, geom, n_blas_simt);

#ifdef PRINT_MAPPING
  std::vector<int> spin_touches(geom.getPxyz() * geom.Nt(), 0);
  std::vector<int> clover_touches(geom.getPxyz() * geom.Nt() * soalen / veclen, 0);

  std::cout << std::setw(3) << "x" << std::setw(3) << "y" << std::setw(3) << "z" << std::setw(3)
            << "t"
            << ":" << std::setw(5) << "spin" << std::setw(5) << "clov"
            << "\n";
#endif

  // Iterate through all the block.
  for (int t = 0; t < geom.Nt(); ++t) {
    for (int z = 0; z < geom.Nz(); ++z) {
      for (int y = 0; y < geom.Ny(); ++y) {
        for (int x = 0; x < geom.Nxh(); ++x) {
          // First element in the current XY plane at desired Z and T.
          auto const xyBase = t * geom.getPxyz() + z * geom.getPxy();
          // Index of the SoA along the X direction.
          auto const xb = x / soalen;
          // Index within the SoA.
          auto const xi = x % soalen;
          // Global spin block index.
          auto const spin_block_idx = xb + geom.Nxh() / soalen * y + xyBase;
          // Global clover/gauge block index.
          auto const clov_block_idx =
              xb + (y / geom.nGY()) * geom.Nxh() / soalen + xyBase / geom.nGY();
          // Index of the SoA structure within the current tile.
          // auto const tile = (geom.Nxh() / soalen * y + xyBase) % geom.nGY();
          auto const tile = y % geom.nGY();
          // Vector index for clover/gauge. The SoA index only runs to
          // `soalen`, this index needs to run to `veclen`, that is across the
          // various SoA within the tile.
          auto const veclen_idx = soalen * tile + xi;

#ifdef PRINT_MAPPING
          ++spin_touches[spin_block_idx];
          ++clover_touches[clov_block_idx];

          std::cout << std::setw(3) << x << std::setw(3) << y << std::setw(3) << z << std::setw(3)
                    << t << ":" << std::setw(5) << spin_block_idx << std::setw(5) << clov_block_idx
                    << "\n";
#endif

          assert(xi + xb * soalen == x);

          // References to the objects at desired block.
          auto const &clov_block = clover[clov_block_idx];
          auto const &spinor_in = in[spin_block_idx];
          auto &spinor_out = out[spin_block_idx];

          InnerCloverProduct<FT, veclen, soalen, compress12, Clover>::multiply(
              spinor_out, spinor_in, clov_block, xi, veclen_idx);
        }
      }
    }
  }

#ifdef PRINT_MAPPING
  std::cout << std::flush;

  // Make sure that each block got touched the correct number of times.
  for (int i = 0; i != spin_touches.size(); ++i) {
    if (spin_touches[i] != soalen) {
      std::cout << "Spin missmatch: Block " << std::setw(4) << i << " accessed " << std::setw(4)
                << spin_touches[i] << " times instead of " << soalen << "\n";
    }
  }

  for (int i = 0; i != clover_touches.size(); ++i) {
    if (clover_touches[i] != veclen) {
      std::cout << "Clover missmatch: Block " << std::setw(4) << i << " accessed " << std::setw(4)
                << clover_touches[i] << " times instead of " << veclen << "\n";
    }
  }

  std::cout << std::flush;
#endif
}

template <typename FT, int veclen, int soalen, bool compress12>
void full_clover_product(
    typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock *const out,
    typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock const *const in,
    typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FullCloverBlock *local_clover,
    ::QPhiX::Geometry<FT, veclen, soalen, compress12> &geom) {
  ::QPhiX::zeroSpinor<FT, veclen, soalen, compress12>(out, geom, n_blas_simt);

  // Iterate through all the block.
  auto const num_blocks =  (geom.Nt() * geom.getPxyz() + geom.Nz() * geom.getPxy() ) / geom.nGY() + (geom.Ny() / geom.nGY()) * geom.nVecs();
  QPhiX::masterPrintf("num_blocks=%d\n",num_blocks);
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
#if 0
                QPhiX::masterPrintf("block_in[sc1:%d][sc2:%d][reim:%d][v:%d]: %lf\n", sc_out, sc_in, re, v ,block_in[sc_out][sc_in][re][v]);
#endif
                // FIXME this cannot be correct, spinors and gauge-fields have different blockings!
                cplx_mul_acc(out[block][c_out][four_s_out][re][v],
                             out[block][c_out][four_s_out][im][v], block_in[sc_out][sc_in][re][v],
                             block_in[sc_out][sc_in][im][v], in[block][c_in][four_s_in][re][v],
                             in[block][c_in][four_s_in][im][v]);
              }
            }
          }
        }
      }
    }
  }
}

/**
  Abstract base class for all single-flavor Dslash variants.

  There are four Dslash operators which are implemented in QPhiX:

  - Wilson
  - Wilson clover
  - Wilson twisted mass
  - Wilson clover with twisted mass

  Each of these has a the actual Dslash operation and a so-called “achimbdpsi” operation. These act
  on four-spinors given a gauge field. This base class provides a uniform interface to all four
  kinds.

  This code should eventually be migrated into the QPhiX repository. Currently these classes are
  mere delegators. In the QPhiX repository, the actual classes there should be used as concrete
  classes.
  */
template <typename FT, int veclen, int soalen, bool compress12>
class Dslash {
 public:
  typedef ::QPhiX::Geometry<FT, veclen, soalen, compress12> Geom;
  typedef typename Geom::FourSpinorBlock Spinor;
  typedef typename Geom::SU3MatrixBlock SU3MatrixBlock;

  explicit Dslash(Geom *geom, double const t_boundary_, double const aniso_coeff_S_,
                  double const aniso_coeff_T_, double const mass_)
      : geom(geom),
        t_boundary(t_boundary_),
        aniso_coeff_S(aniso_coeff_S_),
        aniso_coeff_T(aniso_coeff_T_),
        mass(mass_) {}

  /**
    Computes \f$ \psi_\mathrm o = A_\mathrm{oo} \chi_\mathrm o \f$.

    The actual definition of the matrix \f$ A_\mathrm{oo} \f$ is
    implementation dependent and can be the mass factor \f$ \alpha = 4 + m
    \f$ for plain Wilson or something more complicated for twisted mass.

    \param[out] out Output spinor \f$ \psi \f$.
    \param[in] in Input spinor \f$ \chi \f$.
    */
  virtual void A_chi(Spinor *const out, Spinor const *const in, int const isign, int const cb) = 0;

  /**
    Computes \f$ \psi_\mathrm e = A_\mathrm{ee}^{-1} \chi_\mathrm e \f$.

    \param[out] out Output spinor \f$ \psi \f$.
    \param[in] in Input spinor \f$ \chi \f$.
    */
  virtual void A_inv_chi(Spinor *const out, Spinor const *const in, int const isign,
                         int const cb) = 0;

  /**
    Forwarder for the `dslash`.

    This will call the `dslash` function of the respective QPhiX dslash class. There is a subtle
    difference between the Wilson and all other cases. The Wilson dslash is just the hopping matrix,
    just the operator \f$ D \f$. For every other case (clover, twisted mass, twisted mass clover),
    the `dslash` member function will compute \f$ A^{-1} D \f$. In the Wilson case, this \f$ A =
    \alpha = 4 + m = 1/(2 \kappa) \f$. Since that is _not_ included in the Wilson `dslash`, you will
    obtain different results when using WilsonDslash::dslash and WilsonTMDslash::dslash with \f$
    \mu = 0 \f$.

    \todo Make this member function `const`. For this the member function in
    QPhiX that is called internally must be marked `const` as well.
    */
  virtual void dslash(Spinor *const res, const Spinor *const psi, const SU3MatrixBlock *const u,
                      int const isign, int const cb) = 0;

  /**
    Always plain Wilson dslash.

    In contrast to the \ref dslash member function which just forwards the implementation of QPhiX,
    this will always give you the “naked” plain Wilson dslash without any factors of \f$ A^{-1} \f$
    applied.
    */
  virtual void plain_dslash(Spinor *const res, const Spinor *const psi,
                            const SU3MatrixBlock *const u, int const isign, int const cb) {
    // XXX Perhaps rather implement this with an instance of the WilsonDslash instead?

    auto tmp = QPhiX::makeFourSpinorHandle(*geom);
    dslash(tmp.get(), psi, u, isign, cb);
    A_chi(res, tmp.get(), isign, cb);
  };

  /**
    Always “dressed” dslash.

    This computes \f$ A^{-1} D \f$ for all variants. In the Wilson case, this will give \f$
    \alpha^{-1} D \f$.
    */
  virtual void A_inv_dslash(Spinor *const res, const Spinor *const psi,
                            const SU3MatrixBlock *const u, int const isign, int const cb) {
    dslash(res, psi, u, isign, cb);
  };

  /**
    Forwarder for the `achimbdpsi`.

    \todo Make this member function `const`. For this the member function in QPhiX that is called
    internally must be marked `const` as well.
    */
  virtual void achimbdpsi(Spinor *const res, const Spinor *const psi, const Spinor *const chi,
                          const SU3MatrixBlock *const u, double const alpha, double const beta,
                          int const isign, int const cb) = 0;

  /**
    Prepares the sources on the odd checkerboard.

    This computes
    \f[
        \tilde b_o = \frac 12 D_{oe} M_{ee}^{-1} b_e + b_o \,.
    \f]

    \param[out] tilde_b_odd Prepared source
    \param[in] b_even Source (right hand side) on the even lattice sites
    \param]in] b_odd Source on the odd lattice sites
    \param[in] u Gauge field on the odd lattice sites
    */
  virtual void prepare_source(Spinor *const tilde_b_odd, Spinor const *const b_even,
                              Spinor const *const b_odd, SU3MatrixBlock const *const u);

  /**
    Reconstructs the solution on the even lattices sites.

    This computes
    \f[
        x_e = M_{ee}^{-1} \left( b_e - \frac 12 D_{eo} x_o \right) \,.
    \f]

    \param[out] x_even Solution on the even lattices sites
    \param[in] b_even Source (right hand side) on the even lattice sites
    \param[in] x_odd Solution on the odd lattices sites
    \param[in] u Gauge field on the even lattice sites
    */
  virtual void reconstruct_solution(Spinor *const x_even, Spinor const *const b_even,
                                    Spinor const *const x_odd, SU3MatrixBlock const *const u);

  Geom *getGeometry() const { return geom; }

 private:
  Geom *const geom;

  double const t_boundary;
  double const aniso_coeff_S;
  double const aniso_coeff_T;
  double const mass;
};

template <typename FT, int veclen, int soalen, bool compress12>
class WilsonDslash : public Dslash<FT, veclen, soalen, compress12> {
 public:
  typedef typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock Spinor;
  typedef typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::SU3MatrixBlock SU3MatrixBlock;

  WilsonDslash(::QPhiX::Geometry<FT, veclen, soalen, compress12> *geom_, double const t_boundary_,
               double const aniso_coeff_S_, double const aniso_coeff_T_, double const mass_)
      : Dslash<FT, veclen, soalen, compress12>(geom_, t_boundary_, aniso_coeff_S_, aniso_coeff_T_,
                                               mass_),
        upstream_dslash(geom_, t_boundary_, aniso_coeff_S_, aniso_coeff_T_),
        mass_factor_alpha(4.0 + mass_),
        mass_factor_beta(1.0 / (4.0 * mass_factor_alpha)) {}

  void A_chi(Spinor *const out, Spinor const *const in, int const isign_ignored,
             int const cb_ignored) override {
    int const n_blas_simt = 1;
    ::QPhiX::axy(mass_factor_alpha, in, out, upstream_dslash.getGeometry(), n_blas_simt);
  }

  void A_inv_chi(Spinor *const out, Spinor const *const in, int const isign_ignored,
                 int const cb_ignored) override {
    int const n_blas_simt = 1;
    ::QPhiX::axy(1.0 / mass_factor_alpha, in, out, upstream_dslash.getGeometry(), n_blas_simt);
  }

  void dslash(Spinor *const res, const Spinor *const psi, const SU3MatrixBlock *const u,
              int const isign, int const cb) override {
    upstream_dslash.dslash(res, psi, u, isign, cb);
  }

  void plain_dslash(Spinor *const res, const Spinor *const psi, const SU3MatrixBlock *const u,
                    int const isign, int const cb) override {
    dslash(res, psi, u, isign, cb);
  };

  void A_inv_dslash(Spinor *const res, const Spinor *const psi, const SU3MatrixBlock *const u,
                    int const isign, int const cb) override {
    auto tmp = QPhiX::makeFourSpinorHandle(upstream_dslash.getGeometry());
    dslash(tmp.get(), psi, u, isign, cb);
    A_inv_chi(res, tmp.get(), isign, cb);
  };

  void achimbdpsi(Spinor *const res, const Spinor *const psi, const Spinor *const chi,
                  const SU3MatrixBlock *const u, double const alpha, double const beta,
                  int const isign, int const cb) override {
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

  WilsonTMDslash(::QPhiX::Geometry<FT, veclen, soalen, compress12> *geom_, double const t_boundary_,
                 double const aniso_coeff_S_, double const aniso_coeff_T_, double const mass_,
                 double const twisted_mass_)
      : Dslash<FT, veclen, soalen, compress12>(geom_, t_boundary_, aniso_coeff_S_, aniso_coeff_T_,
                                               mass_),
        upstream_dslash(geom_, t_boundary_, aniso_coeff_S_, aniso_coeff_T_, mass_, twisted_mass_),
        mass_factor_alpha(4.0 + mass_),
        mass_factor_beta(0.25),
        derived_mu(twisted_mass_ / mass_factor_alpha),
        derived_mu_inv(mass_factor_alpha /
                       (mass_factor_alpha * mass_factor_alpha + twisted_mass_ * twisted_mass_)) {}

  void A_chi(Spinor *const out, Spinor const *const in, int const isign,
             int const cb_ignored) override {
    helper_A_chi(out, in, -derived_mu * isign, mass_factor_alpha);
  }

  void A_inv_chi(Spinor *const out, Spinor const *const in, int const isign,
                 int const cb_ignored) override {
    helper_A_chi(out, in, derived_mu * isign, derived_mu_inv);
  }

  void dslash(Spinor *const res, const Spinor *const psi, const SU3MatrixBlock *const u,
              int const isign, int const cb) override {
    upstream_dslash.dslash(res, psi, u, isign, cb);
  }

  void achimbdpsi(Spinor *const res, const Spinor *const psi, const Spinor *const chi,
                  const SU3MatrixBlock *const u, double const alpha, double const beta,
                  int const isign, int const cb) override {
    upstream_dslash.dslashAChiMinusBDPsi(res, psi, chi, u, alpha, beta, isign, cb);
  }

 private:
  void helper_A_chi(Spinor *const out, Spinor const *const in, double const factor_a,
                    double const factor_b);

  ::QPhiX::TMDslash<FT, veclen, soalen, compress12> upstream_dslash;

  double const mass_factor_alpha;
  double const mass_factor_beta;
  double const derived_mu;
  double const derived_mu_inv;
};

template <typename FT, int veclen, int soalen, bool compress12>
class WilsonClovDslash : public Dslash<FT, veclen, soalen, compress12> {
 public:
  typedef typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock Spinor;
  typedef typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::SU3MatrixBlock SU3MatrixBlock;
  typedef typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::CloverBlock CloverBlock;

  WilsonClovDslash(::QPhiX::Geometry<FT, veclen, soalen, compress12> *geom_,
                   double const t_boundary_, double const aniso_coeff_S_,
                   double const aniso_coeff_T_, double const mass_, CloverBlock * const(&clover_)[2],
                   CloverBlock * const(&inv_clover_)[2])
      : Dslash<FT, veclen, soalen, compress12>(geom_, t_boundary_, aniso_coeff_S_, aniso_coeff_T_,
                                               mass_),
        upstream_dslash(geom_, t_boundary_, aniso_coeff_S_, aniso_coeff_T_),
        mass_factor_alpha(4.0 + mass_),
        mass_factor_beta(1.0 / (4.0 * mass_factor_alpha)) 
        {
          for(int cb : {0, 1} ){
            clover[cb] = clover_[cb];
            inv_clover[cb] = inv_clover_[cb];
          }
        }

  void A_chi(Spinor *const out, Spinor const *const in, int const isign_ignored,
             int const cb) override {
    clover_product(out, in, clover[cb], upstream_dslash.getGeometry());
  }

  void A_inv_chi(Spinor *const out, Spinor const *const in, int const isign_ignored,
                 int const cb) override {
    clover_product(out, in, inv_clover[cb], upstream_dslash.getGeometry());
  }

  void dslash(Spinor *const res, const Spinor *const psi, const SU3MatrixBlock *const u,
              int const isign, int const cb) override {
    upstream_dslash.dslash(res, psi, u, inv_clover[cb], isign, cb);
  }

  void achimbdpsi(Spinor *const res, const Spinor *const psi, const Spinor *const chi,
                  const SU3MatrixBlock *const u, double const alpha, double const beta,
                  int const isign, int const cb) override {
    upstream_dslash.dslashAChiMinusBDPsi(res, psi, chi, u, clover[cb], mass_factor_beta, isign, cb);
  }

 private:
  ::QPhiX::ClovDslash<FT, veclen, soalen, compress12> upstream_dslash;

  double const mass_factor_alpha;
  double const mass_factor_beta;

  /**
    Reference to the clover term.

    This class has to provide a `dslash` and `achimbdpsi` member function with the prescribed
    argument list which does not contain the clover term. The user of these classes should not have
    to differentiate between non-clover and clover variants. In order to provide the function
    signature, the clover term is a member. This means that the user has to construct a new operator
    each time the clover term has changed.
    */
  CloverBlock * clover[2];

  /// See \ref clover.
  CloverBlock * inv_clover[2];
};

template <typename FT, int veclen, int soalen, bool compress12>
class WilsonClovTMDslash : public Dslash<FT, veclen, soalen, compress12> {
 public:
  typedef typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock Spinor;
  typedef typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::SU3MatrixBlock SU3MatrixBlock;
  typedef
      typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FullCloverBlock FullCloverBlock;

  WilsonClovTMDslash(::QPhiX::Geometry<FT, veclen, soalen, compress12> *geom_,
                     double const t_boundary_, double const aniso_coeff_S_,
                     double const aniso_coeff_T_, double const mass_, double const twisted_mass_,
                     FullCloverBlock * const (&clover_)[2][2], FullCloverBlock * const(&inv_clover_)[2][2] )
      : Dslash<FT, veclen, soalen, compress12>(geom_, t_boundary_, aniso_coeff_S_, aniso_coeff_T_,
                                               mass_),
        upstream_dslash(geom_, t_boundary_, aniso_coeff_S_, aniso_coeff_T_),
        mass_factor_alpha(4.0 + mass_),
        mass_factor_beta(0.25),
        derived_mu(twisted_mass_ / mass_factor_alpha),
        derived_mu_inv(mass_factor_alpha /
                       (mass_factor_alpha * mass_factor_alpha + twisted_mass_ * twisted_mass_))
        {
          for( int cb : { 0, 1 } ){
            for( int fl : {0, 1 } ){
              clover[cb][fl] = clover_[cb][fl];
              inv_clover[cb][fl] = inv_clover_[cb][fl];
            }
          }
        }

  void A_chi(Spinor *const out, Spinor const *const in, int const isign, int const cb) override {
    QPhiX::masterPrintf("clover[%d][0] -> %p, clover[%d][1] -> %p\n", 
                        cb, (void*)clover[cb][0], cb, (void*)clover[cb][1]);
    if( isign == -1 ){
      full_clover_product(out, in, clover[cb][1], upstream_dslash.getGeometry());
    } else {
      full_clover_product(out, in, clover[cb][0], upstream_dslash.getGeometry());
    }
  }

  void A_inv_chi(Spinor *const out, Spinor const *const in, int const isign,
                 int const cb) override {
    if( isign == -1 ){
      full_clover_product(out, in, inv_clover[cb][1], upstream_dslash.getGeometry());
    } else {
      full_clover_product(out, in, inv_clover[cb][0], upstream_dslash.getGeometry());
    }
  }

  void dslash(Spinor *const res, const Spinor *const psi, const SU3MatrixBlock *const u,
              int const isign, int const cb) override {
    upstream_dslash.dslash(res, psi, u, (const FullCloverBlock **)inv_clover[cb], isign, cb);
  }

  void achimbdpsi(Spinor *const res, const Spinor *const psi, const Spinor *const chi,
                  const SU3MatrixBlock *const u, double const alpha, double const beta,
                  int const isign, int const cb) override {
    upstream_dslash.dslashAChiMinusBDPsi(res, psi, chi, u, (const FullCloverBlock **)clover[cb], mass_factor_beta, isign, cb);
  }

 private:
  ::QPhiX::TMClovDslash<FT, veclen, soalen, compress12> upstream_dslash;

  double const mass_factor_alpha;
  double const mass_factor_beta;
  double const derived_mu;
  double const derived_mu_inv;

  FullCloverBlock * clover[2][2];
  FullCloverBlock * inv_clover[2][2];
};

template <typename FT, int veclen, int soalen, bool compress12>
void WilsonTMDslash<FT, veclen, soalen, compress12>::helper_A_chi(Spinor *const out,
                                                                  Spinor const *const in,
                                                                  double const factor_a,
                                                                  double const factor_b) {
  auto const nVecs = upstream_dslash.getGeometry().nVecs();
  auto const Pxy = upstream_dslash.getGeometry().getPxy();
  auto const Pxyz = upstream_dslash.getGeometry().getPxyz();

  for (uint64_t t = 0; t < T; t++)
    for (uint64_t x = 0; x < LX / 2; x++)
      for (uint64_t y = 0; y < LY; y++)
        for (uint64_t z = 0; z < LZ; z++) {
          uint64_t const SIMD_vector = x / soalen;
          uint64_t const x_internal = x % soalen;
          uint64_t const qphix_idx = t * Pxyz + z * Pxy + y * nVecs + SIMD_vector;

          for (int color = 0; color < 3; ++color) {
            for (int spin_block = 0; spin_block < 2; ++spin_block) {
              // Implement the $\gamma_5$ structure.
              auto const signed_factor_a = factor_a * (spin_block == 0 ? 1.0 : -1.0);

              for (int half_spin = 0; half_spin < 2; ++half_spin) {
                auto const four_spin = 2 * spin_block + half_spin;
                for (int v = 0; v < soalen; ++v) {
                  auto &out_bcs = out[qphix_idx][color][four_spin];
                  auto const &in_bcs = in[qphix_idx][color][four_spin];

                  out_bcs[re][v] = factor_b * (in_bcs[re][v] + signed_factor_a * in_bcs[im][v]);
                  out_bcs[im][v] = factor_b * (in_bcs[im][v] - signed_factor_a * in_bcs[re][v]);
                }
              }
            }
          }

        }  // volume
};

template <typename FT, int veclen, int soalen, bool compress12>
void Dslash<FT, veclen, soalen, compress12>::prepare_source(Spinor *const tilde_b_odd,
                                                            Spinor const *const b_even,
                                                            Spinor const *const b_odd,
                                                            SU3MatrixBlock const *const u) {
  auto Mee_be = QPhiX::makeFourSpinorHandle(*geom);
  WilsonDslash<FT, veclen, soalen, compress12> plain_dslash(geom, t_boundary, aniso_coeff_S,
                                                            aniso_coeff_T, mass);

  A_inv_chi(Mee_be.get(), b_even, 1, cb_even);

  plain_dslash.dslash(tilde_b_odd, Mee_be.get(), u, 1, cb_odd);

  // FIXME Perhaps use a variable number of BLAS threads here (last parameter).
  QPhiX::aypx(0.5, Mee_be.get(), tilde_b_odd, *geom, 1);
}

template <typename FT, int veclen, int soalen, bool compress12>
void Dslash<FT, veclen, soalen, compress12>::reconstruct_solution(Spinor *const x_even,
                                                                  Spinor const *const b_even,
                                                                  Spinor const *const x_odd,
                                                                  SU3MatrixBlock const *const u) {
  auto tmp = QPhiX::makeFourSpinorHandle(*geom);
  WilsonDslash<FT, veclen, soalen, compress12> plain_dslash(geom, t_boundary, aniso_coeff_S,
                                                            aniso_coeff_T, mass);

  plain_dslash.dslash(tmp.get(), x_odd, u, 1, cb_even);
  QPhiX::aypx(0.5, b_even, tmp.get(), *geom, 1);
  A_inv_chi(x_even, tmp.get(), 1, cb_even);
}
}
