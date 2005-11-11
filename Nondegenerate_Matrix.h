/* $Id$ */
#ifndef _NONDEGENRATE_MATRIX_H
#define _NONDEGENRATE_MATRIX_H

void QNon_degenerate(spinor * const l_strange, spinor * const l_charm,
                     spinor * const k_strange,  spinor * const k_charm);

void QdaggerNon_degenerate(spinor * const l_strange, spinor * const l_charm,
                           spinor * const k_strange, spinor * const k_charm);

void Q_Qdagger_ND(spinor * const l_strange, spinor * const l_charm,
                  spinor * const k_strange, spinor * const k_charm);

void Q_Qdagger_ND_BI(bispinor * const bisp_l, bispinor * const bisp_k);

void Q_test_epsilon(spinor * const l_strange, spinor * const l_charm,
                    spinor * const k_strange, spinor * const k_charm);

#endif
