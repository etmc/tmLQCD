/* $Id$ */
#ifndef _NONDEGENRATE_MATRIX_H
#define _NONDEGENRATE_MATRIX_H

void QNon_degenerate(spinor * const l_strange, spinor * const l_charm,
                     spinor * const k_strange,  spinor * const k_charm);

void QdaggerNon_degenerate(spinor * const l_strange, spinor * const l_charm,
                           spinor * const k_strange, spinor * const k_charm);

void QNon_degenerate_eo(spinor * const l_strange, spinor * const l_charm,
                        spinor * const k_strange, spinor * const k_charm);

#endif
