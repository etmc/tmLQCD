/***********************************************************************
 * Copyright (C) 2008,2009,2010,2011,2012  
 * Andreas Stathopoulos, Kostas Orginos, Abdou M. Abdel-Rehim
 *
 * This program is based on interfacing the eigCG solver to the tmLQCD code.
 * It was written by Abdou M. Abdel-Rehim based on the original code written
 * by Andreas Stathopoulos and Kostas Orginos and uses functions written in 
 * tmLQCD by Carsten Urbach 
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
/***********************************************************************
 * Subroutine restart_X - This subroutine computes X*hVecs and places 
 *    the result in X.
 ***********************************************************************/



#ifndef _RESTART_X_H
#define _RESTART_X_H


/*double precision version */
void Zrestart_X(_Complex double  *X, int ldx, _Complex double  *hVecs, int nLocal, 
                int basisSize, int restartSize, _Complex double  *rwork, int rworkSize);


#endif
