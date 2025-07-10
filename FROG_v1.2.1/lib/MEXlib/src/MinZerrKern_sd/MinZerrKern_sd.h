/**	\file
 *
 *	$Author: rana $
 *	$Date: 2019-09-24  $
 *	$Revision: 1.1 $
 */

#ifndef _MinZerrKern_sdH_
#define _MinZerrKern_sdH_

#include "mexArray.h"

//!	Calculates the coefficients.
void MinZerrKern_sd(const TmexArray &Esig, const TmexArray &Et, const TmexArray &dZ, TmexArray &X);

#endif //_MinZerrKern_sdH_
