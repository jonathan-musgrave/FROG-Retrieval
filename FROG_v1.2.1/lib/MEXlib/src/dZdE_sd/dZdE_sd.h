/**	\file
 *
 *	$Author: rana $
 *	$Date: 2019-09-24 $
 *	$Revision: 1.1 $
 */

#ifndef _dZdE_sdH_
#define _dZdE_sdH_

#include "mexArray.h"

//! Calculate the gradient of the Z error.
void dZdE_sd(const TmexArray &Esigp, const TmexArray &Et, TmexArray &dZ);

#endif //_dZdE_sdH_
