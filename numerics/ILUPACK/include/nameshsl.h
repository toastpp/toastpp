#ifndef _NAMESHSL_H
#define _NAMESHSL_H

#include "f2c.h"

/* on several architectures names of fortran routines are passed to C in 
   different ways. To cover this different architectures use in C only lower
   letters for the fortran names. Dependent on the switch you use they are
   replaced by the correct function name
*/

/* only use capital letters */
#if defined __CAPS__ && !defined __UNDERSCORE__ && !defined __2UNDERSCORES__
#define mc64ad          MC64AD
#define mc64id          MC64ID
#define fd64ad          FD64AD
#define mc64as          MC64A
#define mc64is          MC64I
#define fd64as          FD64A

#define mc64az          MC64AZ
#define mc64iz          MC64IZ
#define mc64ac          MC64AC
#define mc64ic          MC64IC

/* no capital letters */
#elif defined __UNDERSCORE__ && !defined __CAPS__ && !defined __2UNDERSCORES__
#define mc64ad          mc64ad_
#define mc64id          mc64id_
#define fd64ad          fd64ad_
#define mc64as          mc64a_
#define mc64is          mc64i_
#define fd64as          fd64a_

#define mc64az          mc64az_
#define mc64iz          mc64iz_
#define mc64ac          mc64ac_
#define mc64ic          mc64ic_

/* both are defined */
#elif defined __CAPS__ && defined __UNDERSCORE__ && !defined __2UNDERSCORES__
#define mc64ad          MC64AD_
#define mc64id          MC64ID_
#define fd64ad          FD64AD_
#define mc64as          MC64A_
#define mc64is          MC64I_
#define fd64as          FD64A_

#define mc64az          MC64AZ_
#define mc64iz          MC64IZ_
#define mc64ac          MC64AC_
#define mc64ic          MC64IC_


/* CAPS and 2 underscores are defined */
#elif defined __CAPS__ && defined __2UNDERSCORES__
#define mc64ad          MC64AD__
#define mc64id          MC64ID__
#define fd64ad          FD64AD__
#define mc64as          MC64A__
#define mc64is          MC64I__
#define fd64as          FD64A__

#define mc64az          MC64AZ__
#define mc64iz          MC64IZ__
#define mc64ac          MC64AC__
#define mc64ic          MC64IC__


/* no capital letters but 2 underscores */
#elif defined __2UNDERSCORES__ && !defined __CAPS__
#define mc64ad          mc64ad__
#define mc64id          mc64id__
#define fd64ad          fd64ad__
#define mc64as          mc64a__
#define mc64is          mc64i__
#define fd64as          fd64a__

#define mc64az          mc64az__
#define mc64iz          mc64iz__
#define mc64ac          mc64ac__
#define mc64ic          mc64ic__

#else
#define mc64ad          mc64ad
#define mc64id          mc64id
#define fd64ad          fd64ad
#define mc64as          mc64a
#define mc64is          mc64i
#define fd64as          fd64a

#define mc64az          mc64az
#define mc64iz          mc64iz
#define mc64ac          mc64ac
#define mc64ic          mc64ic

#endif /* defined __CAPS__ && ... */

#endif /* _NAMESHSL_H */


