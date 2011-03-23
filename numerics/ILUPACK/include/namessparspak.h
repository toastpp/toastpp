#ifndef _NAMESSPARSPAK_H
#define _NAMESSPARSPAK_H

#include "f2c.h"


/* on several architectures names of fortran routines are passed to C in 
   different ways. To cover this different architectures use in C only lower
   letters for the fortran names. Dependent on the switch you use they are
   replaced by the correct function name
*/

/* only use capital letters */
#if defined __CAPS__ && !defined __UNDERSCORE__ && !defined __2UNDERSCORES__
#define bshufl          BSHUFL
#define copysi          COPYSI
#define degree          DEGREE
#define elslv           ELSLV
#define esfct           ESFCT
#define euslv           EUSLV
#define fnbenv          FNBENV
#define fndsep          FNDSEP
#define fnenv           FNENV
#define fnlvls          FNLVLS
#define fnofnz          FNOFNZ
#define fnroot          FNROOT
#define fnspan          FNSPAN
#define fntadj          FNTADJ
#define fntenv          FNTENV
#define fn1wd           FN1WD
#define gennd           GENND
#define genqmd          GENQMD
#define genrcm          GENRCM
#define genrqt          GENRQT
#define gen1wd          GEN1WD
#define gsfct           GSFCT
#define gsslv           GSSLV
#define guslv           GUSLV
#define glslv           GLSLV
#define qmdmrg          QMDMRG
#define qmdqt           QMDQT
#define qmdrch          QMDRCH
#define qmdupd          QMDUPD
#define rcm             RCM
#define reach           REACH
#define rootls          ROOTLS
#define rqtree          RQTREE
#define smbfct          SMBFCT
#define sorts1          SORTS1
#define subrcm          SUBRCM
#define tsfct           TSFCT
#define tsslv           TSSLV
#define invrse          INVRSE
#define revrse          REVRSE
#define permrv          PERMRV
#define chlsetup        CHLSETUP
#define suslv           SUSLV
#define slslv           SLSLV
#define fuslv           FUSLV
#define flslv           FLSLV

#define spkops          SPKOPS


/* no capital letters but underscores */
#elif defined __UNDERSCORE__ && !defined __CAPS__ && !defined __2UNDERSCORES__
#define bshufl          bshufl_
#define copysi          copysi_
#define degree          degree_
#define elslv           elslv_
#define esfct           esfct_
#define euslv           euslv_
#define fnbenv          fnbenv_
#define fndsep          fndsep_
#define fnenv           fnenv_
#define fnlvls          fnlvls_
#define fnofnz          fnofnz_
#define fnroot          fnroot_
#define fnspan          fnspan_
#define fntadj          fntadj_
#define fntenv          fntenv_
#define fn1wd           fn1wd_
#define gennd           gennd_
#define genqmd          genqmd_
#define genrcm          genrcm_
#define genrqt          genrqt_
#define gen1wd          gen1wd_
#define gsfct           gsfct_
#define gsslv           gsslv_
#define guslv           guslv_
#define glslv           glslv_
#define qmdmrg          qmdmrg_
#define qmdqt           qmdqt_
#define qmdrch          qmdrch_
#define qmdupd          qmdupd_
#define rcm             rcm_
#define reach           reach_
#define rootls          rootls_
#define rqtree          rqtree_
#define smbfct          smbfct_
#define sorts1          sorts1_
#define subrcm          subrcm_
#define tsfct           tsfct_
#define tsslv           tsslv_
#define invrse          invrse_
#define revrse          revrse_
#define permrv          permrv_
#define chlsetup        chlsetup_
#define suslv           suslv_
#define slslv           slslv_
#define fuslv           fuslv_
#define flslv           flslv_

#define spkops          spkops_


/* both are defined */
#elif defined __CAPS__ && defined __UNDERSCORE__ && !defined __2UNDERSCORES__
#define bshufl          BSHUFL_
#define copysi          COPYSI_
#define degree          DEGREE_
#define elslv           ELSLV_
#define esfct           ESFCT_
#define euslv           EUSLV_
#define fnbenv          FNBENV_
#define fndsep          FNDSEP_
#define fnenv           FNENV_
#define fnlvls          FNLVLS_
#define fnofnz          FNOFNZ_
#define fnroot          FNROOT_
#define fnspan          FNSPAN_
#define fntadj          FNTADJ_
#define fntenv          FNTENV_
#define fn1wd           FN1WD_
#define gennd           GENND_
#define genqmd          GENQMD_
#define genrcm          GENRCM_
#define genrqt          GENRQT_
#define gen1wd          GEN1WD_
#define gsfct           GSFCT_
#define gsslv           GSSLV_
#define guslv           GUSLV_
#define glslv           GLSLV_
#define qmdmrg          QMDMRG_
#define qmdqt           QMDQT_
#define qmdrch          QMDRCH_
#define qmdupd          QMDUPD_
#define rcm             RCM_
#define reach           REACH_
#define rootls          ROOTLS_
#define rqtree          RQTREE_
#define smbfct          SMBFCT_
#define sorts1          SORTS1_
#define subrcm          SUBRCM_
#define tsfct           TSFCT_
#define tsslv           TSSLV_
#define invrse          INVRSE_
#define revrse          REVRSE_
#define permrv          PERMRV_
#define chlsetup        CHLSETUP_
#define suslv           SUSLV_
#define slslv           SLSLV_
#define fuslv           FUSLV_
#define flslv           FLSLV_

#define spkops          SPKOPS_


/* CAPS and 2 underscores are defined */
#elif defined __CAPS__ && defined __2UNDERSCORES__
#define bshufl          BSHUFL__
#define copysi          COPYSI__
#define degree          DEGREE__
#define elslv           ELSLV__
#define esfct           ESFCT__
#define euslv           EUSLV__
#define fnbenv          FNBENV__
#define fndsep          FNDSEP__
#define fnenv           FNENV__
#define fnlvls          FNLVLS__
#define fnofnz          FNOFNZ__
#define fnroot          FNROOT__
#define fnspan          FNSPAN__
#define fntadj          FNTADJ__
#define fntenv          FNTENV__
#define fn1wd           FN1WD__
#define gennd           GENND__
#define genqmd          GENQMD__
#define genrcm          GENRCM__
#define genrqt          GENRQT__
#define gen1wd          GEN1WD__
#define gsfct           GSFCT__
#define gsslv           GSSLV__
#define guslv           GUSLV__
#define glslv           GLSLV__
#define qmdmrg          QMDMRG__
#define qmdqt           QMDQT__
#define qmdrch          QMDRCH__
#define qmdupd          QMDUPD__
#define rcm             RCM__
#define reach           REACH__
#define rootls          ROOTLS__
#define rqtree          RQTREE__
#define smbfct          SMBFCT__
#define sorts1          SORTS1__
#define subrcm          SUBRCM__
#define tsfct           TSFCT__
#define tsslv           TSSLV__
#define invrse          INVRSE__
#define revrse          REVRSE__
#define permrv          PERMRV__
#define chlsetup        CHLSETUP__
#define suslv           SUSLV__
#define slslv           SLSLV__
#define fuslv           FUSLV__
#define flslv           FLSLV__

#define spkops          SPKOPS__


/* no capital letters but 2 underscores */
#elif defined __2UNDERSCORES__ && !defined __CAPS__
#define bshufl          bshufl__
#define copysi          copysi__
#define degree          degree__
#define elslv           elslv__
#define esfct           esfct__
#define euslv           euslv__
#define fnbenv          fnbenv__
#define fndsep          fndsep__
#define fnenv           fnenv__
#define fnlvls          fnlvls__
#define fnofnz          fnofnz__
#define fnroot          fnroot__
#define fnspan          fnspan__
#define fntadj          fntadj__
#define fntenv          fntenv__
#define fn1wd           fn1wd__
#define gennd           gennd__
#define genqmd          genqmd__
#define genrcm          genrcm__
#define genrqt          genrqt__
#define gen1wd          gen1wd__
#define gsfct           gsfct__
#define gsslv           gsslv__
#define guslv           guslv__
#define glslv           glslv__
#define qmdmrg          qmdmrg__
#define qmdqt           qmdqt__
#define qmdrch          qmdrch__
#define qmdupd          qmdupd__
#define rcm             rcm__
#define reach           reach__
#define rootls          rootls__
#define rqtree          rqtree__
#define smbfct          smbfct__
#define sorts1          sorts1__
#define subrcm          subrcm__
#define tsfct           tsfct__
#define tsslv           tsslv__
#define invrse          invrse__
#define revrse          revrse__
#define permrv          permrv__
#define chlsetup        chlsetup__
#define suslv           suslv__
#define slslv           slslv__
#define fuslv           fuslv__
#define flslv           flslv__

#define spkops          spkops__

// nothing at all
#else
#define bshufl          bshufl
#define copysi          copysi
#define degree          degree
#define elslv           elslv
#define esfct           esfct
#define euslv           euslv
#define fnbenv          fnbenv
#define fndsep          fndsep
#define fnenv           fnenv
#define fnlvls          fnlvls
#define fnofnz          fnofnz
#define fnroot          fnroot
#define fnspan          fnspan
#define fntadj          fntadj
#define fntenv          fntenv
#define fn1wd           fn1wd
#define gennd           gennd
#define genqmd          genqmd
#define genrcm          genrcm
#define genrqt          genrqt
#define gen1wd          gen1wd
#define gsfct           gsfct
#define gsslv           gsslv
#define guslv           guslv
#define glslv           glslv
#define qmdmrg          qmdmrg
#define qmdqt           qmdqt
#define qmdrch          qmdrch
#define qmdupd          qmdupd
#define rcm             rcm
#define reach           reach
#define rootls          rootls
#define rqtree          rqtree
#define smbfct          smbfct
#define sorts1          sorts1
#define subrcm          subrcm
#define tsfct           tsfct
#define tsslv           tsslv
#define invrse          invrse
#define revrse          revrse
#define permrv          permrv
#define chlsetup        chlsetup
#define suslv           suslv
#define slslv           slslv
#define fuslv           fuslv
#define flslv           flslv

#define spkops          spkops

#endif /* defined __CAPS__ && ... */



#endif /* _NAMESSPARSPAK_H */

