# decide which type of matching should be used
MATCHING=-D_MUMPS_MATCHING_
#MATCHING=-D_PARDISO_MATCHING_
#MATCHING=-D_MC64_MATCHING_


# sources

ifeq ($(MATCHING),-D_MC64_MATCHING_)
NOTDISTRIBUTED=$(DIRNOTDISTRIBUTED)MC21S.o   \
               $(DIRNOTDISTRIBUTED)MC21D.o  \
               $(DIRNOTDISTRIBUTED)MC64S.o  \
               $(DIRNOTDISTRIBUTED)MC64D.o
else
# the codes mentioned above are not distributed
# the definition is only done for convenience
#NOTDISTRIBUTED=$(DIRNOTDISTRIBUTED)AMF4.o # $(DIRNOTDISTRIBUTED)mps_pardiso_dummy.o   
#NOTDISTRIBUTED=$(DIRNOTDISTRIBUTED)AMF4.o 
endif

ifeq ($(MATCHING),-D_MUMPS_MATCHING_)
MUMPSLIBS=-lmumps 
MUMPSDLIB=-lmumps
MUMPSSLIB=-lmumps
endif
ifeq ($(MATCHING),-D_PARDISO_MATCHING_)
PARDISOLIB=$(PARDISO)
endif

LIBPRECONDITIONERS=$(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)AMGsetup.o       \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)AMGextract.o     \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)AMGsol1.o        \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)AMGsol2.o        \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)AMGsol.o         \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)AMGtsol.o        \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)ildlc.o          \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)ildlcsol.o       \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)pildlcdlsol.o    \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)pildlcdusol.o    \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)pildlclsol.o     \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)pildlcusol.o     \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)ildlcs.o         \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)ildlcsols.o      \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)ildlcz.o         \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)ildlcsolz.o      \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)ildlczs.o        \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)ildlcsolzs.o     \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)iluc.o           \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)ilucsol.o        \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)iluctsol.o       \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)ilucdlsol.o      \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)iluctdlsol.o     \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)ilucdusol.o      \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)iluctdusol.o     \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)iluclsol.o       \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)iluctlsol.o      \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)ilucusol.o       \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)iluctusol.o      \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)LUPQ.o           \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)LUPQsol.o        \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)LUPQtsol.o       \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)LDLP.o           \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)LDLPsol.o        \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)symAMGextract.o  \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)ssmAMGextract.o  \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)shrAMGextract.o  \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)csymAMGextract.o \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)spdAMGsol.o      \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)spdAMGsol1.o     \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)spdAMGsol2.o     \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)spdAMGsetup.o    \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)pildlc.o         \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)mpildlc.o        \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)piluc.o          \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)spiluc.o         \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)mpiluc.o         \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)sympiluc.o       \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)sympilucs.o      \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)sympilucsol.o    \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)sympilucsols.o   \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)sympiluclsol.o   \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)sympiluclsols.o  \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)sympilucusol.o   \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)sympilucusols.o  \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)psympiluclsol.o  \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)psympiluclsols.o \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)psympilucusol.o  \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)psympilucusols.o \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)AMGlsol.o        \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)AMGdusol.o       \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)AMGdlsol.o       \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)AMGusol.o        \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)AMGtlsol.o       \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)AMGtusol.o       \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)AMGtdlsol.o      \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)AMGtdusol.o      \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)LUPQlsol.o       \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)LUPQusol.o       \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)LUPQdlsol.o      \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)LUPQdusol.o      \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)LUPQtlsol.o      \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)LUPQtusol.o      \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)LUPQtdlsol.o     \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)LUPQtdusol.o     \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)piluclsol.o      \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)pilucusol.o      \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)pilucdlsol.o     \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)pilucdusol.o     \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)piluctlsol.o     \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)piluctusol.o     \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)piluctdlsol.o    \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)piluctdusol.o    \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)symAMGsetup.o    \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)symAMGsetups.o   \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)symAMGsol.o      \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)symAMGsols.o     \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)symAMGsol1.o     \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)symAMGsol1s.o    \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)symAMGsol2.o     \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)symAMGsol2s.o    \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)symiluc.o        \
                   $(DIRPRECONDITIONERS)$(FLOAT)/$(FLOAT)symilucs.o     

OBJECTPRECONDITIONERS=

PRECONDITIONERS=      $(LIBPRECONDITIONERS) $(OBJECTPRECONDITIONERS)



LIBORDERINGS=   $(DIRORDERINGS)$(FLOAT)/$(FLOAT)qsortR2I.o                \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)indset.o                  \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permindset.o              \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permrcm.o                 \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)perm_rcm_FC.o             \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)perm_rcm_FCv.o            \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_rcm_FC.o          \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_rcm_FCv.o         \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_rcm_FCs.o         \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_rcm_FCvs.o        \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permnd.o                  \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permnd_FC.o               \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permnd_FCv.o              \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_nd_FC.o           \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_nd_FCv.o          \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_nd_FCs.o          \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_nd_FCvs.o         \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permamf.o                 \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)perm_amf_FC.o             \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)perm_amf_FCv.o            \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_amf_FC.o          \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_amf_FCv.o         \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_amf_FCs.o         \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_amf_FCvs.o        \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permmmd.o                 \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)perm_mmd_FC.o             \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)perm_mmd_FCv.o            \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mmd_FC.o          \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mmd_FCv.o         \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mmd_FCs.o         \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mmd_FCvs.o        \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permnull.o                \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permPQ.o                  \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)spdpermrcm.o              \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)spdpermPP.o               \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)indPPF3.o                 \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)indPQF5.o                 \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)indFC.o                   \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permFC.o                  \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symindFC.o                \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympindFC.o               \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympindFCs.o              \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympindFCv.o              \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympindFCvs.o             \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symbuildblock.o           \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symbuildblocks.o          \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermFC.o               \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)indFCv.o                  \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permFCv.o                 \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symindFCv.o               \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symindFCvs.o              \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermFCv.o              \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)perm_metis_E.o            \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)perm_metis_E_FC.o         \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)perm_metis_E_FCv.o        \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_metis_e_FC.o      \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_metis_e_FCv.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_metis_e_FCs.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_metis_e_FCvs.o    \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)perm_metis_N.o            \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)perm_metis_N_FC.o         \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)perm_metis_N_FCv.o        \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_metis_n_FC.o      \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_metis_n_FCv.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_metis_n_FCs.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_metis_n_FCvs.o    \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permamd.o                 \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)perm_amd_FC.o             \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)perm_amd_FCv.o            \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_amd_FC.o          \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_amd_FCv.o         \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_amd_FCs.o         \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_amd_FCvs.o        \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permmwmrcm.o              \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permmwmmmd.o              \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permmwmamf.o              \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permmwm_metis_e.o         \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permmwm_metis_n.o         \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permmwmamd.o              \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermmwm_metis_n.o      \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermmwm_metis_ns.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_metis_n_sp.o  \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_metis_n_sps.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_metis_n_FC.o  \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_metis_n_FCv.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_metis_n_FCs.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_metis_n_FCvs.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermmwm_metis_e.o      \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermmwm_metis_es.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_metis_e_sp.o  \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_metis_e_sps.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_metis_e_FC.o  \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_metis_e_FCv.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_metis_e_FCs.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_metis_e_FCvs.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermmwm_amd.o          \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermmwm_amds.o         \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_amd_sp.o      \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_amd_sps.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_amd_FC.o      \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_amd_FCv.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_amd_FCs.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_amd_FCvs.o    \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermmwm_mmd.o          \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermmwm_mmds.o         \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_mmd_sp.o      \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_mmd_sps.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_mmd_FC.o      \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_mmd_FCv.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_mmd_FCs.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_mmd_FCvs.o    \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermmwm_amf.o          \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermmwm_amfs.o         \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_amf_sp.o      \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_amf_sps.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_amf_FC.o      \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_amf_FCv.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_amf_FCs.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_amf_FCvs.o    \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermmwm_rcm.o          \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermmwm_rcms.o         \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_rcm_sp.o      \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_rcm_sps.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_rcm_FC.o      \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_rcm_FCv.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_rcm_FCs.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_mwm_rcm_FCvs.o    \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permMC64rcm.o             \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permMC64mmd.o             \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permMC64amf.o             \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permMC64_metis_e.o        \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permMC64_metis_n.o        \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)permMC64amd.o             \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermMC64rcm.o          \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermMC64rcms.o         \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_rcm_sp.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_rcm_sps.o    \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_rcm_FC.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_rcm_FCv.o    \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_rcm_FCs.o    \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_rcm_FCvs.o   \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermMC64_metis_n.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermMC64_metis_ns.o    \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_metis_n_sp.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_metis_n_sps.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_metis_n_FC.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_metis_n_FCv.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_metis_n_FCs.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_metis_n_FCvs.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermMC64_metis_e.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermMC64_metis_es.o    \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_metis_e_sp.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_metis_e_sps.o    \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_metis_e_FC.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_metis_e_FCv.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_metis_e_FCs.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_metis_e_FCvs.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermMC64amd.o          \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermMC64amds.o         \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_amd_sp.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_amd_sps.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_amd_FC.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_amd_FCv.o    \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_amd_FCs.o    \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_amd_FCvs.o   \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermMC64mmd.o          \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermMC64mmds.o         \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_mmd_sp.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_mmd_sps.o    \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_mmd_FC.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_mmd_FCv.o    \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_mmd_FCs.o    \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_mmd_FCvs.o   \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermMC64amf.o          \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)sympermMC64amfs.o         \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_amf_sp.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_amf_sps.o    \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_amf_FC.o     \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_amf_FCv.o    \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_amf_FCs.o    \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_MC64_amf_FCvs.o   \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symwms.o                  \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symwm.o                   \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)swm.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)perm_matching_metis_n.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_metis_n.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_metis_ns.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_metis_n_sp.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_metis_n_sps.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_metis_n_FC.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_metis_n_FCs.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_metis_n_FCv.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_metis_n_FCvs.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_metis_n_sp.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_metis_n_sps.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)perm_matching_metis_e.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_metis_e.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_metis_es.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_metis_e_sp.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_metis_e_sps.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_metis_e_FC.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_metis_e_FCs.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_metis_e_FCv.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_metis_e_FCvs.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)perm_matching_amf.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_amf.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_amfs.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_amf_sp.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_amf_sps.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_amf_FC.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_amf_FCs.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_amf_FCv.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_amf_FCvs.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)perm_matching_amd.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_amd.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_amds.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_amd_sp.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_amd_sps.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_amd_FC.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_amd_FCs.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_amd_FCv.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_amd_FCvs.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)perm_matching_mmd.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_mmd.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_mmds.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_mmd_sp.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_mmd_sps.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_mmd_FC.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_mmd_FCs.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_mmd_FCv.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_mmd_FCvs.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)perm_matching_rcm.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_rcm.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_rcms.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_rcm_sp.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_rcm_sps.o\
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_rcm_FC.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_rcm_FCs.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_rcm_FCv.o \
                $(DIRORDERINGS)$(FLOAT)/$(FLOAT)symperm_matching_rcm_FCvs.o

#                $(DIRORDERINGS)$(FLOAT)/permMWA_rcm.o             \
#                $(DIRORDERINGS)$(FLOAT)/permMWA_mmd.o             \
#                $(DIRORDERINGS)$(FLOAT)/permMWA_amf.o             \
#                $(DIRORDERINGS)$(FLOAT)/permMWA_amd.o             \
#                $(DIRORDERINGS)$(FLOAT)/permMWA_metis_e.o         \
#                $(DIRORDERINGS)$(FLOAT)/permMWA_metis_n.o         \

OBJECTORDERINGS=

ORDERINGS=      $(LIBORDERINGS)  $(OBJECTORDERINGS)



LIBSOLVERS=   $(DIRSOLVERS)$(FLOAT)/$(FLOAT)AMGsolver.o       \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)AMGsolver_real_prec.o \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)AMGsolver_real_sym_prec.o \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)AMGsolver_real_spd_prec.o \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)AMGsolver_sym_prec.o \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)AMGsolver_sym_precs.o \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)AMGsolver_spd_prec.o \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)AMGinit.o         \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)AMGgetparams.o    \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)AMGsetparams.o    \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)AMGsetupparameters.o \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)spdAMGsetupparameters.o \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)symAMGsetupparameters.o \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)symAMGsetupparameterss.o \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)spdAMGinit.o      \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)spdAMGgetparams.o \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)spdAMGsetparams.o \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)symAMGinit.o      \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)symAMGgetparams.o \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)symAMGsetparams.o \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)symAMGinits.o     \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)symAMGgetparamss.o\
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)symAMGsetparamss.o\
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)spdAMGsolver.o    \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)spdAMGsolver_real_prec.o    \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)symAMGsolver.o    \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)symAMGsolver_real_prec.o    \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)symAMGsolver_real_spd_prec.o    \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)symAMGsolver_spd_prec.o    \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)symAMGsolvers.o  \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)symAMGsolver_real_precs.o    \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)symAMGsolver_real_spd_precs.o    \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)symAMGdeletes.o  \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)symAMGdelete.o   \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)spdAMGdelete.o   \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)AMGdelete.o \
              $(DIRSOLVERS)$(FLOAT)/$(FLOAT)symspd.o 

OBJECTSOLVERS=


SOLVERS=      $(LIBSOLVERS) $(OBJECTSOLVERS)



LIBSPARSKIT=   $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)gmres.o        \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)pcg.o          \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)fpcg.o         \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)readmtc.o      \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)writemtc.o     \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)lusol.o        \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)lutsol.o       \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)lulsol.o       \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)lutlsol.o      \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)luusol.o       \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)lutusol.o      \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)ludlsol.o      \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)lutdlsol.o     \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)ludusol.o      \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)lutdusol.o     \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)ilutp2.o       \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)bisinit.o      \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)brkdn.o        \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)mgsro.o        \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)stopbis.o      \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)tidycg.o       \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)qsplit.o       \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)qsplit2.o      \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)qsplit3.o      \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)scale.o        \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)spdscale.o     \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)symscale.o     \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)csrcsc.o       \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)fgmres.o       \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)ilut2.o        \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)bcg.o          \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)readvectors.o  \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)writevectors.o \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)sbcg.o         \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)sbcgs.o        \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)sqmr.o         \
               $(DIRSPARSKIT)$(FLOAT)/$(FLOAT)sqmrs.o 

OBJECTSPARSKIT=

SPARSKIT=      $(LIBSPARSKIT) $(OBJECTSPARSKIT)



LIBTOOLS=   $(DIRTOOLS)$(FLOAT)/$(FLOAT)Malloc.o         \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)Realloc.o        \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)CSRmatvec.o      \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)CSRmattvec.o     \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)CSRmathvec.o     \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)geteps.o         \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)globals.o        \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)etree.o          \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)etree_postorder.o \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)qsort.o          \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)qsort2.o         \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)qqsort.o         \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)qqsorti.o        \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)qqsort2.o        \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)qqsorts.o        \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)qqsorts2.o       \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)swapm.o          \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)swapj.o          \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)CSRSparTran.o    \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)CSRSetupGraph.o  \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)CSRSetupGraph_epsilon.o\
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)CSRSetupGraph_epsilon_sp.o\
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)clear.o          \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)symmatvec.o      \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)skewmatvec.o     \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)cpermC.o         \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)rpermR.o         \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)csymmatvec.o     \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)cc_etimes.o      \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)GNLSYM.o         \
            $(DIRTOOLS)$(FLOAT)/$(FLOAT)GNLSYMs.o        

#            $(DIRTOOLS)$(FLOAT)/mps_arms.o 

OBJECTTOOLS=

TOOLS=      $(LIBTOOLS) $(OBJECTTOOLS)


LIBMUMPS=   $(DIRMUMPS)$(FLOAT)/$(FLOAT)mumps_match.o\
            $(DIRMUMPS)$(FLOAT)/$(FLOAT)mumps_203.o\
            $(DIRMUMPS)$(FLOAT)/$(FLOAT)mumps_444.o\
            $(DIRMUMPS)$(FLOAT)/$(FLOAT)mumps_445.o\
            $(DIRMUMPS)$(FLOAT)/$(FLOAT)mumps_446.o\
            $(DIRMUMPS)$(FLOAT)/$(FLOAT)mumps_447.o\
            $(DIRMUMPS)$(FLOAT)/$(FLOAT)mumps_448.o\
            $(DIRMUMPS)$(FLOAT)/$(FLOAT)mumps_450.o\
            $(DIRMUMPS)$(FLOAT)/$(FLOAT)mumps_451.o\
            $(DIRMUMPS)$(FLOAT)/$(FLOAT)mumps_452.o\
            $(DIRMUMPS)$(FLOAT)/$(FLOAT)mumps_453.o\
            $(DIRMUMPS)$(FLOAT)/$(FLOAT)mumps_454.o\
            $(DIRMUMPS)$(FLOAT)/$(FLOAT)mumps_455.o\
            $(DIRMUMPS)$(FLOAT)/$(FLOAT)mumps_457.o\
            $(DIRMUMPS)$(FLOAT)/$(FLOAT)mumps_551.o\
            $(DIRMUMPS)$(FLOAT)/$(FLOAT)mumps_559.o\
            $(DIRMUMPS)$(FLOAT)/$(FLOAT)mumps_562.o\
            $(DIRMUMPS)$(FLOAT)/$(FLOAT)mumps_563.o\
            $(DIRMUMPS)$(FLOAT)/$(FLOAT)mumps_abort.o\
            $(DIRMUMPS)$(FLOAT)/$(FLOAT)mumps_metric2x2.o\
            $(DIRMUMPS)$(FLOAT)/$(FLOAT)mumps_updatescore.o\
            $(DIRMUMPS)$(FLOAT)/$(FLOAT)mumps_update_inverse.o


OBJECTMUMPS=

MUMPS=      $(LIBMUMPS) $(OBJECTMUMPS)


LIBMETIS=   $(DIRMETIS)coarsen.o 	\
            $(DIRMETIS)fm.o 		\
            $(DIRMETIS)initpart.o 	\
            $(DIRMETIS)match.o  	\
            $(DIRMETIS)ccgraph.o 	\
            $(DIRMETIS)memory.o 	\
            $(DIRMETIS)pmetis.o 	\
            $(DIRMETIS)pqueue.o 	\
            $(DIRMETIS)refine.o 	\
            $(DIRMETIS)util.o 		\
            $(DIRMETIS)timing.o 	\
            $(DIRMETIS)debug.o		\
            $(DIRMETIS)bucketsort.o 	\
            $(DIRMETIS)graph.o  	\
            $(DIRMETIS)stat.o 		\
            $(DIRMETIS)kmetis.o 	\
            $(DIRMETIS)kwayrefine.o	\
            $(DIRMETIS)kwayfm.o 	\
            $(DIRMETIS)balance.o 	\
            $(DIRMETIS)ometis.o 	\
            $(DIRMETIS)srefine.o 	\
            $(DIRMETIS)sfm.o 		\
            $(DIRMETIS)separator.o	\
            $(DIRMETIS)mincover.o 	\
            $(DIRMETIS)mmd.o 		\
            $(DIRMETIS)mesh.o 		\
            $(DIRMETIS)meshpart.o 	\
            $(DIRMETIS)frename.o 	\
            $(DIRMETIS)fortran.o	\
            $(DIRMETIS)myqsort.o 	\
            $(DIRMETIS)compress.o 	\
            $(DIRMETIS)parmetis.o 	\
            $(DIRMETIS)estmem.o	        \
            $(DIRMETIS)mpmetis.o 	\
            $(DIRMETIS)mcoarsen.o 	\
            $(DIRMETIS)mmatch.o 	\
            $(DIRMETIS)minitpart.o 	\
            $(DIRMETIS)mbalance.o	\
            $(DIRMETIS)mrefine.o 	\
            $(DIRMETIS)mutil.o 	        \
            $(DIRMETIS)mfm.o 		\
            $(DIRMETIS)mkmetis.o 	\
            $(DIRMETIS)mkwayrefine.o 	\
            $(DIRMETIS)mkwayfmh.o 	\
            $(DIRMETIS)mrefine2.o 	\
            $(DIRMETIS)minitpart2.o 	\
            $(DIRMETIS)mbalance2.o 	\
            $(DIRMETIS)mfm2.o		\
            $(DIRMETIS)kvmetis.o 	\
            $(DIRMETIS)kwayvolrefine.o  \
            $(DIRMETIS)kwayvolfm.o    	\
            $(DIRMETIS)subdomains.o
#LIBMETIS=   

OBJECTMETIS=

METIS=      $(LIBMETIS) $(OBJECTMETIS)



LIBAMD=   $(DIRAMD)$(AMD_SWITCH)amd_1.o              \
          $(DIRAMD)$(AMD_SWITCH)amd_2.o              \
          $(DIRAMD)$(AMD_SWITCH)amd_aat.o            \
          $(DIRAMD)$(AMD_SWITCH)amd_control.o        \
          $(DIRAMD)$(AMD_SWITCH)amd_defaults.o       \
          $(DIRAMD)$(AMD_SWITCH)amd_dump.o           \
          $(DIRAMD)$(AMD_SWITCH)amd_info.o           \
          $(DIRAMD)$(AMD_SWITCH)amd_order.o          \
          $(DIRAMD)$(AMD_SWITCH)amd_postorder.o      \
          $(DIRAMD)$(AMD_SWITCH)amd_post_tree.o      \
          $(DIRAMD)$(AMD_SWITCH)amd_preprocess.o     \
          $(DIRAMD)$(AMD_SWITCH)amd_valid.o     

#LIBAMD=   

OBJECTAMD=

AMD=      $(LIBAMD) $(OBJECTAMD)



ALLAMD=   $(DIRAMD)s_amd_1.o              \
          $(DIRAMD)s_amd_2.o              \
          $(DIRAMD)s_amd_aat.o            \
          $(DIRAMD)s_amd_control.o        \
          $(DIRAMD)s_amd_defaults.o       \
          $(DIRAMD)s_amd_dump.o           \
          $(DIRAMD)s_amd_info.o           \
          $(DIRAMD)s_amd_order.o          \
          $(DIRAMD)s_amd_postorder.o      \
          $(DIRAMD)s_amd_post_tree.o      \
          $(DIRAMD)s_amd_preprocess.o     \
          $(DIRAMD)s_amd_valid.o          \
          $(DIRAMD)l_amd_1.o              \
          $(DIRAMD)l_amd_2.o              \
          $(DIRAMD)l_amd_aat.o            \
          $(DIRAMD)l_amd_control.o        \
          $(DIRAMD)l_amd_defaults.o       \
          $(DIRAMD)l_amd_dump.o           \
          $(DIRAMD)l_amd_info.o           \
          $(DIRAMD)l_amd_order.o          \
          $(DIRAMD)l_amd_postorder.o      \
          $(DIRAMD)l_amd_post_tree.o      \
          $(DIRAMD)l_amd_preprocess.o     \
          $(DIRAMD)l_amd_valid.o     
#ALLAMD=   



SPARSPAK=$(DIRSPARSPAK)genqmd.o \
         $(DIRSPARSPAK)fndsep.o \
         $(DIRSPARSPAK)gennd.o  \
         $(DIRSPARSPAK)qmdrch.o \
         $(DIRSPARSPAK)qmdupd.o \
         $(DIRSPARSPAK)qmdmrg.o \
         $(DIRSPARSPAK)qmdqt.o  \
         $(DIRSPARSPAK)revrse.o \
         $(DIRSPARSPAK)genrcm.o \
         $(DIRSPARSPAK)rcm.o    \
         $(DIRSPARSPAK)rootls.o \
         $(DIRSPARSPAK)degree.o \
         $(DIRSPARSPAK)fnroot.o 
# the codes mentioned above are not distributed
# the definition is only done for convenience
#SPARSPAK=


BLAS=$(DIRBLAS)dnrm2.o  $(DIRBLAS)dznrm2.o $(DIRBLAS)dcopy.o  $(DIRBLAS)zcopy.o \
     $(DIRBLAS)daxpy.o  $(DIRBLAS)zaxpy.o  $(DIRBLAS)drotg.o  $(DIRBLAS)zrotg.o \
     $(DIRBLAS)drot.o   $(DIRBLAS)dcabs1.o $(DIRBLAS)idamax.o $(DIRBLAS)csscal.o\
     $(DIRBLAS)ddot.o   $(DIRBLAS)zdotc.o  $(DIRBLAS)zhpr.o   $(DIRBLAS)cgemv.o \
     $(DIRBLAS)zdotu.o  $(DIRBLAS)dscal.o  $(DIRBLAS)zscal.o  $(DIRBLAS)dspr.o  \
     $(DIRBLAS)snrm2.o  $(DIRBLAS)scnrm2.o $(DIRBLAS)scopy.o  $(DIRBLAS)ccopy.o \
     $(DIRBLAS)saxpy.o  $(DIRBLAS)caxpy.o  $(DIRBLAS)srotg.o  $(DIRBLAS)crotg.o \
     $(DIRBLAS)srot.o   $(DIRBLAS)isamax.o $(DIRBLAS)icamax.o $(DIRBLAS)izamax.o\
     $(DIRBLAS)sdot.o   $(DIRBLAS)cdotc.o  $(DIRBLAS)zgemv.o  $(DIRBLAS)dgemv.o \
     $(DIRBLAS)cdotu.o  $(DIRBLAS)sscal.o  $(DIRBLAS)cscal.o  $(DIRBLAS)zdscal.o\
     $(DIRBLAS)cswap.o  $(DIRBLAS)dswap.o  $(DIRBLAS)sswap.o  $(DIRBLAS)zswap.o \
     $(DIRBLAS)ctrsm.o  $(DIRBLAS)dtrsm.o  $(DIRBLAS)strsm.o  $(DIRBLAS)ztrsm.o \
     $(DIRBLAS)dger.o   $(DIRBLAS)sger.o   $(DIRBLAS)zgeru.o  $(DIRBLAS)cgeru.o \
     $(DIRBLAS)dtpsv.o  $(DIRBLAS)stpsv.o  $(DIRBLAS)ztpsv.o  $(DIRBLAS)ctpsv.o \
     $(DIRBLAS)dgemm.o  $(DIRBLAS)sgemm.o  $(DIRBLAS)zgemm.o  $(DIRBLAS)cgemm.o \
     $(DIRBLAS)dasum.o  $(DIRBLAS)sasum.o  $(DIRBLAS)dzasum.o $(DIRBLAS)scasum.o\
     $(DIRBLAS)lsame.o  $(DIRBLAS)sgemv.o  $(DIRBLAS)sspr.o   $(DIRBLAS)dtrsv.o \
     $(DIRBLAS)zher.o   $(DIRBLAS)dtrmm.o  $(DIRBLAS)dsyr.o   $(DIRBLAS)ssyr.o  \
     $(DIRBLAS)dsyrk.o  $(DIRBLAS)dsyr2.o  $(DIRBLAS)dsyr2.o  $(DIRBLAS)dsymv.o \
     $(DIRBLAS)dsyr2k.o $(DIRBLAS)dtrmv.o
#     $(DIRBLAS).o  $(DIRBLAS).o  $(DIRBLAS).o  $(DIRBLAS).o



BLASLIKE=$(DIRBLASLIKE)zrot.o $(DIRBLASLIKE)ddot2.o  $(DIRBLASLIKE)zdotc2.o   \
         $(DIRBLASLIKE)zdotu2.o \
         $(DIRBLASLIKE)crot.o $(DIRBLASLIKE)sdot2.o  $(DIRBLASLIKE)cdotc2.o   \
         $(DIRBLASLIKE)cdotu2.o \
         $(DIRBLASLIKE)sprivatesptrs.o \
         $(DIRBLASLIKE)dprivatesptrs.o \
         $(DIRBLASLIKE)cprivatehptrs.o \
         $(DIRBLASLIKE)zprivatehptrs.o 


LAPACK=$(DIRLAPACK)dlamch.o  $(DIRLAPACK)slamch.o  \
       $(DIRLAPACK)xerbla.o  $(DIRLAPACK)lsame.o   \
       $(DIRLAPACK)dsptrf.o  $(DIRLAPACK)dsptrs.o  \
       $(DIRLAPACK)ssptrf.o  $(DIRLAPACK)ssptrs.o  \
       $(DIRLAPACK)csptrf.o  $(DIRLAPACK)csptrs.o  \
       $(DIRLAPACK)zsptrf.o  $(DIRLAPACK)zsptrs.o  \
       $(DIRLAPACK)chptrf.o  $(DIRLAPACK)chptrs.o  \
       $(DIRLAPACK)zhptrf.o  $(DIRLAPACK)zhptrs.o  \
       $(DIRLAPACK)cspr.o    $(DIRLAPACK)zspr.o    \
       $(DIRLAPACK)clacgv.o  $(DIRLAPACK)zlacgv.o  \
       $(DIRLAPACK)slapy2.o  $(DIRLAPACK)dlapy2.o  \
       $(DIRLAPACK)spptrf.o  $(DIRLAPACK)spptrs.o  \
       $(DIRLAPACK)dpptrf.o  $(DIRLAPACK)dpptrs.o  \
       $(DIRLAPACK)cpptrf.o  $(DIRLAPACK)cpptrs.o  \
       $(DIRLAPACK)zpptrf.o  $(DIRLAPACK)zpptrs.o  \
       $(DIRLAPACK)sgetrf.o  $(DIRLAPACK)sgetrs.o  \
       $(DIRLAPACK)dgetrf.o  $(DIRLAPACK)dgetrs.o  \
       $(DIRLAPACK)cgetrf.o  $(DIRLAPACK)cgetrs.o  \
       $(DIRLAPACK)zgetrf.o  $(DIRLAPACK)zgetrs.o  \
       $(DIRLAPACK)dgetf2.o  $(DIRLAPACK)sgetf2.o  \
       $(DIRLAPACK)zgetf2.o  $(DIRLAPACK)cgetf2.o  \
       $(DIRLAPACK)dlaswp.o  $(DIRLAPACK)slaswp.o  \
       $(DIRLAPACK)zlaswp.o  $(DIRLAPACK)claswp.o  \
       $(DIRLAPACK)ilaenv.o  $(DIRLAPACK)ieeeck.o  \
       $(DIRLAPACK)dgeqrf.o  $(DIRLAPACK)dgetri.o  \
       $(DIRLAPACK)dormqr.o  $(DIRLAPACK)dormqr.o  \
       $(DIRLAPACK)dggev.o   $(DIRLAPACK)dsytrf.o  \
       $(DIRLAPACK)dsytri.o  $(DIRLAPACK)dpotrf.o  \
       $(DIRLAPACK)dlartg.o  $(DIRLAPACK)zpotf2.o  \
       $(DIRLAPACK)dgeev.o   $(DIRLAPACK)dsyev.o   \
       $(DIRLAPACK)zsyr.o    $(DIRLAPACK)dsyevr.o  \
       $(DIRLAPACK)zlaev2.o  $(DIRLAPACK)dpotf2.o  \
       $(DIRLAPACK)dlansy.o  $(DIRLAPACK)dstebz.o  \
       $(DIRLAPACK)dormtr.o  $(DIRLAPACK)dstein.o  \
       $(DIRLAPACK)dstegr.o  $(DIRLAPACK)dsterf.o  \
       $(DIRLAPACK)dlascl.o  $(DIRLAPACK)dsytrd.o  \
       $(DIRLAPACK)dlanst.o  $(DIRLAPACK)dlaset.o  \
       $(DIRLAPACK)dlae2.o   $(DIRLAPACK)dlasrt.o  \
       $(DIRLAPACK)dlatrd.o  $(DIRLAPACK)dsytd2.o  \
       $(DIRLAPACK)dlatrd.o  $(DIRLAPACK)dlarre.o  \
       $(DIRLAPACK)dlarrv.o  $(DIRLAPACK)dlagts.o  \
       $(DIRLAPACK)dlarrf.o  $(DIRLAPACK)dlar1v.o  \
       $(DIRLAPACK)dlarrb.o  $(DIRLAPACK)dlasq2.o  \
       $(DIRLAPACK)dlarfg.o  $(DIRLAPACK)dlarnv.o  \
       $(DIRLAPACK)dlassq.o  $(DIRLAPACK)dlasq3.o  \
       $(DIRLAPACK)dlagtf.o  $(DIRLAPACK)dlaebz.o  \
       $(DIRLAPACK)dgebak.o  $(DIRLAPACK)dtrevc.o  \
       $(DIRLAPACK)dlacpy.o  $(DIRLAPACK)dlange.o  \
       $(DIRLAPACK)dlabad.o  $(DIRLAPACK)dlaln2.o  \
       $(DIRLAPACK)dlasq4.o  $(DIRLAPACK)dlasq5.o  \
       $(DIRLAPACK)dlasq6.o  $(DIRLAPACK)dlaev2.o  \
       $(DIRLAPACK)dormql.o  $(DIRLAPACK)dorghr.o  \
       $(DIRLAPACK)dlaruv.o  $(DIRLAPACK)dorm2l.o  \
       $(DIRLAPACK)dlarft.o  $(DIRLAPACK)dlarfb.o  \
       $(DIRLAPACK)dladiv.o  $(DIRLAPACK)dorgtr.o  \
       $(DIRLAPACK)dsteqr.o  $(DIRLAPACK)dsytf2.o  \
       $(DIRLAPACK)dlasyf.o  $(DIRLAPACK)dgebal.o  \
       $(DIRLAPACK)dgehrd.o  $(DIRLAPACK)dhseqr.o  \
       $(DIRLAPACK)dorgqr.o  $(DIRLAPACK)dorg2r.o  \
       $(DIRLAPACK)dlahqr.o  $(DIRLAPACK)dlanhs.o  \
       $(DIRLAPACK)dlarfx.o  $(DIRLAPACK)dlasr.o   \
       $(DIRLAPACK)dgehd2.o  $(DIRLAPACK)dlahrd.o  \
       $(DIRLAPACK)dggbak.o  $(DIRLAPACK)dtgevc.o  \
       $(DIRLAPACK)dhgeqz.o  $(DIRLAPACK)dggbal.o  \
       $(DIRLAPACK)dlapy3.o  $(DIRLAPACK)dlag2.o   \
       $(DIRLAPACK)dlasv2.o  $(DIRLAPACK)dlarf.o   \
       $(DIRLAPACK)dorgql.o  $(DIRLAPACK)dgghrd.o  \
       $(DIRLAPACK)dorm2r.o  $(DIRLAPACK)dgeqr2.o  \
       $(DIRLAPACK)dtrtri.o  $(DIRLAPACK)dlanv2.o  \
       $(DIRLAPACK)dorg2l.o  $(DIRLAPACK)dtrti2.o  \
       $(DIRLAPACK)dgels.o   $(DIRLAPACK)dlartv.o  \
       $(DIRLAPACK)dgeqlf.o  $(DIRLAPACK)dormlq.o  \
       $(DIRLAPACK)dgelqf.o  $(DIRLAPACK)dorml2.o  \
       $(DIRLAPACK)dgelq2.o #  $(DIRLAPACK).o
#       $(DIRLAPACK).o  $(DIRLAPACK).o



LIBFORTRAN=$(DIRFORTRAN)$(FLOAT)/$(FLOAT)ilupackinit.o          $(DIRFORTRAN)$(FLOAT)/$(FLOAT)spdilupackinit.o \
           $(DIRFORTRAN)$(FLOAT)/$(FLOAT)symilupackinit.o       $(DIRFORTRAN)$(FLOAT)/$(FLOAT)symilupackinits.o \
           $(DIRFORTRAN)$(FLOAT)/$(FLOAT)ilupackfactor.o        $(DIRFORTRAN)$(FLOAT)/$(FLOAT)spdilupackfactor.o \
           $(DIRFORTRAN)$(FLOAT)/$(FLOAT)symilupackfactor.o     $(DIRFORTRAN)$(FLOAT)/$(FLOAT)symilupackfactors.o \
           $(DIRFORTRAN)$(FLOAT)/$(FLOAT)ilupacksolver.o        $(DIRFORTRAN)$(FLOAT)/$(FLOAT)spdilupacksolver.o \
           $(DIRFORTRAN)$(FLOAT)/$(FLOAT)symilupacksolver.o     $(DIRFORTRAN)$(FLOAT)/$(FLOAT)symilupacksolvers.o \
           $(DIRFORTRAN)$(FLOAT)/$(FLOAT)ilupacksol.o           $(DIRFORTRAN)$(FLOAT)/$(FLOAT)spdilupacksol.o \
           $(DIRFORTRAN)$(FLOAT)/$(FLOAT)symilupacksol.o        $(DIRFORTRAN)$(FLOAT)/$(FLOAT)symilupacksols.o \
           $(DIRFORTRAN)$(FLOAT)/$(FLOAT)ilupackdelete.o        $(DIRFORTRAN)$(FLOAT)/$(FLOAT)spdilupackdelete.o \
           $(DIRFORTRAN)$(FLOAT)/$(FLOAT)symilupackdelete.o     $(DIRFORTRAN)$(FLOAT)/$(FLOAT)symilupackdeletes.o\
           $(DIRFORTRAN)$(FLOAT)/$(FLOAT)spdilupackinfo.o       $(DIRFORTRAN)$(FLOAT)/$(FLOAT)ilupackinfo.o \
           $(DIRFORTRAN)$(FLOAT)/$(FLOAT)symilupackinfo.o       $(DIRFORTRAN)$(FLOAT)/$(FLOAT)symilupackinfos.o\
           $(DIRFORTRAN)$(FLOAT)/$(FLOAT)spdilupacknnz.o        $(DIRFORTRAN)$(FLOAT)/$(FLOAT)ilupacknnz.o\
           $(DIRFORTRAN)$(FLOAT)/$(FLOAT)symilupacknnz.o        $(DIRFORTRAN)$(FLOAT)/$(FLOAT)symilupacknnzs.o\
           $(DIRFORTRAN)$(FLOAT)/$(FLOAT)symspdilupackconvert.o $(DIRFORTRAN)$(FLOAT)/$(FLOAT)unscale.o   

OBJECTFORTRAN=

FORTRAN=      $(LIBFORTRAN) $(OBJECTFORTRAN)



SAMPLES=$(FLOAT)/$(FLOAT)$(MAIN).o



STARTDIR=$(PWD)

MYSTARTDIR=$(STARTDIR)







# where are the headers
INCDIR=$(MYSTARTDIR)/include

# where are the libraries
LIBDIR=$(MYSTARTDIR)/lib/$(MYPLATFORM)


.SUFFIXES: .c .f .F .mod .o .a
.DEFAULT: main

$(FLOAT)/$(FLOAT)%.o: %.c
	$(CC)  $(CCFLAGS)  -I$(INCDIR)  $(MATCHING) $(FORTRANNAMES) $(ARITHMETIC) $(LONGINTEGER) -c -o $@ $<


$(FLOAT)/$(FLOAT)%.o: %.F
	$(FCOMPILE)





