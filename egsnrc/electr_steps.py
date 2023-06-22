from typing import Tuple
from egsnrc.randoms import randomset
from egsnrc.params import *
from egsnrc.commons import *
from egsnrc.constants import *

from egsnrc.calcfuncs import (
    calc_tstep_from_demfp, calculate_xi, compute_drange,
    compute_eloss_g
)
from egsnrc.angles import uphi

from math import log

import logging
logger = logging.getLogger("egsnrc")


# EMPTY CALLBACKS ----
calculate_elastic_scattering_mfp = None
call_howfar_in_electr = None
check_negative_ustep = None

de_fluctuation = None
"""Allows the user to change the ionization loss.
(Provides a user hook for Landau/Vavilov processes)
"""

add_work_em_field = None
add_work_em_field2 = None
electron_region_change = None
emfield_initiate_set_tustep = None
evaluate_sig0 = None
evaluate_sigf = None
range_discard = None
scale_sig0 = None
select_electron_mfp = None
set_angles_em_field = None
set_skindepth = None

set_tvstep = None
set_tvstep_em_field = None
set_ustep = None
set_ustep_em_field = None
update_demfp = None
user_controls_tstep_recursion = None
user_range_discard = None


ierust = 0 # To count negative ustep's


# @profile  # for line_profiler
def tstep_ustep(
    lelec:int, medium: int, irl:int,
    eie: float, peie:float,
    hownear, howfar, ausgab
) -> Tuple[int, int, float, float]:
    """Follow one particle until is cut

    Returns
    -------
    int
        Flag for reason returning, user_electron_discard or ecut_discard,
        or None for no reason (rdict > sigratio in tstep check)
    """
    # Important definitions:
    # tstep  = total pathlength to the next discrete interaction
    # vacdst = infinity (actually 10^8)
    # tustep = total pathlength of the electron step
    # ustep  = projected transport distance in the
    #          direction of motion at the start of the step

    global costhe, sinthe, ierust

    next_tstep = False  # flag for jumping out of ustep to top of tstep loop

    # Define output variables from egsfortran functions, to avoid error message
    # For f2py, need 'scalar' which are rank-0 arrays
    uscat = numpy.array(0, dtype=numpy.float64)  # just need to define something
    vscat = numpy.array(0, dtype=numpy.float64)
    wscat = numpy.array(0, dtype=numpy.float64)
    xtrans = numpy.array(0, dtype=numpy.float64)
    ytrans = numpy.array(0, dtype=numpy.float64)
    ztrans = numpy.array(0, dtype=numpy.float64)

    # Define local variables for increased lookup speed
    call_tranausb = True if iausfl[TRANAUSB] !=0 else False
    call_tranausa = True if iausfl[TRANAUSA] !=0 else False
    smallest_electron_mfp = EPSEMFP

    qel = 0 if lelec == -1 else 1  # for array indexing
    medium_m1 = medium - 1  # ** 0-based
    irl_m1 = irl - 1
    np_m1 = np - 1

    while True:  # :TSTEP: LOOP
        # Go through this loop each time we recompute distance to an interaction
        # ******* trying to save evaluation of range.
        # do_range = True # compute the range in $ COMPUTE-RANGE below
        # ********/
        compute_tstep = True  # MFP resampled => calculate distance to the
                                # interaction in the USTEP loop
        epcont.eke = eie - rm  # moved here so that kinetic energy will be known
                        # to user even for a vacuum step, IK January 2000
        # logger.debug(f'New TSTEP, eke={eke} ===================')
        if medium != 0:
            # Not vacuum. Must sample to see how far to next interaction.
            # --- Inline replace: $ SELECT_ELECTRON_MFP; -----
            if select_electron_mfp:
                select_electron_mfp()
            else:
                rnne1 = randomset()
                # logger.info(f'random rnne1={rnne1:0.8e}')
                if rnne1 == 0.0:
                    rnne1=1.e-30
                demfp = max(-log(rnne1), smallest_electron_mfp)
            # End inline replace: $ SELECT_ELECTRON_MFP; ----
            # demfp = differential electron mean free path
            # logger.debug(f"** Random: {rnne1} **")
            # logger.debug(f'non-vac, demfp={demfp}')
            epcont.elke = log(eke)
            # (eke = kinetic energy, rm = rest mass, all in units of MeV)

            # Prepare to approximate cross section
            # $ SET INTERVAL elke,eke;
            lelke = int(eke1[medium_m1]*elke + eke0[medium_m1])
            lelke_m1 = lelke - 1

            # --- Inline replace: $ EVALUATE_SIG0; -----
            if evaluate_sig0:
                evaluate_sig0()
            else:
                if sig_ismonotone[qel,medium_m1]:
                    # --- Inline replace: $ EVALUATE_SIGF; -----
                    if evaluate_sigf:
                        evaluate_sigf()
                    else:
                        if lelec < 0:
                            # EVALUATE sigf USING esig(elke)
                            sigf = esig1[lelke_m1, medium_m1]*elke+ esig0[lelke_m1, medium_m1]
                            # EVALUATE dedx0 USING ededx(elke)
                            dedx0 = ededx1[lelke_m1, medium_m1]*elke+ ededx0[lelke_m1, medium_m1]
                            sigf = sigf/dedx0
                        else:
                            # EVALUATE sigf USING psig(elke)
                            sigf = psig1[lelke_m1, medium_m1]*elke+ psig0[lelke_m1, medium_m1]
                            # EVALUATE dedx0 USING pdedx(elke)
                            dedx0 = pdedx1[lelke_m1, medium_m1]*elke+ pdedx0[lelke_m1, medium_m1]
                            sigf = sigf/dedx0
                    # End inline replace: $ EVALUATE_SIGF; ----
                    sig0 = sigf
                else:
                    if lelec < 0:
                        sig0 = esig_e[medium_m1]
                    else:
                        sig0 = psig_e[medium_m1]

            # End inline replace: $ EVALUATE_SIG0; ----
            # The fix up of the fictitious method uses cross section per
            # energy loss. Therefore, demfp/sig is sub-threshold energy loss
            # until the next discrete interaction occures (see below)
            # As this quantity is a single constant for a material,
            # $ SET INTERVAL is not necessary at this point. However, to not
            # completely alter the logic of the TSTEP and USTEP loops,
            # this is left for now
            # logger.debug(f'sig0={sig0}')
        # end non-vacuum test
        # ----------------------------------------------------------------
        while True:  # :USTEP: LOOP
            # Here for each check with user geometry.
            # Compute size of maximum acceptable step, which is limited
            # by multiple scattering or other approximations.
            # logger.debug("New USTEP -----------------------------")
            if medium == 0:  # vacuum
                epcont.tstep = vacdst
                epcont.tustep = tstep
                epcont.ustep = tstep
                callhowfar = True # Always call HOWFAR for vacuum steps!

                # (Important definitions:
                #  tstep  = total pathlength to the next discrete interaction
                #  vacdst = infinity (actually 10^8)
                #  tustep = total pathlength of the electron step
                #  ustep  = projected transport distance in the
                #           direction of motion at the start of the step
                #  Note that tustep and ustep are modified below.
                #  The above provide defaults.)

                #  EM field step size restriction in vacuum

                epcont.ustep = tustep
            else:
                # non-vacuum
                epcont.rhof = rhor[irl_m1] / rho[medium_m1] # density ratio scaling template
                            # EGS allows the density to vary
                            # continuously (user option)

                # --- Inline replace: $ SCALE_SIG0; -----
                if scale_sig0:
                    scale_sig0()
                else:
                    sig = sig0
                # End inline replace: $ SCALE_SIG0; ----
                if sig <= 0:
                    # This can happen if the threshold for brems,
                    # (ap + rm), is greater than ae.  Moller threshold is
                    # 2*ae - rm. If sig is zero, we are below the
                    # thresholds for both bremsstrahlung and Moller.
                    # In this case we will just lose energy by
                    # ionization loss until we go below cut-off. Do not
                    # assume range is available, so just ask for step
                    # same as vacuum.  Electron transport will reduce
                    # into little steps.
                    # (Note: ae is the lower threshold for creation of a
                    #        secondary Moller electron, ap is the lower
                    #        threshold for creation of a brem.)
                    epcont.tstep = vacdst
                    sig0 = 1.E-15
                    # logger.debug("sig <= 0")
                else:
                    # --- Inline replace: $ CALCULATE_TSTEP_FROM_DEMFP; -----
                    if compute_tstep:
                        total_de = demfp / sig
                        total_tstep = calc_tstep_from_demfp(
                            qel,lelec,medium,lelke,demfp,sig,eke,elke,total_de
                        )
                        compute_tstep = False

                    epcont.tstep = total_tstep / rhof #  non-default density scaling
                    # End inline replace: $ CALCULATE_TSTEP_FROM_DEMFP; ----
                # end sig if-else
                # logger.debug(f'ustep non-vac, calc tstep={tstep}')
                # calculate stopping power
                if lelec < 0:
                    # EVALUATE dedx0 USING ededx(elke)] # e-
                    dedx0 = ededx1[lelke_m1, medium_m1]*elke+ ededx0[lelke_m1, medium_m1]
                else:
                    # EVALUATE dedx0 USING pdedx(elke) # e+
                    dedx0 = pdedx1[lelke_m1, medium_m1]*elke+ pdedx0[lelke_m1, medium_m1]
                dedx = rhof*dedx0
                # logger.debug(f'dedx={dedx}')
                # Determine maximum step-size (Formerly $ SET-TUSTEP)
                # EVALUATE tmxs USING tmxs(elke)
                tmxs = tmxs1[lelke_m1, medium_m1]*elke+ tmxs0[lelke_m1, medium_m1]
                tmxs = tmxs/rhof
                # logger.debug(f'tmxs={tmxs}')
                # Compute the range to E_min[medium_m1] (e_min is the first
                # energy in the table). Do not go more than range.
                # Don't replace this macro and don't override range_, because
                # the energy loss evaluation below relies on the accurate
                # (and self-consistent) evaluation of range_!
                # --- Inline replace: $ COMPUTE_RANGE; -----
                #         ===============
                ekei = e_array[lelke-1, medium_m1]  # ** 0-based
                elkei = (lelke - eke0[medium_m1]) / eke1[medium_m1]
                range_ = compute_drange(
                    lelec, medium, eke, ekei, lelke, elke, elkei
                )
                range_ = (range_ + range_ep[qel, lelke_m1, medium_m1]) / rhof
                # End inline replace: $ COMPUTE_RANGE; ----
                # logger.debug(f'range_={range_}')
                # The RANDOMIZE-TUSTEP option as coded by AFB forced the
                # electrons to approach discrete events (Moller,brems etc.)
                # only in a single scattering mode => waste of CPU time.
                # Moved here and changed by IK Oct 22 1997
                if RANDOMIZE_TUSTEP:
                    rnnotu = randomset()
                    # logger.info(f'random rnnotu={rnnotu:0.8e}')
                    tmxs = rnnotu*min(tmxs,smaxir[irl_m1])
                else:
                    tmxs = min(tmxs,smaxir[irl_m1])

                epcont.tustep = min(tstep,tmxs,range_)
                # logger.debug(f'tustep={tustep}')
                # optional tustep restriction in EM field


                # --- Inline replace: $ CALL_HOWNEAR(tperp); -----
                tperp = hownear(x[np_m1], y[np_m1], z[np_m1], irl)
                # End inline replace: $ CALL_HOWNEAR(tperp); ----
                # logger.debug(f'tperp={tperp} = hownear({x[np_m1]}, {y[np_m1]}, {z[np_m1]}, {irl})')
                dnear[np_m1] = tperp

                # --- Inline replace: $ RANGE_DISCARD; -----
                # optional regional range rejection for
                # particles below e_max_rr if i_do_rr set
                if range_discard:
                    range_discard()
                elif i_do_rr[irl_m1] == 1 and e[np_m1] < e_max_rr[irl_m1]:
                    if tperp >= range_:
                        # particle cannot escape local region
                        # set idisc 1 for electrons, 99 for positrons
                        epcont.idisc = 50 + 49*iq[np_m1]
                        return USER_ELECTRON_DISCARD, lelke, eie, peie
                # End inline replace: $ RANGE_DISCARD; ----

                if user_range_discard:
                    user_range_discard()

                # --- Inline replace: $ SET_SKINDEPTH(eke,elke); -----
                # This macro sets the minimum step size for a condensed
                # history (CH) step. When the exact BCA is used, the minimum
                # CH step is determined by efficiency considerations only
                # At about 3 elastic MFP's single scattering becomes more
                # efficient than CH and so the algorithm switches off CH
                # If one of the various inexact BCA's is invoked, this macro
                # provides a simple way to include more sophisticated
                # decisions about the maximum acceptable approximated CH step
                if set_skindepth:
                    skindepth = set_skindepth(eke, elke)
                else:
                    # --- Inline replace: $ CALCULATE_ELASTIC_SCATTERING_MFP(ssmfp,eke,elke); -----
                    if calculate_elastic_scattering_mfp:
                        ssmfp = calculate_elastic_scattering_mfp(ssmfp, eke, elke)
                    else:
                        blccl = rhof*blcc[medium_m1]
                        xccl  = rhof*xcc[medium_m1]
                        p2 = eke*(eke+rmt2); beta2 = p2/(p2 + rmsq)
                        if spin_effects :
                            if lelec < 0:
                                # EVALUATE etap USING etae_ms(elke)
                                etap = etae_ms1[lelke_m1, medium_m1]*elke+ etae_ms0[lelke_m1, medium_m1]
                            else:
                                # EVALUATE etap USING etap_ms(elke)
                                etap = etap_ms1[lelke_m1, medium_m1]*elke+ etap_ms0[lelke_m1, medium_m1]
                            # EVALUATE ms_corr USING blcce(elke)
                            ms_corr = blcce1[lelke_m1, medium_m1]*elke+ blcce0[lelke_m1, medium_m1]
                            blccl = blccl/etap/(1+0.25*etap*xccl/blccl/p2)*ms_corr
                        ssmfp=beta2/blccl
                    # End inline replace: $ CALCULATE_ELASTIC_SCATTERING_MFP(ssmfp,eke,elke); ----
                    skindepth = skindepth_for_bca*ssmfp
                # End inline replace: $ SET_SKINDEPTH(eke,elke); ----

                epcont.tustep = min(tustep,max(tperp,skindepth))

                if emfield_initiate_set_tustep:
                    emfield_initiate_set_tustep()

                # The transport logic below is determined by the logical
                # variables callhhowfar, domultiple and dosingle
                #
                # There are the following possibilities:
                #
                #    callhowfar = False  This indicates that the
                #    ====================  intended step is shorter than tperp
                #                          independent of BCA used
                #   - domultiple = False dosingle = False and
                #                          callmsdist = True
                #        ==> everything has been done in msdist
                #   - domultiple = True and dosingle = False
                #        ==> should happen only if exact_bca  is False
                #            indicates that MS remains to be done
                #   - domultiple = False and dosingle = True
                #        ==> should happen only if exact_bca  is True
                #            sampled distance to a single scattering event is
                #            shorter than tperp ==> do single scattering at the
                #            end of the step
                #   - domultiple = True and dosingle = True
                #        ==> error condition, something with the logic is wrong!
                #
                #    callhowfar = True This indicates that the intended step
                #    =================== is longer than tperp and forces a
                #                        call to hawfar which returns the
                #                        straight line distance to the boundary
                #                        in the initial direction of motion
                #                        (via a modification of ustep)
                #   - domultiple = False and dosingle = False
                #        ==> should happen only of exact_bca=True
                #            simply put the particle on the boundary
                #   - domultiple = False and dosingle = True
                #        ==> should happen only of exact_bca=True
                #            single elastic scattering has to be done
                #   - domultiple = True and dosingle = False
                #        ==> should happen only of exact_bca=False
                #            indicates that MS remains to be done
                #   - domultiple = True and dosingle = True
                #        ==> error condition, something with the logic is wrong!

                # IF(tustep <= tperp and tustep > skindepth)
                # This statement changed to be consistent with PRESTA-I
                ch_steps.count_all_steps += 1
                ch_steps.is_ch_step = False

                # if iausfl[TUSTEPB]:  # extra
                #     ausgab(TUSTEPB, tustep=tustep, tperp=tperp, skindepth=skindepth)

                if tustep <= tperp and (not exact_bca or tustep > skindepth):
                    # We are further way from a boundary than a skindepth, so
                    # perform a normal condensed-history step
                    callhowfar = False  # Do not call HAWFAR
                    domultiple = False  # Multiple scattering done here
                    dosingle   = False  # MS => no single scattering
                    callmsdist = True  # Remember that msdist has been called

                    # Fourth order technique for de
                    # --- Inline replace: $ COMPUTE_ELOSS_G(tustep,eke,elke,lelke,de); -----
                    de = compute_eloss_g(lelec, medium, tustep, eke, elke, lelke, range_)
                    # End inline replace: $ COMPUTE_ELOSS_G(tustep,eke,elke,lelke,de); ----

                    epcont.tvstep = tustep

                    if transport_algorithm == PRESTA_II:
                        # logger.debug(
                        #     f'PrestaII xyz,uvw={x[np_m1]},{y[np_m1]},{z[np_m1]},'
                        #     f'{u[np_m1]},{v[np_m1]},{w[np_m1]}'
                        # )
                        egsfortran.msdist_pii(  # msdist_pII but for case issues
                            # Inputs
                            eke,de,tustep,rhof,medium,qel,spin_effects,
                            u[np_m1],v[np_m1],w[np_m1],x[np_m1],y[np_m1],z[np_m1],
                            # Outputs
                            uscat,vscat,wscat,xtrans,ytrans,ztrans,ustep
                    )
                        # if iausfl[PRESTAIIA]:  # extra for debugging
                        #     ausgab(PRESTAIIA)
                    else:
                        egsfortran.msdist_pi(  # msdist_pI but for case issues
                            # Inputs
                            eke,de,tustep,rhof,medium,qel,spin_effects,
                            u[np_m1],v[np_m1],w[np_m1],x[np_m1],y[np_m1],z[np_m1],
                            # Outputs
                            uscat,vscat,wscat,xtrans,ytrans,ztrans,ustep
                        )
                        # if iausfl[PRESTAIA]:  # extra for debugging
                        #     ausgab(PRESTAIA)

                    # logger.debug(
                    #     f'presta out: uvwcat,xyztrans,ustep {uscat},{vscat},{wscat}'
                    #     f',{xtrans},{ytrans},{ztrans},{ustep}')
                else:
                    # We are within a skindepth from a boundary, invoke
                    # one of the various boundary-crossing algorithms
                    callmsdist = False  # Remember that msdist has not been called
                    if exact_bca:
                        # Cross the boundary in a single scattering mode
                        domultiple = False  # Do not do multiple scattering
                        # Sample the distance to a single scattering event
                        rnnoss = randomset()
                        # logger.info(f'random rnnoss={rnnoss:0.8e}')
                        if rnnoss < 1.e-30:
                            rnnoss = 1.e-30

                        lambda_ = -log(1 - rnnoss)
                        lambda_max = 0.5*blccl*rm/dedx*(eke/rm+1)**3
                        if lambda_ >= 0 and lambda_max > 0:
                            if lambda_ < lambda_max :
                                tuss=lambda_*ssmfp*(1-0.5*lambda_/lambda_max)
                            else:
                                tuss = 0.5 * lambda_ * ssmfp

                            if tuss < tustep:
                                epcont.tustep = tuss
                                dosingle = True
                            else:
                                dosingle = False

                        else:
                            logger.warning(
                                f' lambda_ > lambda_max: {lambda_},{lambda_max}'
                                f' eke dedx: {eke},{dedx}'
                                f' ir medium blcc: {ir[np_m1]},{medium},{blcc[medium_m1]}'
                                f' position = {x[np_m1]},{y[np_m1]},{z[np_m1]}'
                            )
                            dosingle = False
                            stack.np -= 1
                            return LAMBDA_WARNING, ircode, eie, peie
                        epcont.ustep = tustep
                    else:
                        # Boundary crossing a la EGS4/PRESTA-I but using
                        # exact PLC
                        dosingle = False
                        domultiple = True
                        # --- Inline replace: $ SET_USTEP; -----
                        if set_ustep:
                            set_ustep()
                        else:
                            ekems = eke - 0.5*tustep*dedx # Use mid-point energy to calculate
                                                            # energy dependent quantities
                            # --- Inline replace: $ CALCULATE_XI(tustep); -----
                            xi, blccl = calculate_xi(lelec, medium, ekems, rmt2, rmsq, xccl, blccl, tustep)
                            # End inline replace: $ CALCULATE_XI(tustep); ----
                            if xi < 0.1 :
                                epcont.ustep = tustep*(1 - xi*(0.5 - xi*0.166667))
                            else:
                                epcont.ustep = tustep*(1 - Exp(-xi))/xi
                        # End inline replace: $ SET_USTEP; ----

                    if ustep < tperp:
                        callhowfar = False
                    else:
                        callhowfar = True
            # end non-vacuum test

            # additional ustep restriction in em field
            if set_ustep_em_field:
                set_ustep_em_field()

            epcont.irold = ir[np_m1] # save current region
            epcont.irnew = ir[np_m1] # default new region is current region
            epcont.idisc  = 0 # default is no discard (this flag is initialized here)
            ustep0 = ustep # Save the intended ustep.

            # IF(callhowfar) [ call howfar; ]
            # --- Inline replace: $ CALL_HOWFAR_IN_ELECTR; -----
            if call_howfar_in_electr:
                call_howfar_in_electr()
            else:
                if callhowfar or wt[np_m1] <= 0:
                    # if iausfl[HOWFARB]:
                    #     ausgab(HOWFARB, idisc=idisc, ustep=ustep, irnew=irnew)
                    howfar()
                    # if iausfl[HOWFARA]:
                    #     ausgab(HOWFARA, idisc=idisc, ustep=ustep, irnew=irnew)

            # End inline replace: $ CALL_HOWFAR_IN_ELECTR; ----

            # Now see if user requested discard
            if idisc > 0: # idisc is returned by howfar:
                # User requested immediate discard
                return USER_ELECTRON_DISCARD, None, eie, peie

            # --- Inline replace: $ CHECK_NEGATIVE_USTEP; -----
            if check_negative_ustep:
                check_negative_ustep()
            elif ustep <= 0:
                # Negative ustep---probable truncation problem at a
                # boundary, which means we are not in the region we think
                # we are in.  The default macro assumes that user has set
                # irnew to the region we are really most likely to be
                # in.  A message is written out whenever ustep is less than -1.e-4
                if ustep < -1e-4:
                    ierust += 1
                    logger.info(
                        f"{ierust:4d} Negative ustep = {ustep:12.5e}"
                        f" dedx={dedx:8.4f} ke={e[np_m1]-prm:8.4f}"
                        f" ir,irnew,irold ={ir[np_m1]:4d},{irnew:4d},{irold:4d}"
                        f" x,y,z ={x[np_m1]:10.3e},{y[np_m1]:10.3e},{z[np_m1]:10.3e}"
                    )
                    if ierust > 1000:
                        logger.critical(
                            '\n\n\n\n Called exit---too many ustep errors\n\n\n'
                        )
                        sys.exit(1)

                epcont.ustep = 0
            # End inline replace: $ CHECK_NEGATIVE_USTEP; ----

            if ustep == 0 or medium == 0:
                # Do fast step in vacuum
                if ustep != 0:
                    if EM_MACROS_ACTIVE:
                        epcont.edep = pzero # no energy loss in vacuum
                        # transport in EMF in vacuum:
                        # only a B or and E field can be active
                        # (not both at the same time)
                    else:
                        # Step in vacuum
                        epcont.vstep = ustep
                        epcont.tvstep = vstep
                        # ( vstep is ustep truncated (possibly) by howfar
                        #  tvstep is the total curved path associated with vstep)
                        epcont.edep = pzero # no energy loss in vacuum

                        # additional vacuum transport in em field
                        epcont.e_range =  vacdst
                        # logger.info("vacuum step")
                        if call_tranausb:
                            ausgab(TRANAUSB)
                        # Transport the particle
                        x[np_m1] += u[np_m1]*vstep
                        y[np_m1] += v[np_m1]*vstep
                        z[np_m1] += w[np_m1]*vstep
                        dnear[np_m1] -= vstep
                            # (dnear is distance to the nearest boundary
                            #  that goes along with particle stack and
                            #  which the user's howfar can supply (option)

                            # default for $ SET-ANGLES-EM-FIELD; is ; (null)
                            # (allows for EM field deflection
                    # end of EM_MACROS_ACTIVE block
                # end of vacuum step

                if irnew != irold:
                    # --- Inline replace: $ electron_region_change; -----
                    if electron_region_change:
                        electron_region_change()
                    else:
                        ir[np_m1] = irnew
                        irl = irnew
                        irl_m1 = irl - 1  # ** 0-based
                        useful.medium =med[irl_m1]
                        medium_m1 = medium - 1  # ** 0-based
                    # End inline replace: $ electron_region_change; ----

                if ustep != 0:
                    if call_tranausa:
                        ausgab(TRANAUSA)
                if eie <= ecut[irl_m1]:
                    return ECUT_DISCARD, None, eie, peie
                if ustep != 0 and idisc < 0:
                    return USER_ELECTRON_DISCARD, None, eie, peie

                next_tstep = True
                break  # out of ustep (next :TSTEP:)

            # Go try another big step in (possibly) new medium
            epcont.vstep = ustep
            if callhowfar:
                if exact_bca:
                    # if callhowfar is True and exact_bca=True we are
                    # in a single scattering mode
                    epcont.tvstep = vstep
                    if tvstep != tustep:
                        # Boundary was crossed. Shut off single scattering
                        dosingle = False
                else:
                    # callhowfar=True and exact_bca=False
                    # =>we are doing an approximate CH step
                    # calculate the average curved path-length corresponding
                    # to vstep
                    # --- Inline replace: $ SET_TVSTEP; -----
                    if set_tvstep:
                        set_tvstep()
                    elif vstep < ustep0:
                        ekems = eke - 0.5*tustep*vstep/ustep0*dedx
                            # This estimates the energy loss to the boundary.
                            # tustep was the intended curved path-length,
                            # ustep0 is the average transport distance in the initial direction
                            #        resulting from tustep
                            # vstep = ustep is the reduced average transport distance in the
                            #               initial direction due to boundary crossing
                        # --- Inline replace: $ CALCULATE_XI(vstep); -----
                        xi, blccl = calculate_xi(lelec, medium, ekems, rmt2, rmsq, xccl, blccl, vstep)
                        # End inline replace: $ CALCULATE_XI(vstep); ----
                        if xi < 0.1:
                            epcont.tvstep = vstep*(1 + xi*(0.5 + xi*0.333333))
                        elif xi < 0.999999 :
                            epcont.tvstep = -vstep*log(1 - xi)/xi
                        else:
                            # This is an error condition because the average transition
                            # in the initial direction of motion is always smaller than 1/Q1
                            logger.info(
                                ' Stopped in SET-TVSTEP because xi > 1! \n'
                                f' Medium: {medium}\n'
                                f' Initial energy: {eke}\n'
                                f' Average step energy: {ekems}\n'
                                f' tustep: {tustep}\n'
                                f' ustep0: {ustep0}\n'
                                f' vstep:  {vstep}\n'
                                f' ==> xi = {xi}\n'
                            )
                            logger.critical(
                                '***************** Error: '
                                'This is a fatal error condition'
                                '***************** Quitting now.'
                            )
                            sys.exit(1)
                    else:
                        epcont.tvstep = tustep
                    # End inline replace: $ SET_TVSTEP; ----

                # Fourth order technique for dedx
                # Must be done for an approx. CH step or a
                # single scattering step.
                # --- Inline replace: $ COMPUTE_ELOSS_G(tvstep,eke,elke,lelke,de); -----
                de = compute_eloss_g(lelec, medium, tvstep, eke, elke, lelke, range_)
            else:
                # callhowfar=False => step has not been reduced due to
                #                       boundaries
                epcont.tvstep = tustep
                if not callmsdist:
                    # Second order technique for dedx
                    # Already done in a normal CH step with call to msdist
                    # --- Inline replace: $ COMPUTE_ELOSS_G(tvstep,eke,elke,lelke,de); -----
                    de = compute_eloss_g(lelec, medium, tvstep, eke, elke, lelke, range_)

            # logger.debug(f'Setting TVSTEP tvstep, de= {tvstep},{de}')
            if set_tvstep_em_field:
                set_tvstep_em_field()
                # additional path length correction in em field
                # ( Calculates tvstep given vstep
                #  default for $ SET-TVSTEP-EM-FIELD; is ; (null)

            # the energy loss is used to calculate the number
            # of MFP gone up to now. If energy loss
            # fluctuations are implemented, de will be
            # changed in $ DE-FLUCTUATION; => save
            save_de = de

            if de_fluctuation:
                de_fluctuation()
            epcont.edep = de # energy deposition variable for user
            if add_work_em_field:
                add_work_em_field  # e-loss or gain in em field
            if add_work_em_field2:
                add_work_em_field2  # EEMF implementation
            ekef = eke - de  # (final kinetic energy)
            epcont.eold = eie  # save old value
            epcont.enew = eold - de  # energy at end of transport

            # Now do multiple scattering
            if not callmsdist:
                # everything done if callmsdist  is True
                if domultiple:
                    # Approximated CH step => do multiple scattering
                    #
                    # ekems, elkems, beta2 have been set in either $ SET-TUSTEP
                    # or $ SET-TVSTEP if spin_effects is True, they are
                    # not needed if spin_effects is False
                    #
                    # chia2,etap,xi,xi_corr are also set in the above macros
                    #
                    # qel (0 for e-, 1 for e+) and medium are now also required
                    # (for the spin rejection loop)
                    #
                    lambda_ = blccl*tvstep/beta2/etap/(1+chia2)
                    xi = xi/xi_corr
                    findindex = True; spin_index = True
                    mscat(lambda_,chia2,xi,elkems,beta2,qel,medium,
                            spin_effects,findindex,spin_index,
                            costhe,sinthe)
                    # logger.debug('Multiple scattering called')
                elif dosingle:
                    # Single scattering
                    ekems = max(ekef,ecut[irl_m1]-rm)
                    p2 = ekems*(ekems + rmt2)
                    beta2 = p2/(p2 + rmsq)
                    chia2 = xcc[medium_m1]/(4*blcc[medium_m1]*p2)
                    if spin_effects:
                        elkems = log(ekems)
                        lelkems=int(eke1[medium_m1]*elkems+eke0[medium_m1])
                        lelkems_m1 = lelkems - 1  # ** 0-based
                        if lelec < 0:
                            # EVALUATE etap USING etae_ms(elkems)]
                            etap = etae_ms1[lelkems_m1,medium_m1]*elkems+ etae_ms0[lelkems_m1,medium_m1]
                        else:
                            # EVALUATE etap USING etap_ms(elkems)
                            etap = etap_ms1[lelkems_m1,medium_m1]*elkems+ etap_ms0[lelkems_m1,medium_m1]
                        chia2 = chia2*etap

                    egsfortran.sscat(chia2,elkems,beta2,qel,medium,
                                spin_effects,costhe,sinthe)
                    # logger.debug('Single scattering called')
                else:
                    uphiot.theta  = 0 # No deflection in single scattering model
                    uphiot.sinthe = 0
                    uphiot.costhe = 1

            # We now know distance and amount of energy loss for this step,
            # and the angle by which the electron will be scattered. Hence,
            # it is time to call the user and inform him of this transport,
            # after which we will do it.

            # Now transport, deduct energy loss, and do multiple scatter.
            epcont.e_range =  range_
            # ******* trying to save evaluation of range_.
            # range_ = range_ - tvstep*rhof
            # ********/


            # Put expected final position and direction in common
            # block variables so that they are available to the
            # user for things such as scoring on a grid that is
            # different from the geometry grid

            if callmsdist:
                # Deflection and scattering have been calculated/sampled in msdist
                # if iausfl[UVWAUSB]:
                #     ausgab(
                #         UVWAUSB, msg='callmsdist=True',
                #         us=uscat, vscat=vscat, wscat=wscat
                #     )
                epcont.u_final = uscat
                epcont.v_final = vscat
                epcont.w_final = wscat
                epcont.x_final = xtrans
                epcont.y_final = ytrans
                epcont.z_final = ztrans
            else:
                if not EM_MACROS_ACTIVE:
                    epcont.x_final = x[np_m1] + u[np_m1]*vstep
                    epcont.y_final = y[np_m1] + v[np_m1]*vstep
                    epcont.z_final = z[np_m1] + w[np_m1]*vstep
                    # if iausfl[XYZAUSA]:
                    #     ausgab(
                    #         XYZAUSA, msg="callmsdist=False",
                    #         xf=x_final, yf=y_final, zf=z_final
                    #     )

                if domultiple or dosingle:
                    uvw_tmp = (u[np_m1], v[np_m1], w[np_m1])
                    # logger.debug(
                    #         f"domultiple/dosingle before uvw_tmp={uvw_tmp}"
                    #     )

                    uphi(2,1) # Apply the deflection, save call to uphi if
                                    # no deflection in a single scattering mode
                    np_m1 = np - 1
                    epcont.u_final = u[np_m1]
                    epcont.v_final = v[np_m1]
                    epcont.w_final = w[np_m1]
                    # if iausfl[UVWAUSA]:
                    #     ausgab(
                    #         UVWAUSA, msg="domultiple/dosingle after", uf=u_final, vf=v_final, wf=w_final
                    #     )

                    # logger.debug(f'Called UPHI: final uvw={u_final},{v_final},{w_final}')
                    # logger.debug(f'Called UPHI: uvw={u[np_m1]},{v[np_m1]},{w[np_m1]}')
                    u[np_m1], v[np_m1], w[np_m1] = uvw_tmp
                else:
                    epcont.u_final = u[np_m1]; epcont.v_final = v[np_m1]; epcont.w_final = w[np_m1]
                    # if iausfl[UVWAUSA]:
                    #     ausgab(
                    #         UVWAUSA, msg="NOT domultiple/dosingle", uf=u_final, vf=v_final, wf=w_final
                    #     )

            if call_tranausb:
                ausgab(TRANAUSB)

            # Transport the particle
            x[np_m1] = x_final; y[np_m1] = y_final; z[np_m1] = z_final
            u[np_m1] = u_final; v[np_m1] = v_final; w[np_m1] = w_final

            # iarg=TRANAUSB
            # if iausfl[iarg-1+1] != 0:  # ** 0-based
            #     ausgab(iarg)

            # logger.debug(
            #     f"Transport: new xyz/uvw={x[np_m1]} {y[np_m1]} {z[np_m1]}"
            #     f" {u[np_m1]} {v[np_m1]} {w[np_m1]}"
            # )
            dnear[np_m1] -= vstep
            epcont.irold = ir[np_m1] # save previous region

            if set_angles_em_field:
                set_angles_em_field()

            # Now done with multiple scattering,
            # update energy and see if below cut
            # below subtracts only energy deposited
            peie -= edep
            # below subtracts energy deposited + work due to E field
            # peie = peie - de
            eie   = peie
            e[np_m1] = peie

            # if( irnew != irl and eie <= ecut[irl_m1]) [
            # IK: the above is clearly a bug. If the particle energy falls
            #     below ecut, but the particle is actually entering a new
            #     region, the discard will happen in the current region
            #     instead the next. If the particle is a positron, all
            #     resulting annihilation photons will have the new position
            #     but the old region => confusion in the geometry routine
            #     is very likely.      Jan 27 2004
            if irnew == irl and eie <= ecut[irl_m1]:
                return ECUT_DISCARD, None, eie, peie

            useful.medold = medium
            if medium != 0:
                epcont.eke = eie - rm # update kinetic energy
                epcont.elke = log(eke)
                lelke = int(eke1[medium_m1]*elke + eke0[medium_m1])
                lelke_m1 = lelke - 1

            if irnew != irold:
                # --- Inline replace: $ electron_region_change; -----
                if electron_region_change:
                    electron_region_change()
                else:
                    ir[np_m1] = irnew
                    irl = irnew
                    irl_m1 = irl - 1  # ** 0-based
                    useful.medium = med[irl_m1]
                    medium_m1 = medium - 1
                # End inline replace: $ electron_region_change; ---- ]

            # After transport call to user scoring routine
            if call_tranausa:
                ausgab(TRANAUSA)

            if eie <= ecut[irl_m1]:
                return ECUT_DISCARD, None, eie, peie

            # Now check for deferred discard request.  May have been set
            # by either howfar, or one of the transport ausgab calls
            if idisc < 0:
                return USER_ELECTRON_DISCARD, None, eie, peie

            if medium != medold:
                # logger.info('medium != medold, next tstep')
                next_tstep = True
                break  # leave ustep loop (NEXT :TSTEP:)

            if user_controls_tstep_recursion:
                user_controls_tstep_recursion()

            # --- Inline replace: $ UPDATE_DEMFP; -----
            if update_demfp:
                update_demfp()
            else:
                demfp -= save_de * sig
                total_de -= save_de
                total_tstep -= tvstep * rhof
                if total_tstep < 1e-9 :
                    demfp = 0
            # End inline replace: $ UPDATE_DEMFP; ----

            if demfp < smallest_electron_mfp:
                break  # end ustep loop


            # loop on ustep ----------------------------------------------
            # loop on ustep ----------------------------------------------

        # some `break`s out of ustep go to top of loop immediately
        if next_tstep:
            next_tstep = False
            continue  # back to start of TSTEP loop

        # Compute final sigma to see if resample is needed.
        # this will take the energy variation of the sigma into
        # account using the fictitious sigma method.

        # --- Inline replace: $ EVALUATE_SIGF; -----
        if evaluate_sigf:
            evaluate_sigf()
        else:
            if lelec < 0:
                # EVALUATE sigf USING esig(elke)
                sigf = esig1[lelke_m1, medium_m1]*elke+ esig0[lelke_m1, medium_m1]
                # EVALUATE dedx0 USING ededx(elke)
                dedx0 = ededx1[lelke_m1, medium_m1]*elke+ ededx0[lelke_m1, medium_m1]
                sigf = sigf/dedx0
            else:
                # EVALUATE sigf USING psig(elke)
                sigf = psig1[lelke_m1, medium_m1]*elke+ psig0[lelke_m1, medium_m1]
                # EVALUATE dedx0 USING pdedx(elke)
                dedx0 = pdedx1[lelke_m1, medium_m1]*elke+ pdedx0[lelke_m1, medium_m1]
                sigf = sigf/dedx0

        # End inline replace: $ EVALUATE_SIGF; ----

        sigratio = sigf / sig0
        rfict = randomset()
        # logger.info(f'random rfict={rfict:0.8e}')

        if rfict <= sigratio:
            return None, lelke, eie, peie   # end tstep loop, else continue looping

