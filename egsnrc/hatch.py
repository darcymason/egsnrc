from pathlib import Path
import numpy as np

DATA_DIR = Path(__file__).resolve().parent / "data"



def get_xsection_table(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    i_line = 0
    data = []
    for z in range(1, 101):
        count = int(lines[i_line])
        i_line += 1
        # 2 values per item, 8 values stored per line, so 4 data point pairs
        # Calc number of lines needed
        data_lines = count // 4 + (1 if count % 4 else 0)
        z_data = np.loadtxt(
            x
            for i in range(data_lines)
            for x in lines[i_line + i].strip().split()
        )
        # for i in range(data_lines):
        #     z_data.append(float(x) for x in lines[i_line+i].strip().split())
        # z_data = np.loadtxt(lines, skiprows=i_line, max_rows=data_lines)
        z_data = z_data.reshape((-1, 2)).transpose()
        data.append(z_data)
        print(f"Count {count}, len(data): {len(z_data)}")
        i_line += data_lines

    return data

# egsi_get_data copied here for reference on how post-processing is done

# subroutine egsi_get_data(flag,iunit,n,ne,zsorted,pz_sorted,ge1,ge0,data);
# "=========================================================================="
# implicit none;
# COMIN/EGS-IO/;
# $REAL    eth;
# $INTEGER flag,iunit,n,ne;
# $REAL    ge1,ge0,zsorted(*),pz_sorted(*),data(*);
# $REAL    etmp($MXINPUT),ftmp($MXINPUT);
# $REAL    gle,sig,p,e;
# $INTEGER i,j,k,kk,iz,iz_old,ndat,iiz;

# ;COMIN/USEFUL/;

# "Ali:photonuc. The whole routine is commented out and re-written
# "to accommodate reading photonuclear cross sections. A copy of the
# "commented original routine is at the bottom.
# " flag = 0: photoelectric, Rayleigh, Compton
# " flag = 1: pair
# " flag = 2: triplet
# " flag = 3: photonuclear

# rewind(iunit);
# iz_old = 0;
# DO k=1,n [ data(k) = 0; ]
# DO i=1,ne [
#     iiz = int(zsorted(i)+0.5);
#     DO iz=iz_old+1,iiz [
#         read(iunit,*,err=:user-data-failure:) ndat;
#         IF( ndat > $MXINPUT ) [
#             $egs_fatal(*,'Too many input data points. Max. is ',$MXINPUT);
#         ]
#         IF( flag = 0 | flag = 3) [
#             read(iunit,*,err=:user-data-failure:) (etmp(k),ftmp(k),k=1,ndat);
#         ]
#         ELSE [
#             read(iunit,*,err=:user-data-failure:) (etmp(k+1),ftmp(k+1),
#                 k=1,ndat);
#             IF( flag = 1 ) [ eth = 2*rm; ] ELSE [ eth = 4*rm; ]
#             ndat = ndat + 1;
#             DO k=2,ndat [
#                 ftmp(k) = ftmp(k) - 3*log(1-eth/exp(etmp(k)));
#             ]
#             ftmp(1) = ftmp(2); etmp(1) = log(eth);
#         ]
#     ]
#     iz_old = iiz;
#     DO k=1,n [
#         gle = (k - ge0)/ge1; e = exp(gle);
#         IF( gle < etmp(1) | gle >= etmp(ndat) ) [
#             IF( flag = 0 ) [
#                 $egs_fatal(*,'Energy ',exp(gle),
#                    ' is outside the available data range of ',
#                    exp(etmp(1)),exp(etmp(ndat)));
#             ]
#             ELSEIF (flag = 1 | flag = 2) [
#                 IF( gle < etmp(1) ) [ sig = 0; ]
#                 ELSE [ sig = exp(ftmp(ndat)); ]
#             ]
#             ELSE [ "photonuclear, zero it before and after
#              sig = 0;
#             ]
#         ] ELSE [
#             DO kk=1,ndat-1 [
#                 IF( gle >= etmp(kk) & gle < etmp(kk+1) ) EXIT;
#             ]
#             IF( flag ~= 3) ["log/log interpolation"
#                p = (gle - etmp(kk))/(etmp(kk+1) - etmp(kk));
#                sig = exp(p*ftmp(kk+1) + (1-p)*ftmp(kk));
#             ]
#             ELSE["lin/lin interpolation for photonuc"
#                p = (e - exp(etmp(kk)))/(exp(etmp(kk+1)) - exp(etmp(kk)));
#                sig = p*exp(ftmp(kk+1)) + (1-p)*exp(ftmp(kk));
#             ]
#         ]
#         IF( (flag = 1 | flag = 2) & e > eth ) sig = sig*(1-eth/e)**3;
#         data(k) = data(k) + pz_sorted(i)*sig;
#     ]
# ]

# return;

# :user-data-failure:;
# $egs_fatal(*,'Error while reading user photon cross sections from unit ',
#      iunit);

# return; end egsi_get_data;



def egs_init_user_photon(prefix,comp_prefix,photonuc_prefix,out):
    """
    $INTEGER      out;
    ;COMIN/BREMPR,EDGE,EGS-IO,MEDIA,PHOTIN,THRESH,COMPTON-DATA,X-OPTIONS/;
    $INTEGER   lnblnk1,egs_get_unit,medium,
            photo_unit,pair_unit,rayleigh_unit,triplet_unit,
            ounit,egs_open_file,compton_unit,
    "Ali:photonuc, 1 line"
            photonuc_unit;
    $INTEGER   nge,sorted($MXEL),i,j,k,iz,iz_old,ndat;
    $REAL      z_sorted($MXEL),pz_sorted($MXEL);
    $REAL      sig_photo($MXGE),sig_pair($MXGE),sig_triplet($MXGE),
            sig_rayleigh($MXGE),sig_compton($MXGE);
    $REAL      sigma,cohe,gmfp,gbr1,gbr2,sig_KN,gle,e,sig_p;
    $REAL      cohe_old,gmfp_old,gbr1_old,gbr2_old,
    "Ali:photonuc, 3 lines"
            sig_photonuc($MXGE),
            photonuc,
            photonuc_old;

    $REAL      etmp($MXINPUT),ftmp($MXINPUT);
    $REAL      sumZ,sumA,con1,con2,egs_KN_sigma0;
    $REAL      bc_emin,bc_emax,bc_dle,bc_data($MXBCINP),bc_tmp($MXBCINP),bcf,aj;
    $INTEGER   bc_ne;
    $LOGICAL   input_compton_data,
    "Ali:photonuc, 1 line"
            input_photonuc_data;
    character  data_dir*128,photo_file*140,pair_file*140,rayleigh_file*144,
            triplet_file*142,tmp_string*144,compton_file*144,
    "Ali:photonuc, 1 line"
            photonuc_file*144;
    """

    print('(Re)-initializing photon cross sections')
    print(' with files from the series: ', prefix)

    print(' Compton cross sections: ', comp_prefix)

    # "Ali:photonuc, 1 block"
    # IF(iphotonuc = 1) [
    # $egs_info('(a,a)',' Photonuclear cross sections: ',
    # $cstring(photonuc_prefix));
    # input_photonuc_data = .false.;
    # IF(lnblnk1(photonuc_prefix) > 0 & photonuc_prefix(1:7) ~= 'default') [
    # input_photonuc_data = .true.;
    # ]
    # ]

    input_compton_data = False
    # IF( ibcmp(1) > 1 & lnblnk1(comp_prefix) > 0 ) [
    #     IF( comp_prefix(1:7) ~= 'default' ) input_compton_data = .true.;
    # ]
    photo_file = DATA_DIR / f"{prefix}_photo.data"
    pair_file = DATA_DIR / f"{prefix}_pair.data"
    triplet_file = DATA_DIR / f"{prefix}_triplet.data"
    rayleigh_file = DATA_DIR / f"{prefix}_rayleigh.data"
    if input_compton_data:
        compton_file = DATA_DIR / f"{prefix}_compton.data"
    else:
        compton_file = DATA_DIR / 'compton_sigma.data'

    print(f" Using Compton cross sections from {compton_file}")

    # "Ali:photonuc, 1 block"
    # IF(iphotonuc = 1) [
    # IF( input_photonuc_data ) [
    #     photonuc_file = $cstring(data_dir) // $cstring(photonuc_prefix) //
    #                     '_photonuc.data';
    # ]
    # ELSE [
    #     photonuc_file = $cstring(data_dir) // 'iaea_photonuc.data';
    # ]
    # $egs_info('(a,a)',' Using photonuclear cross sections from ',
    # $cstring(photonuc_file));
    # ]

    # IF( ibcmp(1) > 1 ) [ $OPEN-UNIT(compton_unit,88,compton_file); ]
    # " Note: ibcmp > 1 means the user wants to use Bound Compton scattering "
    # "       without rejections. For this we have to use the actual bound   "
    # "       Compton scattering cross section, which is now available in a  "
    # "       file called bound_compton.data (the file actually contains the "
    # "       ratio of the Bound Compton to the KN cross section).           "
    # "       Because this option is not available on a region by region     "
    # "       basis, we just need to check ibcmp(1)                          "
    # "Ali:photonuc, 1 line"
    # IF( iphotonuc = 1 ) [ $OPEN-UNIT(photonuc_unit,89,photonuc_file); ]

    # IF( out = 1 ) [
    #     ounit = egs_open_file(87,0,1,'.xsections');
    #     write(ounit,'(/a,a,a)') 'Photon cross sections initialized from ',
    #     $cstring(prefix),' data files';
    #     write(ounit,'(a,/)')
    # '============================================================================';
    #     write(ounit,'(a,/)') 'Grid energies and cross sections are output';
    # "Ali:photonuc, 1 block"
    #     IF(iphotonuc = 1) [
    #     write(ounit,'(5x,a,t19,a,t34,a,t49,a,t64,a,t79,a)')
    #         'Energy',' GMFP(cm) ',' Pair ','Compton',' GMFP(cm) ',
    #         ' GMFP(cm) ';
    #     write(ounit,'(5x,a,t19,a,t34,a,t49,a,t64,a,t79,a/)')
    #         '(MeV)','no Rayleigh','(fraction)','(fraction)','with Rayleigh',
    #         'w/ Ray + photnuc';
    #     ]
    #     ELSE[
    #     write(ounit,'(5x,a,t19,a,t34,a,t49,a,t64,a)')
    #             'Energy',' GMFP(cm) ',' Pair ','Compton',' GMFP(cm) ';
    #     write(ounit,'(5x,a,t19,a,t34,a,t49,a,t64,a/)')
    #             '(MeV)','no Rayleigh','(fraction)','(fraction)','with Rayleigh';
    # ]
    # ]
    # /* Replace binding energies with the edges in the photo-absorption file */
    # DO iz=1,100 [
    #     read(photo_unit,*) ndat;
    #     read(photo_unit,*) (etmp(k),ftmp(k),k=1,ndat);
    #     k = 0;
    #     DO j=ndat,2,-1 [
    #         IF( etmp(j)-etmp(j-1) < 1e-5 ) [
    #             k = k+1;
    #             IF ( k <= $MXSHXSEC )[
    #             binding_energies(k,iz) = exp(etmp(j));
    #             ]
    #             ELSE[
    #             $egs_fatal('(i3,a,i3,//a)',
    #                         k,' binding energies read exceeding array size of',
    #                         $MXSHXSEC,'Increase $MXSHXSEC in egsnrc.macros!');
    #             ]
    #             IF( ~eadl_relax & k >= 4 ) EXIT;
    #         ]
    #     ]
    # ]

    # IF (mcdf_pe_xsections)[call egs_read_shellwise_pe();]

    # DO medium = 1,nmed [

    #     mge(medium) = $MXGE; nge = $MXGE;
    #     ge1(medium) = nge-1; ge1(medium) = ge1(medium)/log(up(medium)/ap(medium));
    #     ge0(medium) = 1 - ge1(medium)*log(ap(medium));

    #     $egs_info('(a,i3,a,$)',' Working on medium ',medium,' ... ');
    #     IF( out = 1 ) [
    #         write(ounit,'(/,2x,a,i3,a,24a1/)') 'Medium ',medium,': ',
    #         (media(k,medium),k=1,24);
    #     ]
    #     /sumZ,sumA/ = 0;
    #     DO i=1,nne(medium) [
    #         z_sorted(i) = zelem(medium,i);
    #         sumZ = sumZ + pz(medium,i)*zelem(medium,i);
    #         sumA = sumA + pz(medium,i)*wa(medium,i);
    #     ]
    #     con1 = sumZ*rho(medium)/(sumA*1.6605655);
    #     con2 = rho(medium)/(sumA*1.6605655);
    #     call egs_heap_sort(nne(medium),z_sorted,sorted);
    #     DO i=1,nne(medium) [ pz_sorted(i) = pz(medium,sorted(i)); ]

    #     IF (mcdf_pe_xsections)[
    #     call egsi_get_shell_data(medium,nge,nne(medium),z_sorted,pz_sorted,
    #                                 ge1(medium),ge0(medium),sig_photo);
    #     ]
    #     ELSE[
    #     call egsi_get_data(0,photo_unit,nge,nne(medium),z_sorted,pz_sorted,
    #                     ge1(medium),ge0(medium),sig_photo);
    #     ]
    #     call egsi_get_data(0,rayleigh_unit,nge,nne(medium),z_sorted,pz_sorted,
    #                     ge1(medium),ge0(medium),sig_rayleigh);
    #     call egsi_get_data(1,pair_unit,nge,nne(medium),z_sorted,pz_sorted,
    #                     ge1(medium),ge0(medium),sig_pair);
    #     call egsi_get_data(2,triplet_unit,nge,nne(medium),z_sorted,pz_sorted,
    #                     ge1(medium),ge0(medium),sig_triplet);
    # "Ali:photonuc, 1 block"
    #     IF( iphotonuc = 1 ) [
    #     call egsi_get_data(3,photonuc_unit,nge,nne(medium),z_sorted,pz_sorted,
    #                         ge1(medium),ge0(medium),sig_photonuc);
    #     ]

    #     IF( ibcmp(1) > 1 ) [
    #         "Get the bound compton cross section data"
    #         IF( input_compton_data ) [
    #             call egsi_get_data(0,compton_unit,nge,nne(medium),
    #                     z_sorted,pz_sorted,ge1(medium),ge0(medium),
    #                     sig_compton);
    #         ]
    #         ELSE [
    #             rewind(compton_unit);
    #             read(compton_unit,*) bc_emin,bc_emax,bc_ne;
    #             IF( bc_ne > $MXBCINP ) [
    #             $egs_fatal(*,'Number of input Compton data exceeds array size');
    #             ]
    #             "write(6,*) 'bc emin,emax,ne = ',bc_emin,bc_emax,bc_ne;
    #             bc_dle = log(bc_emax/bc_emin)/(bc_ne-1);
    #             DO j=1,bc_ne [ bc_data(j) = 0; ]
    #             iz_old = 1;
    #             DO i=1,nne(medium) [
    #                 iz = int(z_sorted(i)+0.5);
    #                 "write(6,*) ' reading bc data for ',iz;
    #                 DO j=iz_old,iz [ read(compton_unit,*) (bc_tmp(k),k=1,bc_ne); ]
    #                 DO j=1,bc_ne [
    #                     bc_data(j)=bc_data(j)+pz_sorted(i)*z_sorted(i)*bc_tmp(j);
    #                 ]
    #                 iz_old = iz+1;
    #             ]
    #             DO j=1,bc_ne [ bc_data(j)=log(bc_data(j)/sumZ); ]
    #         ]
    #     ]

    #     /* prepare data needed for Rayleigh scattering sampling */
    #     call egs_init_rayleigh(medium,sig_rayleigh);

    #     DO i=1,nge [

    #         gle = (i - ge0(medium))/ge1(medium); e = exp(gle);
    #         sig_KN = sumZ*egs_KN_sigma0(e);
    #         IF( ibcmp(1) > 1 ) [
    #             IF( input_compton_data ) [
    #                 sig_KN = sig_compton(i);
    #             ]
    #             ELSE [
    #                 "Apply the bound Compton correction to sig_KN"
    #                 IF( e <= bc_emin ) [ bcf = exp(bc_data(1)); ]
    #                 ELSE IF( e < bc_emax ) [
    #                     aj = 1 + log(e/bc_emin)/bc_dle;
    #                     j = int(aj); aj = aj - j;
    #                     bcf = exp(bc_data(j)*(1-aj) + bc_data(j+1)*aj);
    #                 ]
    #                 ELSE [ bcf = 1; ]
    #                 sig_KN = sig_KN*bcf;
    #                 "write(6,*) 'e = ',e,' bcf = ',bcf;
    #             ]
    #         ]
    #         sig_p  = sig_pair(i) + sig_triplet(i);
    #         sigma  = sig_KN + sig_p + sig_photo(i);
    #         gmfp   = 1/(sigma*con2);
    #         gbr1   = sig_p/sigma;
    #         gbr2   = gbr1 + sig_KN/sigma;
    #         cohe   = sigma/(sig_rayleigh(i) + sigma);
    # "Ali:photonuc, 1 line"
    #         photonuc = sigma/(sig_photonuc(i) + sigma);

    #         IF( out = 1 ) [
    # "Ali:photonuc, 1 block"
    #         IF(iphotonucm(medium) = 1) [
    #             write(ounit,'(6(1pe15.6))') e,gmfp,gbr1,gbr2-gbr1,
    #             gmfp*cohe,gmfp*cohe*photonuc;
    #         ]
    #         ELSE[
    #             write(ounit,'(5(1pe15.6))') e,gmfp,gbr1,gbr2-gbr1,gmfp*cohe;
    #         ]
    #         ]
    #         IF( i > 1 ) [
    #             gmfp1(i-1,medium) = (gmfp - gmfp_old)*ge1(medium);
    #             gmfp0(i-1,medium) =  gmfp - gmfp1(i-1,medium)*gle;
    #             gbr11(i-1,medium) = (gbr1 - gbr1_old)*ge1(medium);
    #             gbr10(i-1,medium) =  gbr1 - gbr11(i-1,medium)*gle;
    #             gbr21(i-1,medium) = (gbr2 - gbr2_old)*ge1(medium);
    #             gbr20(i-1,medium) =  gbr2 - gbr21(i-1,medium)*gle;
    #             cohe1(i-1,medium) = (cohe - cohe_old)*ge1(medium);
    #             cohe0(i-1,medium) =  cohe - cohe1(i-1,medium)*gle;
    # "Ali:photonuc, 2 lines"
    #             photonuc1(i-1,medium) = (photonuc - photonuc_old)*ge1(medium);
    #             photonuc0(i-1,medium) =  photonuc - photonuc1(i-1,medium)*gle;
    #         ]
    #         gmfp_old = gmfp; gbr1_old = gbr1; gbr2_old = gbr2; cohe_old = cohe;
    # "Ali:photonuc, 1 line"
    #         photonuc_old = photonuc;
    #     ]

    #     gmfp1(nge,medium) = gmfp1(nge-1,medium);
    #     gmfp0(nge,medium) = gmfp - gmfp1(nge,medium)*gle;
    #     gbr11(nge,medium) = gbr11(nge-1,medium);
    #     gbr10(nge,medium) = gbr1 - gbr11(nge,medium)*gle;
    #     gbr21(nge,medium) = gbr21(nge-1,medium);
    #     gbr20(nge,medium) = gbr2 - gbr21(nge,medium)*gle;
    #     cohe1(nge,medium) = cohe1(nge-1,medium);
    #     cohe0(nge,medium) = cohe - cohe1(nge,medium)*gle;
    # "Ali:photonuc, 2 lines"
    #     photonuc1(nge,medium) = photonuc1(nge-1,medium);
    #     photonuc0(nge,medium) = photonuc - photonuc1(nge,medium)*gle;

    #     $egs_info('(a)','OK');
    # ]

    # close(photo_unit); close(pair_unit);
    # close(trip
    #
    #
    # let_unit); close(rayleigh_unit);
    # "Ali:photonuc, 1 line"
    # IF( iphotonuc = 1 ) [ close(photonuc_unit);]
    # IF( ibcmp(1) > 1 ) [ close(compton_unit);]
    # IF( out = 1 )      [ close(ounit); ]
    # return;

    # :no-user-data-file:;
    # $egs_fatal('(//a,a)','Failed to open data file ',$cstring(tmp_string));

    # # return; end;

if __name__ == "__main__":
    # data_filename = DATA_DIR / "xcom_photo.data"
    data_filename = DATA_DIR / "xcom_pair.data"
    data = get_xsection_table(data_filename)
    for i in range(10):
        print(data[i][0][:8])
