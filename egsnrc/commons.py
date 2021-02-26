from egsnrc import egsfortran

#  COMMON block randomm --------
randomm = egsfortran.randomm
rng_seed = randomm.rng_seed
rng_array = randomm.rng_array
seeds = randomm.seeds

#  COMMON block stack --------
stack = egsfortran.stack
e = stack.e
x = stack.x
y = stack.y
z = stack.z
u = stack.u
v = stack.v
w = stack.w
dnear = stack.dnear
wt = stack.wt
iq = stack.iq
ir = stack.ir
latch = stack.latch
latchi = stack.latchi
np = stack.np
npold = stack.npold

#  COMMON block media --------
media_ = egsfortran.media  # underscore because `media`* is also an attribute
rlc = media_.rlc
rldu = media_.rldu
rho = media_.rho
msge = media_.msge
mge = media_.mge
mseke = media_.mseke
meke = media_.meke
mleke = media_.mleke
mcmfp = media_.mcmfp
mrange = media_.mrange
iraylm = media_.iraylm
iphotonucm = media_.iphotonucm
media = media_.media  # * see note for `media_`
iphotonuc = media_.iphotonuc
nmed = media_.nmed
eii_xfile = media_.eii_xfile
photon_xsections = media_.photon_xsections
comp_xsections = media_.comp_xsections
photonuc_xsections = media_.photonuc_xsections

#  COMMON block bounds --------
bounds = egsfortran.bounds
ecut = bounds.ecut
pcut = bounds.pcut
vacdst = bounds.vacdst

#  COMMON block elecin --------
elecin = egsfortran.elecin
esige_max = elecin.esige_max
psige_max = elecin.psige_max
esig_e = elecin.esig_e
psig_e = elecin.psig_e
range_ep = elecin.range_ep
e_array = elecin.e_array
etae_ms0 = elecin.etae_ms0
etae_ms1 = elecin.etae_ms1
etap_ms0 = elecin.etap_ms0
etap_ms1 = elecin.etap_ms1
q1ce_ms0 = elecin.q1ce_ms0
q1ce_ms1 = elecin.q1ce_ms1
q1cp_ms0 = elecin.q1cp_ms0
q1cp_ms1 = elecin.q1cp_ms1
q2ce_ms0 = elecin.q2ce_ms0
q2ce_ms1 = elecin.q2ce_ms1
q2cp_ms0 = elecin.q2cp_ms0
q2cp_ms1 = elecin.q2cp_ms1
blcce0 = elecin.blcce0
blcce1 = elecin.blcce1
eke0 = elecin.eke0
eke1 = elecin.eke1
xr0 = elecin.xr0
teff0 = elecin.teff0
blcc = elecin.blcc
xcc = elecin.xcc
esig0 = elecin.esig0
esig1 = elecin.esig1
psig0 = elecin.psig0
psig1 = elecin.psig1
ededx0 = elecin.ededx0
ededx1 = elecin.ededx1
pdedx0 = elecin.pdedx0
pdedx1 = elecin.pdedx1
ebr10 = elecin.ebr10
ebr11 = elecin.ebr11
pbr10 = elecin.pbr10
pbr11 = elecin.pbr11
pbr20 = elecin.pbr20
pbr21 = elecin.pbr21
tmxs0 = elecin.tmxs0
tmxs1 = elecin.tmxs1
expeke1 = elecin.expeke1
iunrst = elecin.iunrst
epstfl = elecin.epstfl
iaprim = elecin.iaprim
sig_ismonotone = elecin.sig_ismonotone

#  COMMON block thresh --------
thresh = egsfortran.thresh
rmt2 = thresh.rmt2
rmsq = thresh.rmsq
ap = thresh.ap
ae = thresh.ae
up = thresh.up
ue = thresh.ue
te = thresh.te
thmoll = thresh.thmoll

#  COMMON block uphiot --------
uphiot = egsfortran.uphiot
theta = uphiot.theta
sinthe = uphiot.sinthe
costhe = uphiot.costhe
sinphi = uphiot.sinphi
cosphi = uphiot.cosphi
pi = uphiot.pi
twopi = uphiot.twopi
pi5d2 = uphiot.pi5d2

#  COMMON block useful --------
#  from tutor1.f: DATA RM,PRM,PRMT2,PZERO/0.5109989461,0.5109989461,1.0219978922,0.D0/
useful = egsfortran.useful
pzero = useful.pzero = 0.0
prm = useful.prm = 0.51099896430969238
prmt2 = useful.prmt2 = 1.0219979286193848  # from fortran gdb, was 1.0219978922
rm = useful.rm = 0.51099896430969238  # was 0.5109989461, updated from fortr gdb
medium = useful.medium
medold = useful.medold

#  COMMON block misc --------
misc = egsfortran.misc
dunit = misc.dunit
kmpi = misc.kmpi
kmpo = misc.kmpo
rhor = misc.rhor
med = misc.med
iraylr = misc.iraylr
iphotonucr = misc.iphotonucr

#  COMMON block epcont --------
epcont = egsfortran.epcont
edep = epcont.edep
edep_local = epcont.edep_local
tstep = epcont.tstep
tustep = epcont.tustep
ustep = epcont.ustep
tvstep = epcont.tvstep
vstep = epcont.vstep
rhof = epcont.rhof
eold = epcont.eold
enew = epcont.enew
eke = epcont.eke
elke = epcont.elke
gle = epcont.gle
e_range = epcont.e_range
x_final = epcont.x_final
y_final = epcont.y_final
z_final = epcont.z_final
u_final = epcont.u_final
v_final = epcont.v_final
w_final = epcont.w_final
idisc = epcont.idisc
irold = epcont.irold
irnew = epcont.irnew
iausfl = epcont.iausfl

#  COMMON block et_control --------
et_control = egsfortran.et_control
estepe = et_control.estepe
ximax = et_control.ximax
skindepth_for_bca = et_control.skindepth_for_bca
transport_algorithm = et_control.transport_algorithm
bca_algorithm = et_control.bca_algorithm
exact_bca = et_control.exact_bca
spin_effects = et_control.spin_effects
smaxir = et_control.smaxir


#  COMMON block egs_vr --------
egs_vr = egsfortran.egs_vr
prob_rr = egs_vr.prob_rr
nbr_split = egs_vr.nbr_split
i_play_rr = egs_vr.i_play_rr
i_survived_rr = egs_vr.i_survived_rr
n_rr_warning = egs_vr.n_rr_warning
e_max_rr = egs_vr.e_max_rr
i_do_rr = egs_vr.i_do_rr


#  COMMON block ch_steps --------
ch_steps = egsfortran.ch_steps
count_pii_steps = ch_steps.count_pii_steps
count_all_steps = ch_steps.count_all_steps
is_ch_step = ch_steps.is_ch_step


#  COMMON block nist_brems --------
nist_brems = egsfortran.nist_brems
nb_fdata = nist_brems.nb_fdata
nb_xdata = nist_brems.nb_xdata
nb_wdata = nist_brems.nb_wdata
nb_idata = nist_brems.nb_idata
nb_emin = nist_brems.nb_emin
nb_emax = nist_brems.nb_emax
nb_lemin = nist_brems.nb_lemin
nb_lemax = nist_brems.nb_lemax
nb_dle = nist_brems.nb_dle
nb_dlei = nist_brems.nb_dlei
log_ap = nist_brems.log_ap


#  COMMON block brempr --------
brempr = egsfortran.brempr
ibrdst = brempr.ibrdst
iprdst = brempr.iprdst
ibr_nist = brempr.ibr_nist
pair_nrc = brempr.pair_nrc
itriplet = brempr.itriplet
dl1 = brempr.dl1
dl2 = brempr.dl2
dl3 = brempr.dl3
dl4 = brempr.dl4
dl5 = brempr.dl5
dl6 = brempr.dl6
alphi = brempr.alphi
bpar = brempr.bpar
delpos = brempr.delpos
wa = brempr.wa
pz = brempr.pz
zelem = brempr.zelem
rhoz = brempr.rhoz
pwr2i = brempr.pwr2i
delcm = brempr.delcm
zbrang = brempr.zbrang
lzbrang = brempr.lzbrang
nne = brempr.nne
asym = brempr.asym


#  COMMON block eii_data --------
eii_data = egsfortran.eii_data
eii_l_factor = eii_data.eii_l_factor
eii_flag = eii_data.eii_flag
eii_xsection_a = eii_data.eii_xsection_a
eii_xsection_b = eii_data.eii_xsection_b
eii_cons = eii_data.eii_cons
eii_a = eii_data.eii_a
eii_b = eii_data.eii_b
eii_z = eii_data.eii_z
eii_sh = eii_data.eii_sh
eii_nshells = eii_data.eii_nshells
eii_nsh = eii_data.eii_nsh
eii_first = eii_data.eii_first
eii_no = eii_data.eii_no


#  COMMON block edge --------
edge = egsfortran.edge
binding_energies = edge.binding_energies
interaction_prob = edge.interaction_prob
relaxation_prob = edge.relaxation_prob
edge_energies = edge.edge_energies
edge_number = edge.edge_number
edge_a = edge.edge_a
edge_b = edge.edge_b
edge_c = edge.edge_c
edge_d = edge.edge_d
iedgfl = edge.iedgfl
iphter = edge.iphter


#  COMMON block photin --------
photin = egsfortran.photin
dpmfp = photin.dpmfp
ebinda = photin.ebinda
ge0 = photin.ge0
ge1 = photin.ge1
gmfp0 = photin.gmfp0
gbr10 = photin.gbr10
gbr20 = photin.gbr20
gmfp1 = photin.gmfp1
gbr11 = photin.gbr11
gbr21 = photin.gbr21
rco0 = photin.rco0
rco1 = photin.rco1
rsct0 = photin.rsct0
rsct1 = photin.rsct1
cohe0 = photin.cohe0
cohe1 = photin.cohe1
photonuc0 = photin.photonuc0
photonuc1 = photin.photonuc1
mpgem = photin.mpgem
ngr = photin.ngr