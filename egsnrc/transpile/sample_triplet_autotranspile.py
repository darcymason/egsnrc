
# EMPTY CALLBACKS ----
check_stack = None


# CALLBACKS ---- 

def randomset():
    global rng_seed

    if rng_seed > 24:
        ranlux(rng_array)
        rng_seed = 1

    random_num = rng_array[rng_seed-1]
    rng_seed += 1

    return random_num



# ***************************************************************************
#                                                                            
#  Sampling of triplet production events.                                    
#                                                                            
#  The treatment is based on Borsellino's first Born approximation           
#  result (see Eq. 4B-3002 in the pair article of Motz, Olsen and Koch)        
#  As the kinematic of the process is already complicated enough and the     
#  cross section itself is not simple either, a Markov-chain method is used  
#  to sample triplet events from the Borsellino equation without any         
#  additional approximations (other then the use of the first Born           
#  approximation and the assumption of free electrons implied by             
#  Borsellino's derivation)                                                  
#                                                                            
#  Iwan Kawrakow, April 2005.                                                
# ***************************************************************************

subroutine sample_triplet

# ***************************************************************************
implicit none

# ;comin/#,#/


#  We use double precision throughout as in many cases the kinematically 
#  permitted angular interval is too small to be resolved accurately enough 
#  in single precision 

real*8 fmax_array(MAX_TRIPLET), eta_p_array(MAX_TRIPLET),
       eta_Ep_array(MAX_TRIPLET), eta_costp_array(MAX_TRIPLET),
       eta_costm_array(MAX_TRIPLET), ebin_array(MAX_TRIPLET),
       wp_array(MAX_TRIPLET), qmin_array(MAX_TRIPLET)

real*8 kmin, kmax, dlogki, alogkm, prmi, tiny_eta

real*8 ai,rnno,k,qmin,qmax,aux,a1,a2,a3,D,px1,px2,pp_min,pp_max,
       Ep_min,Ep_max,k2p2,k2p2x,peig,b,aux1,aux12,D1,aux3,xmin,xmax,
       aux6,aux7,uu,cphi,sphi,cphi_factor,aux5,phi,tmp
real*8 Er,pr,pr2,eta_pr
real*8 Ep,pp,pp2,wEp,cost_p,sint_p,eta_Ep,mup_min,wmup,
       eta_costp,Epp,pp_sintp,pp_sntp2
real*8 Em,pm,pm2,cost_m,sint_m,Emm,wmum,pm_sintm,
       eta_costm
real*8 k2,k3,s2,s3,k2k3i,k22,k32,q2,aux4,S_1,S_2,sigma
real*8 ppx, ppy, ppz, pmx, pmy, pmz, prx, pry, prz,
       a,c,sindel,cosdel,sinpsi

;integer*4 i
;logical use_it
;integer*4 iscore #  needed for BEAM 

;logical is_initialized
data is_initialized/False/
save is_initialized,fmax_array,eta_p_array,eta_Ep_array,eta_costp_array,
     eta_costm_array,ebin_array,wp_array,qmin_array,
     kmin,kmax,dlogki,alogkm,prmi,tiny_eta

if  ~is_initialized :

    is_initialized = True
    tiny_eta = 1e-6
    #  Set current cross section value to -1 in each energy bin 
    DO i=1,MAX_TRIPLET [ fmax_array(i) = -1; ]
    #  Find the maximum energy of the cross section data 
    kmax = 0; kmin = 4.1*prm
    DO i=1,nmed [ if( up(i) > kmax ) kmax == UP(i); ]
    if  kmax <= kmin :
         return
    dlogki = MAX_TRIPLET - 1; dlogki = dlogki/log(kmax/kmin)
    alogkm = 1 - dlogki*log(kmin)
    prmi = 1/prm
    DO i=1,MAX_TRIPLET [
        k = 4.1*exp((i-1.)/dlogki); ebin_array(i) = k
        qmin = 4*k/(k*(k-1)+(k+1)*sqrt(k*(k-4)))
        qmax = (k*(k-1) + (k+1)*sqrt(k*(k-4)))/(2*k+1)
        qmin_array(i) = qmin; wp_array(i) = log(qmax/qmin)


peig = e[np]
if  peig <= 4*prm :
     return
# --- Inline replace: $ CHECK_STACK(np+2,'sample_triplet'); -----
if check_stack:
    <XXX> = check_stack(np+2, 'sample_triplet')
else:
    
      if  np+2 > MXSTACK :

          $egs_fatal('(//,3a,/,2(a,i9))',' In subroutine ','sample_triplet',
              ' stack size exceeded! ',' $MAXSTACK = ',MXSTACK,' np = ',np+2)

# End inline replace: $ CHECK_STACK(np+2,'sample_triplet'); ----

#  Determine energy bin 
if  peig <= kmin :
     i = 1; 
else:
     i = MAX_TRIPLET; 
else:
    ai = alogkm + dlogki*gle; i = ai; ai = ai - i
    rnno = randomset()
    if  rnno < ai :
         i = i+1; 

#  First use the bin energy to sample the random numbers 
#  that determine recoil momentum and electron/postron angles 
k = ebin_array(i)

/*
   In the following:  k is incident photon energy in units of m*c^2
                      (all energies are in units of m*c^2, momenta in
                       units of m*c)
                      Er,pr is energy, momentum of the recoil electron
                      Ep,pp is energy, momentum of the pair positron
                      Em,pm is energy, momentum of the pair electron
                      cost_p, sint_p is cos, sin of the positron angle
                                     with respect to k
                      cost_m, sint_m same but for the electron
                      cphi is cos of azimuthal angle between positron
                              and pair electron directions.
*/

if goto_retry_triplet:  XXX

#  Pick the recoil electron momentum from 1/p.
eta_pr = randomset() if( eta_pr < tiny_eta ) eta_pr == tiny_eta
pr = qmin_array(i)*exp(eta_pr*wp_array(i))
pr2 = pr*pr; Er = sqrt(1+pr2)

#  Determine min./max. kinematically permitted postron energy for 
#  this k and p 
aux = Er-pr-1; a1=(k-pr)*(1-Er-k*aux); a2=1+k-Er; a3=1/(aux*(pr+Er-2*k-1))
D = a2*sqrt(aux*(2*k*Er+k*k*aux-pr*(Er+pr+1)/2))
px1 = (a1 + D)*a3; px2 = (a1 - D)*a3
if  px1 < px2 :
     pp_min = px1; pp_max = px2; 
else:
     pp_min = px2; pp_max = px1; 
Ep_min = sqrt(1 + pp_min*pp_min); Ep_max = sqrt(1 + pp_max*pp_max)

#  Pick the positron energy 
eta_Ep = randomset() if( eta_Ep < tiny_eta ) eta_Ep == tiny_eta
wEp = Ep_max - Ep_min; Ep = Ep_min + eta_Ep*wEp
pp2 = Ep*Ep - 1; pp = sqrt(pp2); k2p2 = k*k + pp2

#  Now we can determine the pair electron energy from energy conservation 
Em = k + 1 - Er - Ep
pm2 = Em*Em-1; pm = sqrt(pm2)

#  The minimum cosine of the positron angle follows from the kinematics. 
mup_min = (k2p2 - (pr + pm)*(pr + pm))/(2*k*pp)

#  Now pick the positron direction from 1/(Ep-pp*cost_p) 
eta_costp = randomset() if( eta_costp < tiny_eta ) eta_costp == tiny_eta
Epp = Ep/pp; wmup = log((Epp-1)/(Epp-mup_min))
cost_p = Epp - (Epp - mup_min)*exp(wmup*eta_costp)
wmup = wmup*(cost_p - Epp)
sint_p = 1-cost_p*cost_p
if  sint_p > 1e-20 :
     sint_p = sqrt(sint_p);  ELSE [ sint_p = 1e-10; ]
k2p2x = k2p2 - 2*k*pp*cost_p

#  The minimum amd maximum cosine of the pair electron angle follows from 
#  the kinematics 
b = pr2-k2p2x-pm2; aux1 = k - pp*cost_p; aux12 = aux1*aux1
pp_sintp = pp*sint_p; pp_sntp2 = pp_sintp*pp_sintp
D1 = pm2*(aux12+pp_sntp2)-b*b/4
if( D1 <= 0 ) [ goto_retry_triplet  is True
 break # XXX; ]
D = 2*pp_sintp*sqrt(D1)
aux3 = 0.5/(aux12+pp_sntp2)
xmin = (-b*aux1-D)*aux3; xmax = (-b*aux1+D)*aux3

#  Now pick the electron direction from 
#   1/(Em-pm*cost_m)/sqrt((cost_m_max-cost_m)*(cost_m-cost_m_min)) 
#  We have to take into account the 
#  1/sqrt((cost_m_max-cost_m)*(cost_m-cost_m_min)) factor in the sampling 
#  otherwise we end up with 1/sqrt() singularities near the ends of the 
#  allowed cost_m range                                                 
eta_costm = randomset() if( eta_costm < tiny_eta ) eta_costm == tiny_eta
aux6 = sqrt((Em-xmin)/(Em-xmax))
aux7 = aux6*tan(1.570796326794897*eta_costm)
uu = (aux7-1)/(aux7+1)
cost_m = 0.5*(xmax + xmin + 2*uu*(xmax-xmin)/(1+uu*uu))
wmum = sqrt((xmax-cost_m)*(cost_m-xmin))
wmum = wmum*aux6*(Em-cost_m)/(Em-xmin)
cost_m = cost_m/pm
sint_m = sqrt(1-cost_m*cost_m); pm_sintm = pm*sint_m

#  Now we have selected all independent kinematic variables. 
#  Determine the azimuthal angle between the pair electrons 
cphi = (b + 2*pm*cost_m*aux1)/(2*pp_sintp*pm_sintm)
if  abs(cphi) >= 1 :
     goto_retry_triplet = True
     break # XXX; 
sphi = sqrt(1-cphi*cphi)

#  And now evaluate the Borsellino cross section 
k3 = k*(pp*cost_p - Ep); k2 = k*(pm*cost_m - Em)
k22 = k2*k2; k32 = k3*k3; k2k3i = 1/(k2*k3)
s2 = pp*pm*(cost_p*cost_m + sint_p*sint_m*cphi) - Ep*Em
s3 = k2 - Em + 1 - s2; q2 = 2*(Er-1)
S_1 = k32+k22+(q2-2)*s2-(1-q2/2)*(k32+k22)*k2k3i
aux4 = k3*Ep-k2*Em
S_2 = -q2*(Ep*Ep+Em*Em) + 2*s2 - (2*aux4*aux4 - k22 - k32)*k2k3i
sigma = abs(pp*pm2*pm*k2k3i/(q2*q2*(Em*s3+Er))*(S_1*(1-q2/4)+S_2*(1+q2/4)))

#  We get the following factor due to the transformation from phi to 
#  the recoil momentum pr 
cphi_factor = abs(2*Er*pm2-Em*(k2p2x-pr2-pm2))/(2*pp_sintp*pm_sintm*pm2*sphi)

#  We have to also multiply by the various factors from the sampling of 
#  pr, Ep, cost_p and cost_m 
sigma = sigma*cphi_factor*wEp*wmup*wmum*wp_array(i)*pr2/Er
if  sigma < 0 :

    logger.warning('***************** Warning: 
    'In triplet sigma < 0 ? ',sigma')

#  Now determine if we accept this new event 
use_it = True
if  sigma < fmax_array(i) :

    rnno = randomset()
    if  sigma < fmax_array(i)*rnno :
         use_it = False 

if  use_it :
     [       #  Yes, event accepted 
    fmax_array(i) = sigma
    eta_p_array(i) = eta_pr; eta_Ep_array(i) = eta_Ep
    eta_costp_array(i) = eta_costp; eta_costm_array(i) = eta_costm
else:
    eta_pr = eta_p_array(i); eta_Ep = eta_Ep_array(i)
    eta_costp = eta_costp_array(i); eta_costm = eta_costm_array(i)

#  We now have a set of random number accepted for sampling around 
#  the i'th bin energy. We need to recalculate all variables using 
#  the actual photon energy 

k = peig*prmi
aux5 = k*(k-1)+(k+1)*sqrt(k*(k-4))
qmin = 4*k/aux5; qmax = aux5/(2*k+1)
pr = qmin*exp(eta_pr*log(qmax/qmin))
pr2 = pr*pr; Er = sqrt(1+pr2)

aux = Er-pr-1; a1=(k-pr)*(1-Er-k*aux); a2=1+k-Er; a3=1/(aux*(pr+Er-2*k-1))
D = a2*sqrt(aux*(2*k*Er+k*k*aux-pr*(Er+pr+1)/2))
px1 = (a1 + D)*a3; px2 = (a1 - D)*a3
if  px1 < px2 :
     pp_min = px1; pp_max = px2; 
else:
     pp_min = px2; pp_max = px1; 
Ep_min = sqrt(1 + pp_min*pp_min); Ep_max = sqrt(1 + pp_max*pp_max)

wEp = Ep_max - Ep_min; Ep = Ep_min + eta_Ep*wEp
pp2 = Ep*Ep - 1; pp = sqrt(pp2); k2p2 = k*k + pp2
Em = k + 1 - Er - Ep
pm2 = Em*Em-1; pm = sqrt(pm2)

mup_min = (k2p2 - (pr + pm)*(pr + pm))/(2*k*pp)
Epp = Ep/pp; wmup = log((Epp-1)/(Epp-mup_min))
cost_p = Epp - (Epp - mup_min)*exp(wmup*eta_costp)
sint_p = sqrt(1-cost_p*cost_p)
k2p2x = k2p2 - 2*k*pp*cost_p

b = pr2-k2p2x-pm2; aux1 = k - pp*cost_p; aux12 = aux1*aux1
pp_sintp = pp*sint_p; pp_sntp2 = pp_sintp*pp_sintp
D1 = pm2*(aux12+pp_sntp2)-b*b/4
if( D1 <= 0 ) [ goto_retry_triplet  is True
 break # XXX; ]
D = 2*pp_sintp*sqrt(D1)
aux3 = 0.5/(aux12+pp_sntp2)
xmin = (-b*aux1-D)*aux3; xmax = (-b*aux1+D)*aux3
aux6 = sqrt((Em-xmin)/(Em-xmax))
aux7 = aux6*tan(1.570796326794897*eta_costm)
uu = (aux7-1)/(aux7+1)
cost_m = 0.5*(xmax + xmin + 2*uu*(xmax-xmin)/(1+uu*uu))/pm
sint_m = sqrt(1-cost_m*cost_m); pm_sintm = pm*sint_m

cphi = (b + 2*pm*cost_m*aux1)/(2*pp_sintp*pm_sintm)
if  abs(cphi) >= 1 :
     goto_retry_triplet = True
     break # XXX; 
sphi = sqrt(1-cphi*cphi)

/*
   OK, now the final momenta are
     Positron:     (pp*sint_p,      0,             pp*cost_p)
     Electron:     (pm*sint_m*cphi,pm*sint_m*sphi, pm*cost_m)
 Recoil electron:  k - pp - pm
   This is in a frame where the photon is moving along the z axis.
   We have to pick another azimuthal angle randomly, rotate the
   x- and y- components of pp and pm by that, determine the recoil
   momentum from momentum conservation and then rotate all three
   momenta back into the lab frame.
*/
phi = randomset() phi = phi*6.283185307179586
ppx = pp*sint_p; ppy = 0
pmx = pm*sint_m*cphi; pmy = pm*sint_m*sphi
cphi = cos(phi); sphi = sin(phi)
tmp = ppx*sphi; ppx = ppx*cphi - ppy*sphi; ppy = tmp + ppy*cphi
tmp = pmx*sphi; pmx = pmx*cphi - pmy*sphi; pmy = tmp + pmy*cphi
ppz = pp*cost_p; pmz = pm*cost_m
prx = -ppx-pmx; pry = -ppy-pmy; prz = k - ppz - pmz
#  Set up particles on the stack 
#  We always put the recoil electron on top (even if its energy is higher 
#  then the energies of the pair particles) because                       
#    - that way, we know which particle is the recoil  electron in case   
#      we want to score some quantity related to it                       
#    - its energy is, on average, lower than the pair particle energies   
NPold = np
$TRANSFER PROPERTIES TO (np)   FROM (np)
$TRANSFER PROPERTIES TO (np+1) FROM (np)
$TRANSFER PROPERTIES TO (np+2) FROM (np+1)
pp = 1/pp; pm = 1/pm; pr = 1/pr
a = u[np]; b = v[np]; c = w[np]; sinpsi = a*a + b*b
if  sinpsi > 1e-20 :

    sinpsi = sqrt(sinpsi); sindel = b/sinpsi; cosdel = a/sinpsi
    if  Ep > Em :

        u[np]   = pp*(c*cosdel*ppx - sindel*ppy + a*ppz)
        v[np]   = pp*(c*sindel*ppx + cosdel*ppy + b*ppz)
        w[np]   = pp*(c*ppz - sinpsi*ppx); iq[np] = 1; e[np] = Ep*prm
        u(np+1) = pm*(c*cosdel*pmx - sindel*pmy + a*pmz)
        v(np+1) = pm*(c*sindel*pmx + cosdel*pmy + b*pmz)
        w(np+1) = pm*(c*pmz - sinpsi*pmx); iq(np+1) = -1; E(np+1) = Em*prm
    else:
        u(np+1) = pp*(c*cosdel*ppx - sindel*ppy + a*ppz)
        v(np+1) = pp*(c*sindel*ppx + cosdel*ppy + b*ppz)
        w(np+1) = pp*(c*ppz - sinpsi*ppx); iq(np+1) = 1; E(np+1) = Ep*prm
        u[np]   = pm*(c*cosdel*pmx - sindel*pmy + a*pmz)
        v[np]   = pm*(c*sindel*pmx + cosdel*pmy + b*pmz)
        w[np]   = pm*(c*pmz - sinpsi*pmx); iq[np] = -1; e[np] = Em*prm

    np = np + 2
    u[np]   = pr*(c*cosdel*prx - sindel*pry + a*prz)
    v[np]   = pr*(c*sindel*prx + cosdel*pry + b*prz)
    w[np]   = pr*(c*prz - sinpsi*prx); iq[np] = -1; e[np] = Er*prm
else:
    if  Ep > Em :

        u[np] = pp*ppx; v[np] = pp*ppy; w[np] = c*pp*ppz
        iq[np] = 1; e[np] = Ep*prm
        u(np+1) = pm*pmx; v(np+1) = pm*pmy; w(np+1) = c*pm*pmz
        iq(np+1) = -1; E(np+1) = Em*prm
    else:
        u(np+1) = pp*ppx; v(np+1) = pp*ppy; w(np+1) = c*pp*ppz
        iq(np+1) = 1; E(np+1) = Ep*prm
        u[np] = pm*pmx; v[np] = pm*pmy; w[np] = c*pm*pmz
        iq[np] = -1; e[np] = Em*prm

    np = np + 2
    u[np] = pr*prx; v[np] = pr*pry; w[np] = c*pr*prz
    iq[np] = -1; e[np] = Er*prm

return; end