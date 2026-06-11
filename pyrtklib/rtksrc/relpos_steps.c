/*------------------------------------------------------------------------------
* relpos_steps.c : step-by-step relative positioning (exposing internals)
*
*          Wrapper around RTKLIB's relpos function, broken into individual steps
*          so that Python code can inspect intermediate state (float states,
*          fixed ambiguities, DD residuals, covariances, etc.) at each stage.
*
*-----------------------------------------------------------------------------*/
#include "rtklib.h"
#include <string.h>

/* macros from rtkpos.c (duplicated here since they are file-local there) */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics==0?3:9)
#define NI(opt)     ((opt)->ionoopt!=IONOOPT_EST?0:MAXSAT)
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt<TROPOPT_ESTG?2:6))
#define NL(opt)     ((opt)->glomodear!=2?0:NFREQGLO)
#define NB(opt)     ((opt)->mode<=PMODE_DGPS?0:MAXSAT*NF(opt))
#define NR_OPT(opt) (NP(opt)+NI(opt)+NT(opt)+NL(opt))
#define NX(opt)     (NR_OPT(opt)+NB(opt))
#define IB(s,f,opt) (NR_OPT(opt)+MAXSAT*(f)+(s)-1)

#define TTOL_MOVEB  (1.0+2*DTTOL)

/* time-interpolation of residuals (local copy for step-by-step API) ---------
*  This is a separate copy from the one in rtkpos.c. Each copy owns its own
*  static local buffers so the original relpos() path and the step-by-step
*  path do not corrupt each other's interpolation state.
*---------------------------------------------------------------------------*/
static double intpres_step(gtime_t time, const obsd_t *obs, int n,
                           const nav_t *nav, rtk_t *rtk, double *y)
{
    static obsd_t obsb[MAXOBS];
    static double yb[MAXOBS*NFREQ*2],rs[MAXOBS*6],dts[MAXOBS*2],var[MAXOBS];
    static double e[MAXOBS*3],azel[MAXOBS*2],freq[MAXOBS*NFREQ];
    static int nb=0,svh[MAXOBS*2];
    prcopt_t *opt=&rtk->opt;
    double tt=timediff(time,obs[0].time),ttb,*p,*q;
    int i,j,k,nf=NF(opt);

    if (nb==0||fabs(tt)<DTTOL) {
        nb=n; for (i=0;i<n;i++) obsb[i]=obs[i];
        return tt;
    }
    ttb=timediff(time,obsb[0].time);
    if (fabs(ttb)>opt->maxtdiff*2.0||ttb==tt) return tt;

    satposs(time,obsb,nb,nav,opt->sateph,rs,dts,var,svh);

    if (!zdres(1,obsb,nb,rs,dts,var,svh,nav,rtk->rb,opt,1,yb,e,azel,freq)) {
        return tt;
    }
    for (i=0;i<n;i++) {
        for (j=0;j<nb;j++) if (obsb[j].sat==obs[i].sat) break;
        if (j>=nb) continue;
        for (k=0,p=y+i*nf*2,q=yb+j*nf*2;k<nf*2;k++,p++,q++) {
            if (*p==0.0||*q==0.0) *p=0.0; else *p=(ttb*(*p)-tt*(*q))/(ttb-tt);
        }
    }
    return fabs(ttb)>fabs(tt)?ttb:tt;
}

/* pre-relpos processing (from rtkpos): base setup, SPP, time sync ----------
*  Does everything rtkpos() does before calling relpos().
*  Returns: 1 on success (ready for relpos steps), 0 on failure.
*---------------------------------------------------------------------------*/
int rtkpos_pre_relpos(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    prcopt_t *opt = &rtk->opt;
    sol_t solb = {{0}};
    gtime_t prev_time;
    int i, nu, nr;
    char msg[128] = "";

    /* set base station position */
    if (opt->refpos <= POSOPT_RINEX && opt->mode != PMODE_SINGLE &&
        opt->mode != PMODE_MOVEB) {
        for (i = 0; i < 6; i++) rtk->rb[i] = i < 3 ? opt->rb[i] : 0.0;
    }

    /* count rover/base observations */
    for (nu = 0; nu < n && obs[nu].rcv == 1; nu++) ;
    for (nr = 0; nu + nr < n && obs[nu + nr].rcv == 2; nr++) ;

    prev_time = rtk->sol.time;

    /* SPP for rover */
    if (!pntpos(obs, nu, nav, &rtk->opt, &rtk->sol, NULL, rtk->ssat, msg)) {
        if (!rtk->opt.dynamics) {
            return 0;
        }
    }
    if (prev_time.time != 0) rtk->tt = timediff(rtk->sol.time, prev_time);

    if (opt->mode == PMODE_SINGLE) return 1;
    if (!opt->outsingle) rtk->sol.stat = SOLQ_NONE;
    if (opt->mode >= PMODE_PPP_KINEMA) return 0; /* not supported here */

    if (nr == 0) return 0;

    if (opt->mode == PMODE_MOVEB) {
        if (!pntpos(obs + nu, nr, nav, &rtk->opt, &solb, NULL, NULL, msg)) {
            return 0;
        }
        rtk->sol.age = (float)timediff(rtk->sol.time, solb.time);
        if (fabs(rtk->sol.age) > TTOL_MOVEB) return 0;
        for (i = 0; i < 6; i++) rtk->rb[i] = solb.rr[i];
        for (i = 0; i < 3; i++) rtk->rb[i] += rtk->rb[i + 3] * rtk->sol.age;
    }
    else {
        rtk->sol.age = (float)timediff(obs[0].time, obs[nu].time);
        if (fabs(rtk->sol.age) > opt->maxtdiff) return 0;
    }
    return 1;
}

/* initialize relpos context -------------------------------------------------
*  Sets up all working arrays and counts rover/base observations.
*  Call once per epoch before the other step functions.
*  Returns 1 on success, 0 on failure.
*---------------------------------------------------------------------------*/
int relpos_init(relpos_ctx_t *ctx, rtk_t *rtk, const obsd_t *obs, int n,
                const nav_t *nav)
{
    prcopt_t *opt = &rtk->opt;
    int i, j, nu, nr;

    memset(ctx, 0, sizeof(relpos_ctx_t));
    ctx->rtk = rtk;
    ctx->nav = nav;

    ctx->nf = opt->ionoopt == IONOOPT_IFLC ? 1 : opt->nf;
    ctx->stat = opt->mode <= PMODE_DGPS ? SOLQ_DGPS : SOLQ_FLOAT;

    /* count rover/base observations */
    for (nu = 0; nu < n && obs[nu].rcv == 1; nu++) ;
    for (nr = 0; nu + nr < n && obs[nu + nr].rcv == 2; nr++) ;
    ctx->nu = nu;
    ctx->nr = nr;

    n = nu + nr;

    /* allocate working arrays */
    ctx->rs   = mat(6, n);
    ctx->dts  = mat(2, n);
    ctx->var  = mat(1, n);
    ctx->y    = mat(ctx->nf * 2, n);
    ctx->e    = mat(3, n);
    ctx->azel = zeros(2, n);
    ctx->freq = zeros(ctx->nf, n);

    ctx->ns = 0;
    ctx->nv = 0;
    ctx->ny = 0;

    /* reset satellite status */
    for (i = 0; i < MAXSAT; i++) {
        rtk->ssat[i].sys = satsys(i + 1, NULL);
        for (j = 0; j < NFREQ; j++) rtk->ssat[i].vsat[j] = 0;
        for (j = 1; j < NFREQ; j++) rtk->ssat[i].snr[j]  = 0;
    }

    /* time difference */
    if (nu > 0 && nr > 0) {
        ctx->dt = timediff(obs[0].time, obs[nu].time);
    }

    /* xp, Pp, xa, bias will be allocated after selsat when we know ns */
    ctx->xp = NULL;
    ctx->Pp = NULL;
    ctx->xa = NULL;
    ctx->bias = NULL;
    ctx->v = NULL;
    ctx->H = NULL;
    ctx->R = NULL;

    return 1;
}

/* step 1: compute satellite positions/clocks --------------------------------*/
int relpos_satpos(relpos_ctx_t *ctx, const obsd_t *obs)
{
    int n = ctx->nu + ctx->nr;
    prcopt_t *opt = &ctx->rtk->opt;

    satposs(obs[0].time, obs, n, ctx->nav, opt->sateph,
            ctx->rs, ctx->dts, ctx->var, ctx->svh);

    return 1;
}

/* step 2: undifferenced residuals for base station --------------------------*/
int relpos_zdres_base(relpos_ctx_t *ctx, const obsd_t *obs)
{
    rtk_t *rtk = ctx->rtk;
    prcopt_t *opt = &rtk->opt;
    int nu = ctx->nu, nr = ctx->nr;

    if (!zdres(1, obs + nu, nr, ctx->rs + nu * 6, ctx->dts + nu * 2,
               ctx->var + nu, ctx->svh + nu, ctx->nav, rtk->rb, opt, 1,
               ctx->y + nu * ctx->nf * 2, ctx->e + nu * 3,
               ctx->azel + nu * 2, ctx->freq + nu * ctx->nf)) {
        return 0;
    }
    /* time-interpolation of residuals (for post-processing) */
    if (opt->intpref) {
        ctx->dt = intpres_step(obs[0].time, obs + nu, nr, ctx->nav, rtk,
                               ctx->y + nu * ctx->nf * 2);
    }
    return 1;
}

/* step 3: select common satellites between rover and base -------------------*/
int relpos_selsat(relpos_ctx_t *ctx, const obsd_t *obs)
{
    rtk_t *rtk = ctx->rtk;
    prcopt_t *opt = &rtk->opt;

    ctx->ns = selsat(obs, ctx->azel, ctx->nu, ctx->nr, opt,
                     ctx->sat, ctx->iu, ctx->ir);
    if (ctx->ns <= 0) return 0;

    /* now allocate filter working arrays */
    ctx->xp   = mat(rtk->nx, 1);
    ctx->Pp   = zeros(rtk->nx, rtk->nx);
    ctx->xa   = mat(rtk->nx, 1);
    ctx->bias = mat(rtk->nx, 1);

    ctx->ny = ctx->ns * ctx->nf * 2 + 2;
    ctx->v = mat(ctx->ny, 1);
    ctx->H = zeros(rtk->nx, ctx->ny);
    ctx->R = mat(ctx->ny, ctx->ny);

    return ctx->ns;
}

/* step 4: temporal update of states -----------------------------------------*/
void relpos_udstate(relpos_ctx_t *ctx, const obsd_t *obs)
{
    rtk_t *rtk = ctx->rtk;
    prcopt_t *opt = &rtk->opt;

    udstate(rtk, obs, ctx->sat, ctx->iu, ctx->ir, ctx->ns, ctx->nav);

    /* copy initial state for iterative filter */
    matcpy(ctx->xp, rtk->x, rtk->nx, 1);

    ctx->niter = opt->niter +
                 (opt->mode == PMODE_MOVEB && opt->baseline[0] > 0.0 ? 2 : 0);
}

/* step 5: iterative Kalman filter (float solution) --------------------------
*  Runs the KF measurement update iterations.
*  After this, xp/Pp contain the float solution.
*  Post-fit DD residuals are computed and the float solution is validated.
*  Returns: 1 if float solution is valid, 0 otherwise.
*---------------------------------------------------------------------------*/
int relpos_float_filter(relpos_ctx_t *ctx, const obsd_t *obs)
{
    rtk_t *rtk = ctx->rtk;
    prcopt_t *opt = &rtk->opt;
    int i, f, nv, info;

    for (i = 0; i < ctx->niter; i++) {
        /* UD residuals for rover */
        if (!zdres(0, obs, ctx->nu, ctx->rs, ctx->dts, ctx->var, ctx->svh,
                   ctx->nav, ctx->xp, opt, 0, ctx->y, ctx->e, ctx->azel,
                   ctx->freq)) {
            ctx->stat = SOLQ_NONE;
            return 0;
        }
        /* DD residuals and partial derivatives */
        nv = ddres(rtk, ctx->nav, ctx->dt, ctx->xp, ctx->Pp, ctx->sat,
                   ctx->y, ctx->e, ctx->azel, ctx->freq, ctx->iu, ctx->ir,
                   ctx->ns, ctx->v, ctx->H, ctx->R, ctx->vflg);
        if (nv < 1) {
            ctx->stat = SOLQ_NONE;
            return 0;
        }
        /* Kalman filter measurement update */
        matcpy(ctx->Pp, rtk->P, rtk->nx, rtk->nx);
        info = filter(ctx->xp, ctx->Pp, ctx->H, ctx->v, ctx->R, rtk->nx, nv);
        if (info) {
            ctx->stat = SOLQ_NONE;
            return 0;
        }
    }

    /* compute post-fit residuals for float solution */
    if (!zdres(0, obs, ctx->nu, ctx->rs, ctx->dts, ctx->var, ctx->svh,
               ctx->nav, ctx->xp, opt, 0, ctx->y, ctx->e, ctx->azel,
               ctx->freq)) {
        ctx->stat = SOLQ_NONE;
        return 0;
    }
    ctx->nv = ddres(rtk, ctx->nav, ctx->dt, ctx->xp, ctx->Pp, ctx->sat,
                    ctx->y, ctx->e, ctx->azel, ctx->freq, ctx->iu, ctx->ir,
                    ctx->ns, ctx->v, NULL, ctx->R, ctx->vflg);

    /* validate float solution */
    if (!valpos(rtk, ctx->v, ctx->R, ctx->vflg, ctx->nv, 4.0)) {
        ctx->stat = SOLQ_NONE;
        return 0;
    }

    /* update state and covariance */
    matcpy(rtk->x, ctx->xp, rtk->nx, 1);
    matcpy(rtk->P, ctx->Pp, rtk->nx, rtk->nx);

    /* update ambiguity control */
    rtk->sol.ns = 0;
    for (i = 0; i < ctx->ns; i++) for (f = 0; f < ctx->nf; f++) {
        if (!rtk->ssat[ctx->sat[i] - 1].vsat[f]) continue;
        rtk->ssat[ctx->sat[i] - 1].lock[f]++;
        rtk->ssat[ctx->sat[i] - 1].outc[f] = 0;
        if (f == 0) rtk->sol.ns++;
    }
    if (rtk->sol.ns < 4) {
        ctx->stat = SOLQ_NONE;
        return 0;
    }
    return 1;
}

/* step 6: integer ambiguity resolution (LAMBDA) ----------------------------
*  Attempts LAMBDA ambiguity resolution.
*  On success, xa[] contains the fixed state vector (position + resolved
*  ambiguities restored to single-difference form).
*  Returns: number of resolved DD ambiguities (>1 on success, 0 on failure).
*---------------------------------------------------------------------------*/
int relpos_ambiguity_resolution(relpos_ctx_t *ctx, const obsd_t *obs)
{
    rtk_t *rtk = ctx->rtk;
    prcopt_t *opt = &rtk->opt;
    int nb, nv;

    nb = resamb_LAMBDA(rtk, ctx->bias, ctx->xa);
    if (nb <= 1) return 0;

    /* post-fit residuals for fixed solution */
    if (!zdres(0, obs, ctx->nu, ctx->rs, ctx->dts, ctx->var, ctx->svh,
               ctx->nav, ctx->xa, opt, 0, ctx->y, ctx->e, ctx->azel,
               ctx->freq)) {
        return 0;
    }
    nv = ddres(rtk, ctx->nav, ctx->dt, ctx->xa, NULL, ctx->sat, ctx->y,
               ctx->e, ctx->azel, ctx->freq, ctx->iu, ctx->ir, ctx->ns,
               ctx->v, NULL, ctx->R, ctx->vflg);

    if (!valpos(rtk, ctx->v, ctx->R, ctx->vflg, nv, 4.0)) {
        return 0;
    }

    /* hold integer ambiguity */
    if (++rtk->nfix >= rtk->opt.minfix &&
        rtk->opt.modear == ARMODE_FIXHOLD) {
        holdamb(rtk, ctx->xa);
    }
    ctx->stat = SOLQ_FIX;
    return nb;
}

/* step 7: save solution status and clean up per-epoch bookkeeping ----------*/
void relpos_save_solution(relpos_ctx_t *ctx, const obsd_t *obs)
{
    rtk_t *rtk = ctx->rtk;
    int i, j, n = ctx->nu + ctx->nr;
    int nf = ctx->nf;
    int stat = ctx->stat;

    /* save solution status */
    if (stat == SOLQ_FIX) {
        for (i = 0; i < 3; i++) {
            rtk->sol.rr[i] = rtk->xa[i];
            rtk->sol.qr[i] = (float)rtk->Pa[i + i * rtk->na];
        }
        rtk->sol.qr[3] = (float)rtk->Pa[1];
        rtk->sol.qr[4] = (float)rtk->Pa[1 + 2 * rtk->na];
        rtk->sol.qr[5] = (float)rtk->Pa[2];

        if (rtk->opt.dynamics) {
            for (i = 3; i < 6; i++) {
                rtk->sol.rr[i] = rtk->xa[i];
                rtk->sol.qv[i - 3] = (float)rtk->Pa[i + i * rtk->na];
            }
            rtk->sol.qv[3] = (float)rtk->Pa[4 + 3 * rtk->na];
            rtk->sol.qv[4] = (float)rtk->Pa[5 + 4 * rtk->na];
            rtk->sol.qv[5] = (float)rtk->Pa[5 + 3 * rtk->na];
        }
    }
    else {
        for (i = 0; i < 3; i++) {
            rtk->sol.rr[i] = rtk->x[i];
            rtk->sol.qr[i] = (float)rtk->P[i + i * rtk->nx];
        }
        rtk->sol.qr[3] = (float)rtk->P[1];
        rtk->sol.qr[4] = (float)rtk->P[1 + 2 * rtk->nx];
        rtk->sol.qr[5] = (float)rtk->P[2];

        if (rtk->opt.dynamics) {
            for (i = 3; i < 6; i++) {
                rtk->sol.rr[i] = rtk->x[i];
                rtk->sol.qv[i - 3] = (float)rtk->P[i + i * rtk->nx];
            }
            rtk->sol.qv[3] = (float)rtk->P[4 + 3 * rtk->nx];
            rtk->sol.qv[4] = (float)rtk->P[5 + 4 * rtk->nx];
            rtk->sol.qv[5] = (float)rtk->P[5 + 3 * rtk->nx];
        }
        rtk->nfix = 0;
    }
    /* carrier-phase history */
    for (i = 0; i < n; i++) for (j = 0; j < nf; j++) {
        if (obs[i].L[j] == 0.0) continue;
        rtk->ssat[obs[i].sat - 1].pt[obs[i].rcv - 1][j] = obs[i].time;
        rtk->ssat[obs[i].sat - 1].ph[obs[i].rcv - 1][j] = obs[i].L[j];
    }
    /* SNR */
    for (i = 0; i < ctx->ns; i++) for (j = 0; j < nf; j++) {
        rtk->ssat[ctx->sat[i] - 1].snr[j] = obs[ctx->iu[i]].SNR[j];
    }
    /* fix/slip flags */
    for (i = 0; i < MAXSAT; i++) for (j = 0; j < nf; j++) {
        if (rtk->ssat[i].fix[j] == 2 && stat != SOLQ_FIX)
            rtk->ssat[i].fix[j] = 1;
        if (rtk->ssat[i].slip[j] & 1) rtk->ssat[i].slipc[j]++;
    }

    if (stat != SOLQ_NONE) rtk->sol.stat = stat;
}

/* bulk per-satellite data extraction ----------------------------------------
*  Extracts per-satellite data for all common satellites in one C call.
*  Computes geometric range, tropospheric/ionospheric corrections, and
*  base station geometry for all satellites, avoiding Python-level loops.
*
*  Output arrays must be pre-allocated by the caller to at least ns elements
*  (or ns*nf / ns*3 / ns*6 depending on the field).
*
*  flags bitmask controls which fields are computed:
*    bit 0 (1)  : VRS corrections (geodist, tropo, iono, sagnac)
*    bit 1 (2)  : base station geometry (base geodist, base azel)
*    bit 2 (4)  : per-frequency ssat fields (fix, lock, slip, snr, resc, resp)
*    bit 3 (8)  : float ambiguities + wavelengths
*  All bits set (0xF) = extract everything.
*---------------------------------------------------------------------------*/
void relpos_extract_sat_data(
    const relpos_ctx_t *ctx,
    const double *rover_ecef,    /* rover ECEF position [3] */
    const double *base_ecef,     /* base station ECEF [3] */
    int flags,                   /* bitmask controlling which fields to compute */
    /* outputs: per-satellite arrays (all ns-sized or ns*nf / ns*3 etc.) */
    double *out_el_deg,          /* [ns] elevation (deg) */
    double *out_az_deg,          /* [ns] azimuth (deg) */
    double *out_sat_pos,         /* [ns*3] satellite ECEF positions */
    double *out_sat_vel,         /* [ns*3] satellite ECEF velocities */
    double *out_sat_clk,         /* [ns] satellite clock bias (m) */
    double *out_sat_clk_drift,   /* [ns] satellite clock drift (m/s) */
    double *out_los,             /* [ns*3] line-of-sight unit vectors */
    double *out_geom_range,      /* [ns] geometric range (m) (flag bit 0) */
    double *out_sagnac,          /* [ns] Sagnac correction (m) (flag bit 0) */
    double *out_tropo,           /* [ns] tropospheric delay (m) (flag bit 0) */
    double *out_iono,            /* [ns] ionospheric delay L1 (m) (flag bit 0) */
    double *out_phw,             /* [ns] phase wind-up (cycles) */
    double *out_base_geom_range, /* [ns] base geometric range (m) (flag bit 1) */
    double *out_base_el_deg,     /* [ns] base elevation (deg) (flag bit 1) */
    double *out_base_az_deg,     /* [ns] base azimuth (deg) (flag bit 1) */
    double *out_float_amb,       /* [ns*nf] float ambiguities (flag bit 3) */
    double *out_wl,              /* [ns*nf] wavelengths (m) (flag bit 3) */
    double *out_resc,            /* [ns*nf] DD phase residuals (flag bit 2) */
    double *out_resp,            /* [ns*nf] DD code residuals (flag bit 2) */
    double *out_fix,             /* [ns*nf] fix flags (flag bit 2) */
    double *out_lock,            /* [ns*nf] lock counts (flag bit 2) */
    double *out_slip,            /* [ns*nf] slip flags (flag bit 2) */
    double *out_snr,             /* [ns*nf] SNR (dBHz) (flag bit 2) */
    double *out_rover_dant,      /* [ns*nf] rover receiver antenna corr (m) (flag bit 0) */
    double *out_base_dant,       /* [ns*nf] base receiver antenna corr (m) (flag bit 1) */
    double *out_base_tropo,      /* [ns] base tropospheric delay (m) (flag bit 1) */
    double *out_base_iono        /* [ns] base ionospheric delay L1 (m) (flag bit 1) */
)
{
    rtk_t *rtk = ctx->rtk;
    prcopt_t *opt = &rtk->opt;
    int j, f, sat_no, iu_j, ir_j, nx = rtk->nx;
    int ns = ctx->ns, nf = ctx->nf;
    double rover_pos[3], base_pos[3];  /* geodetic (lat,lon,h) */
    double scratch_sv[6], scratch_e[3], azel_buf[2];
    double trp[1], trp_var[1], ion[1], ion_var[1];
    double btrp[1], btrp_var[1], bion[1], bion_var[1];
    int do_vrs   = (flags & 1);
    int do_base  = (flags & 2);
    int do_ssat  = (flags & 4);
    int do_amb   = (flags & 8);

    /* pre-compute geodetic positions for iono/tropo models */
    if (do_vrs) {
        ecef2pos(rover_ecef, rover_pos);
    }
    if (do_base) {
        ecef2pos(base_ecef, base_pos);
    }

    for (j = 0; j < ns; j++) {
        sat_no = ctx->sat[j];
        if (sat_no <= 0) continue;

        iu_j = ctx->iu[j];
        ir_j = ctx->ir[j];

        /* elevation and azimuth from ssat (already set by zdres) */
        out_el_deg[j] = rtk->ssat[sat_no - 1].azel[1] * R2D;
        out_az_deg[j] = rtk->ssat[sat_no - 1].azel[0] * R2D;

        /* satellite position/velocity from ctx->rs */
        {
            int off6 = 6 * iu_j;
            out_sat_pos[j*3+0] = ctx->rs[off6+0];
            out_sat_pos[j*3+1] = ctx->rs[off6+1];
            out_sat_pos[j*3+2] = ctx->rs[off6+2];
            out_sat_vel[j*3+0] = ctx->rs[off6+3];
            out_sat_vel[j*3+1] = ctx->rs[off6+4];
            out_sat_vel[j*3+2] = ctx->rs[off6+5];
        }

        /* satellite clock bias/drift */
        {
            int off2 = 2 * iu_j;
            out_sat_clk[j]       = -CLIGHT * ctx->dts[off2+0];
            out_sat_clk_drift[j] = -CLIGHT * ctx->dts[off2+1];
        }

        /* LOS unit vector from ctx->e */
        {
            int off3 = 3 * iu_j;
            out_los[j*3+0] = ctx->e[off3+0];
            out_los[j*3+1] = ctx->e[off3+1];
            out_los[j*3+2] = ctx->e[off3+2];
        }

        /* phase wind-up */
        out_phw[j] = rtk->ssat[sat_no - 1].phw;

        /* VRS corrections: geodist, sagnac, tropo, iono */
        if (do_vrs) {
            int off6 = 6 * iu_j;
            for (f = 0; f < 6; f++) scratch_sv[f] = ctx->rs[off6+f];

            out_geom_range[j] = geodist(scratch_sv, rover_ecef, scratch_e);

            out_sagnac[j] = OMGE * (scratch_sv[0] * rover_ecef[1]
                                   - scratch_sv[1] * rover_ecef[0]) / CLIGHT;

            azel_buf[0] = rtk->ssat[sat_no - 1].azel[0];
            azel_buf[1] = rtk->ssat[sat_no - 1].azel[1];

            trp[0] = 0.0;
            tropcorr(rtk->sol.time, ctx->nav, rover_pos, azel_buf,
                     opt->tropopt, trp, trp_var);
            out_tropo[j] = trp[0];

            ion[0] = 0.0;
            ionocorr(rtk->sol.time, ctx->nav, sat_no, rover_pos, azel_buf,
                     opt->ionoopt, ion, ion_var);
            out_iono[j] = ion[0];

            /* rover receiver antenna correction (PCO + PCV) */
            if (out_rover_dant) {
                double dant_buf[NFREQ] = {0};
                antmodel(opt->pcvr+0, opt->antdel[0], azel_buf, 1, dant_buf);
                for (f = 0; f < nf; f++) out_rover_dant[j*nf+f] = dant_buf[f];
            }
        }

        /* base station geometry */
        if (do_base) {
            int boff = 6 * ir_j;
            double base_scratch_e[3], base_azel[2];
            for (f = 0; f < 6; f++) scratch_sv[f] = ctx->rs[boff+f];

            out_base_geom_range[j] = geodist(scratch_sv, base_ecef, base_scratch_e);
            satazel(base_pos, base_scratch_e, base_azel);
            out_base_el_deg[j] = base_azel[1] * R2D;
            out_base_az_deg[j] = base_azel[0] * R2D;

            /* base tropospheric / ionospheric delays (same models as rover) */
            if (out_base_tropo) {
                btrp[0] = 0.0;
                tropcorr(rtk->sol.time, ctx->nav, base_pos, base_azel,
                         opt->tropopt, btrp, btrp_var);
                out_base_tropo[j] = btrp[0];
            }
            if (out_base_iono) {
                bion[0] = 0.0;
                ionocorr(rtk->sol.time, ctx->nav, sat_no, base_pos, base_azel,
                         opt->ionoopt, bion, bion_var);
                out_base_iono[j] = bion[0];
            }

            /* base receiver antenna correction (PCO + PCV) */
            if (out_base_dant) {
                double dant_buf[NFREQ] = {0};
                antmodel(opt->pcvr+1, opt->antdel[1], base_azel, 1, dant_buf);
                for (f = 0; f < nf; f++) out_base_dant[j*nf+f] = dant_buf[f];
            }
        }

        /* per-frequency ssat fields */
        if (do_ssat) {
            ssat_t *ss = &rtk->ssat[sat_no - 1];
            for (f = 0; f < nf; f++) {
                out_resc[j*nf+f] = ss->resc[f];
                out_resp[j*nf+f] = ss->resp[f];
                out_fix[j*nf+f]  = (double)ss->fix[f];
                out_lock[j*nf+f] = (double)ss->lock[f];
                out_slip[j*nf+f] = (double)ss->slip[f];
                out_snr[j*nf+f]  = (double)ss->snr[f] * SNR_UNIT;
            }
        }

        /* float ambiguities and wavelengths */
        if (do_amb) {
            for (f = 0; f < nf; f++) {
                int amb_idx = IB(sat_no, f, opt);
                out_float_amb[j*nf+f] = (amb_idx < nx) ? rtk->x[amb_idx] : 0.0/0.0;

                {
                    double freq_hz = ctx->freq[f + nf * iu_j];
                    out_wl[j*nf+f] = (freq_hz > 0.0) ? CLIGHT / freq_hz : 0.0/0.0;
                }
            }
        }
    }
}

/* bulk extraction of fixed-solution ambiguities ----------------------------*/
void relpos_extract_fixed_amb(
    const relpos_ctx_t *ctx,
    double *out_fixed_amb,  /* [ns*nf] fixed ambiguities from ctx->xa */
    double *out_fix_flags   /* [ns*nf] updated fix flags after AR */
)
{
    rtk_t *rtk = ctx->rtk;
    prcopt_t *opt = &rtk->opt;
    int j, f, sat_no, nx = rtk->nx;
    int ns = ctx->ns, nf = ctx->nf;

    if (!ctx->xa) return;

    for (j = 0; j < ns; j++) {
        sat_no = ctx->sat[j];
        if (sat_no <= 0) continue;
        for (f = 0; f < nf; f++) {
            int amb_idx = IB(sat_no, f, opt);
            out_fixed_amb[j*nf+f] = (amb_idx < nx) ? ctx->xa[amb_idx] : 0.0/0.0;
            out_fix_flags[j*nf+f] = (double)rtk->ssat[sat_no - 1].fix[f];
        }
    }
}

/* free relpos context -------------------------------------------------------*/
void relpos_free(relpos_ctx_t *ctx)
{
    if (ctx->rs)   { free(ctx->rs);   ctx->rs   = NULL; }
    if (ctx->dts)  { free(ctx->dts);  ctx->dts  = NULL; }
    if (ctx->var)  { free(ctx->var);  ctx->var  = NULL; }
    if (ctx->y)    { free(ctx->y);    ctx->y    = NULL; }
    if (ctx->e)    { free(ctx->e);    ctx->e    = NULL; }
    if (ctx->azel) { free(ctx->azel); ctx->azel = NULL; }
    if (ctx->freq) { free(ctx->freq); ctx->freq = NULL; }
    if (ctx->xp)   { free(ctx->xp);   ctx->xp   = NULL; }
    if (ctx->Pp)   { free(ctx->Pp);   ctx->Pp   = NULL; }
    if (ctx->xa)   { free(ctx->xa);   ctx->xa   = NULL; }
    if (ctx->bias) { free(ctx->bias); ctx->bias = NULL; }
    if (ctx->v)    { free(ctx->v);    ctx->v    = NULL; }
    if (ctx->H)    { free(ctx->H);    ctx->H    = NULL; }
    if (ctx->R)    { free(ctx->R);    ctx->R    = NULL; }
}
