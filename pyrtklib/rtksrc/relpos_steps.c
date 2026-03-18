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
        ctx->dt = intpres(obs[0].time, obs + nu, nr, ctx->nav, rtk,
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
