# =====================================================================
#  BugCast - interactive R Shiny dashboard (Review 2)
#  Run in RStudio:  shiny::runApp("path/to/bugcast")
#  or from R:       setwd("bugcast"); shiny::runApp()
#  Needs: install.packages("shiny")
# =====================================================================
library(shiny)
# ---- BugCast core functions (built in, no other file needed) ----
# =====================================================================
#  BugCast - core functions
#  BAMAT207 Probability & Statistics | Review 2
#  Goel-Okumoto (GO) and Delayed S-shaped (DSS) NHPP reliability models
#  Used by: bugcast_analysis.R (report + plots) and app.R (Shiny dashboard)
#  Base R only (no extra packages needed for the analysis).
# =====================================================================

# ---------------------------------------------------------------------
# 1. DATA
# ---------------------------------------------------------------------

# Musa SYS1 (DACS system 1): 136 inter-failure times in CPU seconds.
sys1_interfailure <- c(
  3,30,113,81,115,9,2,91,112,15,138,50,77,24,108,88,670,120,26,114,325,55,
  242,68,422,180,10,1146,600,15,36,4,0,8,227,65,176,58,457,300,97,263,452,
  255,197,193,6,79,816,1351,148,21,233,134,357,193,236,31,369,748,0,232,330,
  365,1222,543,10,16,529,379,44,129,810,290,300,529,281,160,828,1011,445,296,
  1755,1064,1783,860,983,707,33,868,724,2323,2930,1461,843,12,261,1800,865,
  1435,30,143,108,0,3110,1247,943,700,875,245,729,1897,447,386,446,122,990,
  948,1082,22,75,482,5509,100,10,1071,371,790,6150,3321,1045,648,5485,1160,
  1864,4116)

# Cumulative failure times in HOURS (88,682 s = 24.63 h of testing).
load_sys1 <- function() {
  t <- cumsum(sys1_interfailure) / 3600
  list(t = t, T = max(t), name = "Musa SYS1 (CPU hours)", unit = "hours")
}

# Simulated GO data - same generator and seed as Review 1 (continuity).
simulate_go <- function(a, b, T) {
  N <- rpois(1, a * (1 - exp(-b * T)))
  u <- runif(N)
  sort(-log(1 - u * (1 - exp(-b * T))) / b)
}
load_simulated <- function(seed = 207, a = 120, b = 0.01, T = 300) {
  set.seed(seed)
  list(t = simulate_go(a, b, T), T = T,
       name = "Simulated GO (a=120, b=0.01)", unit = "time units")
}

# Delayed S-shaped simulator (used to test model selection).
simulate_dss <- function(a, b, T) {
  N <- rpois(1, a * (1 - (1 + b * T) * exp(-b * T)))
  if (N == 0) return(numeric(0))
  F <- function(x) (1 - (1 + b * x) * exp(-b * x)) / (1 - (1 + b * T) * exp(-b * T))
  sort(sapply(runif(N), function(u) uniroot(function(x) F(x) - u, c(0, T))$root))
}

# Generic loader: CSV with one column of cumulative failure times
# (or created_at dates from GitHub issues -> converted to days since first).
load_csv_times <- function(path, end_time = NULL) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  x <- d[[1]]
  if (is.character(x)) {               # ISO dates from GitHub
    tt <- as.POSIXct(x, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    if (all(is.na(tt))) tt <- as.POSIXct(x, tz = "UTC")
    x <- as.numeric(difftime(tt, min(tt), units = "days"))
    unit <- "days"
  } else unit <- "time units"
  x <- sort(x[!is.na(x)])
  if (x[1] == 0) x[1] <- 1e-6         # avoid log(0) for the first failure
  T <- if (is.null(end_time)) max(x) else end_time
  list(t = x, T = T, name = basename(path), unit = unit)
}

# ---------------------------------------------------------------------
# 2. MODELS  (NHPP:  log L = sum(log lambda(t_i)) - m(T))
# ---------------------------------------------------------------------
models <- list(
  GO = list(
    label = "Goel-Okumoto",
    m      = function(t, a, b) a * (1 - exp(-b * t)),
    lambda = function(t, a, b) a * b * exp(-b * t),
    G      = function(T, b) 1 - exp(-b * T)                 # m(T)/a
  ),
  DSS = list(
    label = "Delayed S-shaped",
    m      = function(t, a, b) a * (1 - (1 + b * t) * exp(-b * t)),
    lambda = function(t, a, b) a * b^2 * t * exp(-b * t),
    G      = function(T, b) 1 - (1 + b * T) * exp(-b * T)
  )
)

# MLE via the profile likelihood: for fixed b, a_hat(b) = n / G(T, b).
# GO uses the closed-form score equation (exactly as in Review 1);
# DSS uses a 1-D optimisation of the profile log-likelihood.
fit_model <- function(t, T = max(t), model = "GO") {
  n <- length(t); S <- sum(t); M <- models[[model]]
  if (n < 3) stop("Need at least 3 failures to fit a model")
  if (model == "GO") {
    g <- function(b) n / b - S - n * T * exp(-b * T) / (-expm1(-b * T))
    lo <- 1e-6 / T; hi <- 1000 / T
    if (g(lo) <= 0) stop("No reliability growth: GO model not applicable")
    b <- uniroot(g, c(lo, hi), tol = 1e-12)$root
  } else {
    prof <- function(lb) {
      b <- exp(lb); a <- n / M$G(T, b)
      sum(log(M$lambda(t, a, b))) - a * M$G(T, b)
    }
    o <- optimize(prof, c(log(1e-4 / T), log(100 / T)), maximum = TRUE, tol = 1e-10)
    b <- exp(o$maximum)
  }
  a <- n / M$G(T, b)
  logL <- sum(log(M$lambda(t, a, b))) - M$m(T, a, b)
  list(model = model, a = a, b = b, n = n, T = T, logL = logL,
       AIC = 4 - 2 * logL, BIC = 2 * log(n) - 2 * logL)
}

# Derived quantities -------------------------------------------------
m_fit   <- function(t, f) models[[f$model]]$m(t, f$a, f$b)
lam_fit <- function(t, f) models[[f$model]]$lambda(t, f$a, f$b)
residual_faults <- function(f) f$a - m_fit(f$T, f)          # a - m(T)
reliability <- function(x, f) exp(-(m_fit(f$T + x, f) - m_fit(f$T, f)))
# Under GO/DSS, R(x|T) -> exp(-residual) > 0 as x -> Inf (finite faults), so
# the mean time to next failure is infinite. Report the median instead.
p_no_more_failures <- function(f) exp(-residual_faults(f))
median_next_failure <- function(f) {
  if (p_no_more_failures(f) >= 0.5) return(Inf)
  uniroot(function(x) reliability(x, f) - 0.5, c(0, f$T * 1e3))$root
}

# ---------------------------------------------------------------------
# 3. GOODNESS OF FIT
# ---------------------------------------------------------------------
# Given n failures in (0,T], NHPP failure times are i.i.d. with
# CDF F(t) = m(t)/m(T).  One-sample KS test on that CDF.
ks_check <- function(t, f) {
  cdf <- function(q) m_fit(q, f) / m_fit(f$T, f)
  suppressWarnings(ks.test(t, cdf))       # ties (0-s gaps in SYS1) -> warning only
}

# ---------------------------------------------------------------------
# 4. PARAMETRIC BOOTSTRAP  (closes Research Gap 1)
# ---------------------------------------------------------------------
# Re-simulate failure histories from the fitted NHPP on (0,T], refit,
# and collect a, b, residual faults, lambda(T) and R(x|T).
simulate_from_fit <- function(f) {
  if (f$model == "GO") simulate_go(f$a, f$b, f$T) else simulate_dss(f$a, f$b, f$T)
}
bootstrap_fit <- function(f, x = 1, B = 1000, seed = 2025) {
  set.seed(seed)
  out <- matrix(NA_real_, B, 5,
                dimnames = list(NULL, c("a", "b", "residual", "intensity", "R")))
  for (k in 1:B) {
    tb <- simulate_from_fit(f)
    fb <- tryCatch(fit_model(tb, f$T, f$model), error = function(e) NULL)
    if (is.null(fb) || !is.finite(fb$a) || fb$a > 100 * f$a) next
    out[k, ] <- c(fb$a, fb$b, residual_faults(fb), lam_fit(f$T, fb), reliability(x, fb))
  }
  out[complete.cases(out), , drop = FALSE]
}
ci_table <- function(f, boot, x = 1, level = 0.95) {
  p <- c((1 - level) / 2, 1 - (1 - level) / 2)
  est <- c(a = f$a, b = f$b, residual = residual_faults(f),
           intensity = lam_fit(f$T, f), R = reliability(x, f))
  q <- apply(boot, 2, quantile, probs = p)
  data.frame(parameter = c("a (total faults)", "b (detection rate)",
                           "Residual faults", "Failure intensity at T",
                           sprintf("R(x = %g | T)", x)),
             estimate = est, lower = q[1, ], upper = q[2, ],
             boot_se = apply(boot, 2, sd), row.names = NULL)
}

# ---------------------------------------------------------------------
# 5. HOLDOUT VALIDATION (fit on first 70 %, predict the remaining 30 %)
# ---------------------------------------------------------------------
holdout_validate <- function(t, T = max(t), model = "GO", frac = 0.7, B = 500) {
  n <- length(t); k <- floor(frac * n); Ttr <- t[k]
  f <- fit_model(t[1:k], Ttr, model)
  test_t <- t[(k + 1):n]
  pred <- m_fit(test_t, f)                         # predicted cumulative count
  actual <- (k + 1):n
  exp_new <- m_fit(T, f) - m_fit(Ttr, f)           # expected failures in test window
  obs_new <- n - k
  pi95 <- qpois(c(0.025, 0.975), exp_new)          # naive: Poisson only (params fixed)
  # Bootstrap predictive interval: parameter uncertainty + Poisson noise
  bt <- bootstrap_fit(f, B = B, seed = 99)
  M <- models[[model]]
  mu_b <- M$m(T, bt[, "a"], bt[, "b"]) - M$m(Ttr, bt[, "a"], bt[, "b"])
  set.seed(100); new_b <- rpois(length(mu_b), pmax(mu_b, 0))
  bpi <- quantile(new_b, c(0.025, 0.975), names = FALSE)
  list(fit = f, k = k, T_train = Ttr, test_t = test_t, pred = pred, actual = actual,
       MAE = mean(abs(pred - actual)), RMSE = sqrt(mean((pred - actual)^2)),
       expected_new = exp_new, observed_new = obs_new,
       PI_lower = pi95[1], PI_upper = pi95[2],
       inside_PI = obs_new >= pi95[1] && obs_new <= pi95[2],
       BPI_lower = bpi[1], BPI_upper = bpi[2],
       inside_BPI = obs_new >= bpi[1] && obs_new <= bpi[2], boot = bt,
       RE_end = (m_fit(T, f) - n) / n)             # relative error at end of testing
}

# ---------------------------------------------------------------------
# 6. RELEASE-READINESS RULE  (closes Research Gap 5)
# ---------------------------------------------------------------------
# SHIP  if  R_hat(x|T) >= R_target  AND  bootstrap lower bound >= R_floor
# Otherwise report how much more testing the point estimate needs.
extra_test_time <- function(f, x, R_target) {
  # solve R(x | T') = R_target for T' (closed form for GO, numeric for DSS)
  Rt <- function(Tp) exp(-(m_fit(Tp + x, f) - m_fit(Tp, f)))
  if (Rt(f$T) >= R_target) return(0)
  hi <- f$T * 2
  while (Rt(hi) < R_target && hi < f$T * 1e4) hi <- hi * 2
  if (Rt(hi) < R_target) return(Inf)
  uniroot(function(Tp) Rt(Tp) - R_target, c(f$T, hi))$root - f$T
}
release_decision <- function(f, boot, x = 1, R_target = 0.95, R_floor = 0.90, level = 0.95) {
  R_hat <- reliability(x, f)
  R_low <- unname(quantile(boot[, "R"], (1 - level) / 2))
  ship <- R_hat >= R_target && R_low >= R_floor
  list(R_hat = R_hat, R_low = R_low, ship = ship,
       verdict = if (ship) "SHIP" else "KEEP TESTING",
       extra_time = extra_test_time(f, x, R_target),
       x = x, R_target = R_target, R_floor = R_floor)
}

# Readiness over time: refit at checkpoints (as if testing stopped there)
readiness_trajectory <- function(t, model = "GO", x = 1, checkpoints = 10,
                                 B = 200, start_frac = 0.3) {
  n <- length(t)
  idx <- unique(round(seq(max(10, floor(start_frac * n)), n, length.out = checkpoints)))
  res <- lapply(idx, function(k) {
    f <- tryCatch(fit_model(t[1:k], t[k], model), error = function(e) NULL)
    if (is.null(f)) return(NULL)
    bt <- bootstrap_fit(f, x = x, B = B, seed = k)
    data.frame(k = k, T = t[k], R = reliability(x, f),
               R_low = quantile(bt[, "R"], 0.025), R_up = quantile(bt[, "R"], 0.975),
               residual = residual_faults(f))
  })
  do.call(rbind, res)
}

# ---- end of core functions ----

# ---- look & feel ---------------------------------------------------------
BG <- "#0B0B12"; FG <- "#DDDDE8"; GO_C <- "#FF4D4D"; DSS_C <- "#6C5CFF"
OK_C <- "#3DDC97"; WARN_C <- "#FFD166"; BAND <- adjustcolor("#6C5CFF", 0.30)
dark_par <- function(...) par(bg = BG, fg = FG, col.axis = FG, col.lab = FG, col.main = "white",
                              mar = c(4.5, 4.8, 3, 1.5), ...)
gridd <- function() grid(col = "#2A2A3A", lty = 1)
leg <- function(...) legend(..., bty = "n", text.col = FG)

css <- "
body{background:#06060B;color:#DDDDE8;font-family:'Segoe UI',Helvetica,Arial,sans-serif}
.well{background:#11111B;border:1px solid #23233A;color:#DDDDE8}
.nav-tabs>li>a{color:#9A9AC0}.nav-tabs>li.active>a,.nav-tabs>li.active>a:focus{background:#11111B;color:#fff;border-color:#23233A}
.nav-tabs{border-bottom:1px solid #23233A}
h2{font-family:Georgia,serif;letter-spacing:2px;color:#fff}
.kpi{background:#11111B;border:1px solid #23233A;border-radius:10px;padding:14px 16px;margin-bottom:14px}
.kpi .lab{font-size:12px;letter-spacing:1px;text-transform:uppercase;color:#9A9AC0}
.kpi .val{font-size:28px;font-weight:700;color:#fff}
.kpi .sub{font-size:12px;color:#9A9AC0}
.verdict{border-radius:12px;padding:18px;text-align:center;font-size:30px;font-weight:800;letter-spacing:3px}
.ship{background:rgba(61,220,151,.15);border:2px solid #3DDC97;color:#3DDC97}
.keep{background:rgba(255,77,77,.12);border:2px solid #FF4D4D;color:#FF4D4D}
table{color:#DDDDE8!important}.table>thead>tr>th,.table>tbody>tr>td{border-color:#23233A!important}
.form-control,.selectize-input{background:#0B0B12!important;color:#DDDDE8!important;border-color:#23233A!important}
.irs-bar,.irs-single{background:#6C5CFF!important;border-color:#6C5CFF!important}
.btn{background:#6C5CFF;color:#fff;border:none}
.note{color:#9A9AC0;font-size:13px}
"
kpi <- function(label, value, sub = "") div(class = "kpi", div(class = "lab", label),
                                             div(class = "val", value), div(class = "sub", sub))

# ---- UI --------------------------------------------------------------------
ui <- fluidPage(
  tags$head(tags$style(HTML(css))),
  titlePanel(div(h2("BUGCAST"), div(class = "note",
    "Software reliability growth & release-readiness dashboard  |  BAMAT207 Probability and Statistics"))),
  sidebarLayout(
    sidebarPanel(width = 3,
      selectInput("ds", "Failure data",
                  c("Musa SYS1 (real, 136 failures)" = "sys1",
                    "Simulated GO (Review 1)" = "simgo",
                    "Simulated S-shaped" = "sims",
                    "Upload CSV (times or GitHub dates)" = "csv")),
      conditionalPanel("input.ds == 'csv'",
        fileInput("file", "CSV: first column = cumulative failure times, or ISO dates", accept = ".csv")),
      radioButtons("model", "Model", c("Auto (lowest AIC)" = "auto", "Goel-Okumoto" = "GO",
                                       "Delayed S-shaped" = "DSS")),
      numericInput("x", "Mission time x (same units as data)", 1, min = 0.01, step = 0.5),
      sliderInput("Rt", "Target reliability R*", 0.5, 0.999, 0.95, step = 0.005),
      sliderInput("Rf", "Required lower 95% bound", 0.5, 0.999, 0.90, step = 0.005),
      sliderInput("frac", "Holdout: training fraction", 0.5, 0.9, 0.7, step = 0.05),
      selectInput("B", "Bootstrap replicates", c(200, 500, 1000, 2000), selected = 500),
      div(class = "note", "Change the data or model and every tab refits by maximum likelihood.")
    ),
    mainPanel(width = 9,
      tabsetPanel(
        tabPanel("Overview",
          br(), uiOutput("verdict"), br(), uiOutput("kpis"),
          plotOutput("p_fit_small", height = 330)),
        tabPanel("Model fit",
          br(), plotOutput("p_fit", height = 430), tableOutput("t_cmp"),
          div(class = "note", "AIC = 2k - 2 log L (lower is better). KS: H0 = failure times follow F(t) = m(t)/m(T).")),
        tabPanel("Failure intensity", br(), plotOutput("p_int", height = 460)),
        tabPanel("Goodness of fit", br(), plotOutput("p_ks", height = 420), verbatimTextOutput("ks_txt")),
        tabPanel("Bootstrap CIs",
          br(), plotOutput("p_boot", height = 380), tableOutput("t_ci"),
          downloadButton("dl_ci", "Download CI table")),
        tabPanel("Holdout validation",
          br(), plotOutput("p_ho", height = 440), tableOutput("t_ho")),
        tabPanel("Reliability & release",
          br(), plotOutput("p_rel", height = 360), plotOutput("p_traj", height = 360)),
        tabPanel("Data", br(), verbatimTextOutput("data_txt"), tableOutput("t_data"))
      )
    )
  )
)

# ---- server ------------------------------------------------------------------
server <- function(input, output, session) {

  dat <- reactive({
    switch(input$ds,
      sys1  = load_sys1(),
      simgo = load_simulated(),
      sims  = { set.seed(11); list(t = simulate_dss(120, 0.03, 300), T = 300,
                                   name = "Simulated S-shaped", unit = "time units") },
      csv   = { req(input$file); load_csv_times(input$file$datapath) })
  })
  observeEvent(input$ds, {        # sensible default mission time per dataset
    updateNumericInput(session, "x", value = if (input$ds == "sys1") 1 else 5)
  })

  fits <- reactive({
    d <- dat()
    lapply(c(GO = "GO", DSS = "DSS"), function(m)
      tryCatch(fit_model(d$t, d$T, m), error = function(e) NULL))
  })
  fit <- reactive({
    f <- fits()
    if (input$model != "auto") { validate(need(!is.null(f[[input$model]]),
                                    "This model cannot be fitted to these data (no reliability growth).")); return(f[[input$model]]) }
    ok <- Filter(Negate(is.null), f)
    validate(need(length(ok) > 0, "No model could be fitted."))
    ok[[which.min(sapply(ok, `[[`, "AIC"))]]
  })
  boot <- reactive({
    withProgress(message = "Bootstrapping...", value = 0.5,
                 bootstrap_fit(fit(), x = input$x, B = as.integer(input$B)))
  })
  ci <- reactive(ci_table(fit(), boot(), x = input$x))
  rd <- reactive(release_decision(fit(), boot(), x = input$x, R_target = input$Rt, R_floor = input$Rf))
  ho <- reactive(tryCatch(holdout_validate(dat()$t, dat()$T, fit()$model, frac = input$frac, B = 300),
                          error = function(e) NULL))

  # ---- Overview
  output$verdict <- renderUI({
    r <- rd()
    div(class = paste("verdict", if (r$ship) "ship" else "keep"), r$verdict,
        div(style = "font-size:14px;letter-spacing:0;font-weight:400;margin-top:6px",
            sprintf("R(%g | T) = %.3f, lower 95%% bound = %.3f. Rule: R >= %.3f and bound >= %.3f.%s",
                    r$x, r$R_hat, r$R_low, r$R_target, r$R_floor,
                    if (r$ship) "" else if (is.finite(r$extra_time))
                      sprintf(" About %.1f more %s of testing needed (point estimate).", r$extra_time, dat()$unit)
                    else " Target not reachable under this model.")))
  })
  output$kpis <- renderUI({
    f <- fit(); c <- ci()
    fluidRow(
      column(3, kpi("Model", models[[f$model]]$label, sprintf("AIC %.1f", f$AIC))),
      column(3, kpi("Total faults a", sprintf("%.1f", f$a), sprintf("95%% CI [%.1f, %.1f]", c$lower[1], c$upper[1]))),
      column(3, kpi("Residual faults", sprintf("%.1f", c$estimate[3]), sprintf("95%% CI [%.1f, %.1f]", c$lower[3], c$upper[3]))),
      column(3, kpi(sprintf("R(%g | T)", input$x), sprintf("%.3f", c$estimate[5]),
                    sprintf("95%% CI [%.3f, %.3f]", c$lower[5], c$upper[5]))))
  })

  fit_plot <- function() {
    d <- dat(); f <- fits(); tt <- seq(0, d$T, length.out = 400); dark_par()
    top <- max(length(d$t), sapply(Filter(Negate(is.null), f), `[[`, "a")) * 1.08
    plot(d$t, seq_along(d$t), type = "s", col = "white", lwd = 2, ylim = c(0, top),
         xlab = paste0("Testing time (", d$unit, ")"), ylab = "Cumulative failures", main = d$name); gridd()
    if (!is.null(f$GO))  lines(tt, m_fit(tt, f$GO), col = GO_C, lwd = 3)
    if (!is.null(f$DSS)) lines(tt, m_fit(tt, f$DSS), col = DSS_C, lwd = 3, lty = 2)
    abline(h = fit()$a, lty = 3, col = FG)
    leg("bottomright", legend = c("Observed", "Goel-Okumoto", "Delayed S-shaped", "Estimated a (selected)"),
        col = c("white", GO_C, DSS_C, FG), lwd = c(2, 3, 3, 1), lty = c(1, 1, 2, 3))
  }
  output$p_fit_small <- renderPlot(fit_plot(), bg = BG)
  output$p_fit <- renderPlot(fit_plot(), bg = BG)

  output$t_cmp <- renderTable({
    d <- dat(); do.call(rbind, lapply(Filter(Negate(is.null), fits()), function(f) data.frame(
      Model = models[[f$model]]$label, a = f$a, b = f$b, `log L` = f$logL, AIC = f$AIC, BIC = f$BIC,
      `KS p-value` = ks_check(d$t, f)$p.value, `Residual faults` = residual_faults(f), check.names = FALSE)))
  }, digits = 4)

  # ---- Intensity
  output$p_int <- renderPlot({
    d <- dat(); f <- fits(); tt <- seq(0, d$T, length.out = 400); dark_par()
    nb <- min(30, max(8, round(length(d$t) / 5)))
    h <- hist(d$t, breaks = seq(0, d$T, length.out = nb + 1), plot = FALSE)
    rate <- h$counts / diff(h$breaks)
    plot(h$mids, rate, type = "h", lwd = 10, col = adjustcolor(DSS_C, 0.55), ylim = c(0, max(rate) * 1.1),
         xlab = paste0("Testing time (", d$unit, ")"), ylab = "Failures per unit time",
         main = "Observed failure rate vs model intensity lambda(t)"); gridd()
    if (!is.null(f$GO))  lines(tt, lam_fit(tt, f$GO), col = GO_C, lwd = 3)
    if (!is.null(f$DSS)) lines(tt, lam_fit(tt, f$DSS), col = DSS_C, lwd = 3, lty = 2)
    points(d$T, lam_fit(d$T, fit()), pch = 19, cex = 1.6, col = WARN_C)
    text(d$T, lam_fit(d$T, fit()), sprintf("lambda(T) = %.3f", lam_fit(d$T, fit())), pos = 2, col = WARN_C)
    leg("topright", legend = c("Observed (binned)", "GO", "Delayed S"),
        col = c(adjustcolor(DSS_C, 0.55), GO_C, DSS_C), lwd = c(10, 3, 3), lty = c(1, 1, 2))
  }, bg = BG)

  # ---- KS
  output$p_ks <- renderPlot({
    d <- dat(); f <- fit(); tt <- seq(0, d$T, length.out = 400); dark_par()
    plot(ecdf(d$t), do.points = FALSE, verticals = TRUE, col = "white", lwd = 2,
         main = paste("Empirical vs", models[[f$model]]$label, "CDF"),
         xlab = paste0("Testing time (", d$unit, ")"), ylab = "F(t) = m(t)/m(T)"); gridd()
    lines(tt, m_fit(tt, f) / m_fit(d$T, f), col = GO_C, lwd = 3)
    leg("bottomright", legend = c("Empirical", "Model"), col = c("white", GO_C), lwd = c(2, 3))
  }, bg = BG)
  output$ks_txt <- renderPrint({
    k <- ks_check(dat()$t, fit()); print(k)
    cat(if (k$p.value > 0.05) "Decision: do not reject H0 at 5% - model is adequate.\n"
        else "Decision: reject H0 at 5% - model does not describe these data well.\n")
  })

  # ---- Bootstrap
  output$p_boot <- renderPlot({
    b <- boot(); c <- ci(); dark_par(mfrow = c(1, 2))
    for (j in c(3, 5)) {
      col <- c("a", "b", "residual", "intensity", "R")[j]
      xl <- if (j == 5) c(0, 1) else range(b[, col])
      hist(b[, col], breaks = 40, col = adjustcolor(DSS_C, 0.7), border = BG,
           main = c$parameter[j], xlab = "", xlim = xl)
      abline(v = c$estimate[j], col = GO_C, lwd = 3)
      abline(v = c(c$lower[j], c$upper[j]), col = WARN_C, lwd = 2, lty = 2)
      if (j == 5) abline(v = input$Rt, col = OK_C, lwd = 2, lty = 3)
    }
  }, bg = BG)
  output$t_ci <- renderTable({
    c <- ci(); names(c) <- c("Quantity", "MLE", "Lower 95%", "Upper 95%", "Bootstrap SE"); c
  }, digits = 4)
  output$dl_ci <- downloadHandler("bugcast_bootstrap_ci.csv", function(file) write.csv(ci(), file, row.names = FALSE))

  # ---- Holdout
  output$p_ho <- renderPlot({
    h <- ho(); validate(need(!is.null(h), "Model cannot be fitted on the training portion (no growth yet)."))
    d <- dat(); f <- h$fit; k <- h$k; M <- models[[f$model]]
    tp <- seq(h$T_train, d$T, length.out = 150)
    mu <- m_fit(tp, f) - m_fit(h$T_train, f)
    set.seed(5)
    bb <- sapply(tp, function(x) { m <- M$m(x, h$boot[, "a"], h$boot[, "b"]) - M$m(h$T_train, h$boot[, "a"], h$boot[, "b"])
                                   k + quantile(rpois(length(m), pmax(m, 0)), c(.025, .975)) })
    dark_par()
    plot(d$t, seq_along(d$t), type = "n", ylim = c(0, max(length(d$t), bb[2, ]) * 1.05),
         xlab = paste0("Testing time (", d$unit, ")"), ylab = "Cumulative failures",
         main = sprintf("Fit on first %d%%, predict the rest", round(100 * input$frac))); gridd()
    polygon(c(tp, rev(tp)), c(bb[1, ], rev(bb[2, ])), col = adjustcolor(OK_C, 0.18), border = NA)
    polygon(c(tp, rev(tp)), c(k + qpois(.025, mu), rev(k + qpois(.975, mu))), col = BAND, border = NA)
    lines(d$t[1:k], 1:k, type = "s", col = "white", lwd = 2)
    lines(d$t[(k + 1):length(d$t)], (k + 1):length(d$t), type = "s", col = WARN_C, lwd = 2)
    tf <- seq(0, d$T, length.out = 300); lines(tf, m_fit(tf, f), col = GO_C, lwd = 3)
    abline(v = h$T_train, lty = 2, col = FG)
    leg("bottomright", legend = c("Training", "Held-out", "Model fitted on training",
                                  "95% Poisson band", "95% bootstrap predictive band"),
        col = c("white", WARN_C, GO_C, BAND, adjustcolor(OK_C, 0.35)), lwd = c(2, 2, 3, 10, 10))
  }, bg = BG)
  output$t_ho <- renderTable({
    h <- ho(); req(h)
    data.frame(`Training failures` = h$k, `Predicted new` = h$expected_new, `Observed new` = h$observed_new,
               `Poisson 95% PI` = sprintf("[%d, %d]", h$PI_lower, h$PI_upper),
               `Bootstrap 95% PI` = sprintf("[%.0f, %.0f]", h$BPI_lower, h$BPI_upper),
               `Inside bootstrap PI` = ifelse(h$inside_BPI, "Yes", "No"),
               MAE = h$MAE, RMSE = h$RMSE, `Relative error at T` = h$RE_end, check.names = FALSE)
  }, digits = 3)

  # ---- Reliability & release
  output$p_rel <- renderPlot({
    f <- fit(); b <- boot(); M <- models[[f$model]]; dark_par()
    xs <- seq(0, input$x * 6, length.out = 150)
    q <- sapply(xs, function(x) quantile(exp(-(M$m(f$T + x, b[, "a"], b[, "b"]) - M$m(f$T, b[, "a"], b[, "b"]))),
                                         c(.025, .975)))
    plot(xs, reliability(xs, f), type = "n", ylim = c(0, 1), main = "R(x | T) after release",
         xlab = paste0("Mission time x (", dat()$unit, ")"), ylab = "Reliability"); gridd()
    polygon(c(xs, rev(xs)), c(q[1, ], rev(q[2, ])), col = BAND, border = NA)
    lines(xs, reliability(xs, f), col = GO_C, lwd = 3)
    abline(h = c(input$Rt, input$Rf), col = c(OK_C, WARN_C), lty = c(3, 2), lwd = 2)
    abline(v = input$x, col = FG, lty = 3)
  }, bg = BG)
  output$p_traj <- renderPlot({
    d <- dat()
    tr <- withProgress(message = "Refitting at checkpoints...", value = 0.5,
                       readiness_trajectory(d$t, fit()$model, x = input$x, checkpoints = 10, B = 150))
    validate(need(!is.null(tr) && nrow(tr) > 1, "Not enough checkpoints could be fitted."))
    dark_par()
    plot(tr$T, tr$R, type = "n", ylim = c(0, 1), xlim = c(0, d$T),
         xlab = paste0("Stop-testing time (", d$unit, ")"), ylab = sprintf("R(%g | T)", input$x),
         main = "Release readiness as testing progresses"); gridd()
    polygon(c(tr$T, rev(tr$T)), c(tr$R_low, rev(tr$R_up)), col = BAND, border = NA)
    lines(tr$T, tr$R, type = "b", pch = 19, col = GO_C, lwd = 3)
    lines(tr$T, tr$R_low, col = WARN_C, lty = 2, lwd = 2)
    abline(h = c(input$Rt, input$Rf), col = c(OK_C, WARN_C), lty = 3, lwd = 2)
    ok <- tr$R >= input$Rt & tr$R_low >= input$Rf
    if (any(ok)) abline(v = tr$T[which(ok)[1]], col = OK_C, lwd = 2)
  }, bg = BG)

  # ---- Data
  output$data_txt <- renderPrint({
    d <- dat(); n <- length(d$t)
    U <- (mean(d$t) - d$T / 2) / (d$T * sqrt(1 / (12 * n)))
    cat(sprintf("%s\nn = %d failures, T = %.3f %s\nLaplace trend U = %.2f (%s)\n", d$name, n, d$T, d$unit, U,
                if (U < -1.96) "significant reliability growth" else "no significant growth"))
  })
  output$t_data <- renderTable(data.frame(failure = seq_along(dat()$t), time = dat()$t)[1:min(25, length(dat()$t)), ],
                               digits = 3)
}

shinyApp(ui, server)
