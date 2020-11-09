source("Evidence_Synthesis_and_Probability_Matrix.R")

# ==========================================================================================
# Markov Model ------------------------------------------------------------
# ==========================================================================================
# This is code that runs a Markov Model script for the HPV model develop by Sinanovic, E., et al. 
# 2009): "The potential cost-effectiveness of adding a human papillomavirus vaccine to the 
# cervical cancer screening programme in South Africa." Note that this code sits on top of the 
# code for the evidence synthesis model as well as the probability matrix. 

# m_M_ad_1 is the Status Quo treatment model; m_M_ad_2 is the vaccine treatment model.

# Initial state vector:
v_s_init <- c("Well" = 100000,     "Infection" = 0, "LSIL" = 0, "HSIL" = 0, 
                     "Stage-I Cancer" = 0,           "Stage-II Cancer" = 0, 
                   "Stage-III Cancer" = 0,           "Stage-IV Cancer" = 0, 
            "Detected.Stage-I Year 1" = 0,   "Detected.Stage-I Year 2" = 0, 
            "Detected.Stage-I Year 3" = 0,   "Detected.Stage-I Year 4" = 0,
            "Detected.Stage-I Year 5" = 0,  "Detected.Stage-II Year 1" = 0, 
           "Detected.Stage-II Year 2" = 0,  "Detected.Stage-II Year 3" = 0, 
           "Detected.Stage-II Year 4" = 0,  "Detected.Stage-II Year 5" = 0, 
          "Detected.Stage-III Year 1" = 0, "Detected.Stage-III Year 2" = 0, 
          "Detected.Stage-III Year 3" = 0, "Detected.Stage-III Year 4" = 0,
          "Detected.Stage-III Year 5" = 0,  "Detected.Stage-IV Year 1" = 0,
           "Detected.Stage-IV Year 2" = 0,  "Detected.Stage-IV Year 3" = 0,
           "Detected.Stage-IV Year 4" = 0,  "Detected.Stage-IV Year 5" = 0,
                    "Cancer Survivor" = 0,                     "Death" = 0)

# Initialize cohort trace for age-dependent cSTM
m_M_ad_1 <- array(matrix(0, nrow = n_t + 1, ncol = n_states),
                dim = c(c(n_t + 1, n_states), n.sims), 
                dimnames = list(0:n_t, v_n, 1:n.sims))
m_M_ad_2 <- array(matrix(0, nrow = n_t + 1, ncol = n_states),
                dim = c(c(n_t + 1, n_states), n.sims), 
                dimnames = list(0:n_t, v_n, 1:n.sims))

# Store the initial state vector in the first row of the cohort trace
m_M_ad_1[1, , ] <- v_s_init
m_M_ad_2[1, , ] <- v_s_init
# Initialize transition array for each cycle, for each n.sim:
a_A_1 <- array(0, dim = c(n_states, n_states, n_t + 1, n.sims), 
             dimnames = list(v_n, v_n, 0:n_t, 1:n.sims))
a_A_2 <- array(0, dim = c(n_states, n_states, n_t + 1, n.sims), 
             dimnames = list(v_n, v_n, 0:n_t, 1:n.sims))

# Set first slice of Array with the initial state vector in its diagonal
for (i in 1:n.sims) {
 diag(a_A_1[, , 1, i]) <- v_s_init
}
for (i in 1:n.sims) {
 diag(a_A_2[, , 1, i]) <- v_s_init
}

# Iterative solution of age-dependent cSTM model 1:
for (t in 1:n_t) {
 for(i in 1:n.sims) {
  # Fill in cohort trace
  m_M_ad_1[t + 1, , i] <- m_M_ad_1[t, , i] %*% a_P_1[ , , t, i]
  # Fill in transition dynamics array
  a_A_1[, , t + 1, i]  <- m_M_ad_1[t, , i] * a_P_1[ , , t, i]
  }
}
# Iterative solution of age-dependent cSTM model 2:
for (t in 1:n_t) {
 for(i in 1:n.sims) {
  # Fill in cohort trace
  m_M_ad_2[t + 1, , i] <- m_M_ad_2[t, , i] %*% a_P_2[ , , t, i]
  # Fill in transition dynamics array
  a_A_2[, , t + 1, i]  <- m_M_ad_2[t, , i] * a_P_2[ , , t, i]
  }
}

# Model run time ----------------------------------------------------------
end_time <- Sys.time()
end_time - start_time

# ==========================================================================================
# Visualisation of Markov Models ------------------------------------------
# ==========================================================================================
# Associates each state with a colour:
cols <- c("Well" = "#BD1B00",     "Infection" = "#FA7700", "LSIL" = "#FAD343",
          "HSIL" = "#FAF27F", "Stage-I Cancer" = "#9AFA65",
          "Stage-II Cancer" = "#89E05A", "Stage-III Cancer" = "#4B7A31",
          "Stage-IV Cancer" = "#243B17", "Detected.Stage-I Year 1" = "#00C8FA",
          "Detected.Stage-I Year 2" = "#00B4E0", "Detected.Stage-I Year 3" = "#009BC2",
          "Detected.Stage-I Year 4" = "#00627A", "Detected.Stage-I Year 5" = "#00313D",
          "Detected.Stage-II Year 1" = "#00FAD9",
          "Detected.Stage-II Year 2" = "#00E0C2",
          "Detected.Stage-II Year 3" = "#00C2A8",
          "Detected.Stage-II Year 4" = "#007A6A",
          "Detected.Stage-II Year 5" = "#003D35",
          "Detected.Stage-III Year 1" = "#9800FA",
          "Detected.Stage-III Year 2" = "#8700E0",
          "Detected.Stage-III Year 3" = "#7400C2",
          "Detected.Stage-III Year 4" = "#49007A",
          "Detected.Stage-III Year 5" = "#25003D",
          "Detected.Stage-IV Year 1" = "#E6FAFA",
          "Detected.Stage-IV Year 2" = "#CEE0E0",
          "Detected.Stage-IV Year 3" = "#B2C2C2",
          "Detected.Stage-IV Year 4" = "#717A7A",
          "Detected.Stage-IV Year 5" = "#4C5353",
          "Cancer Survivor" = "#3EFA80",
          "Death" = "#B3BABA")
# Associates similar states with similar line types:
lty <-  c("Well" = 1,     "Infection" = 2, "LSIL" = 3,
          "HSIL" = 3, "Stage-I Cancer" = 4,
          "Stage-II Cancer" = 4, "Stage-III Cancer" = 4,
          "Stage-IV Cancer" = 4, "Detected.Stage-I Year 1" = 5,
          "Detected.Stage-I Year 2" = 5, "Detected.Stage-I Year 3" = 5,
          "Detected.Stage-I Year 4" = 5, "Detected.Stage-I Year 5" = 5,
          "Detected.Stage-II Year 1" = 5,
          "Detected.Stage-II Year 2" = 5,
          "Detected.Stage-II Year 3" = 5,
          "Detected.Stage-II Year 4" = 5,
          "Detected.Stage-II Year 5" = 5,
          "Detected.Stage-III Year 1" = 5,
          "Detected.Stage-III Year 2" = 5,
          "Detected.Stage-III Year 3" = 5,
          "Detected.Stage-III Year 4" = 5,
          "Detected.Stage-III Year 5" = 5,
          "Detected.Stage-IV Year 1" = 5,
          "Detected.Stage-IV Year 2" = 5,
          "Detected.Stage-IV Year 3" = 5,
          "Detected.Stage-IV Year 4" = 5,
          "Detected.Stage-IV Year 5" = 5,
          "Cancer Survivor" = 6,
          "Death" = 1)

# Visualisation of cohort proportions for Markov Model 1:
ggplot(melt(apply(m_M_ad_1, c(1, 2), mean)), aes(x = Var1, y = value, 
                      color = Var2, linetype = Var2)) +
 geom_line(size = 1) +
 scale_colour_manual(name = "Health state", 
                     values = cols) +
  scale_linetype_manual(name = "Health state",
                       values = lty) +
  # scale_y_continuous(labels = scales::percent) +
  xlab("Cycle") +
  ylab("Proportion of the cohort") +
  theme_light(base_size = 14) +
  theme(legend.position = "bottom", 
        legend.background = element_rect(fill = NA),
        text = element_text(size = 15))
# Survival curve for Markov Model 1:
v_S_ad_1 <- rowSums(m_M_ad_1[, -30, ]) / n.sims # vector with survival curve
ggplot(data.frame(Cycle = 0:n_t, Survival = v_S_ad_1), 
       aes(x = Cycle, y = Survival)) +
  geom_line(size = 1.3) +
  xlab("Cycle") +
  ylab("Proportion alive") +
  theme_bw(base_size = 14) +
  theme()

# Visualisation of cohort proportions for Markov Model 2:
ggplot(melt(apply(m_M_ad_2, c(1, 2), mean)), aes(x = Var1, y = value, 
                      color = Var2, linetype = Var2)) +
 geom_line(size = 1) +
 scale_colour_manual(name = "Health state", 
                     values = cols) +
  scale_linetype_manual(name = "Health state",
                       values = lty) +
  # scale_y_continuous(labels = scales::percent) +
  xlab("Cycle") +
  ylab("Proportion of the cohort") +
  theme_light(base_size = 14) +
  theme(legend.position = "bottom", 
        legend.background = element_rect(fill = NA),
        text = element_text(size = 15))
# Survival curve for Markov Model 2:
v_S_ad_2 <- rowSums(m_M_ad_2[, -30, ]) / n.sims # vector with survival curve
v_S_ad_2

ggplot(data.frame(Cycle = 0:n_t, Survival = v_S_ad_2), 
       aes(x = Cycle, y = Survival)) +
  geom_line(size = 1.3) +
  xlab("Cycle") +
  ylab("Proportion alive") +
  theme_bw(base_size = 14) +
  theme()

# Life expectancy for average individual in Markov model 1 cohort:
le_ad_1 <- sum(v_S_ad_1) / 100000
le_ad_1
# Life expectancy for average individual in Markov model  cohort:
le_ad_2 <- sum(v_S_ad_2) / 100000
le_ad_2

# Prevalence of LSIL in Model 1:
v_prev_LSIL_1 <- apply(m_M_ad_1[, "Infection", ], 1, mean) / v_S_ad_1
# Prevalence of LSIL in Model 2:
v_prev_LSIL_2 <- apply(m_M_ad_2[, "Infection", ], 1, mean) / v_S_ad_2
ggplot() +
 geom_line(aes(x = 0:n_t, y = v_prev_LSIL_1), size = 1.3, colour = "skyblue") +
 geom_line(aes(x = 0:n_t, y = v_prev_LSIL_2), size = 1.3, colour = "darkred", 
           alpha = 0.65) + 
 scale_y_continuous(labels = scales::percent) +
 xlab("Cycle") +
 ylab("Prevalence (%): HPV Infection") +
 xlim(15, 85) +
 theme_bw(base_size = 14) +
 theme()
# Prevalence of Stage IV Cancer in Model 1:
v_prev_StageIV_1 <- apply(m_M_ad_1[, "Stage-IV Cancer", ], 1, mean) / v_S_ad_1
# Prevalence of Stage IV Cancer in Model 2:
v_prev_StageIV_2 <- apply(m_M_ad_2[, "Stage-IV Cancer", ], 1, mean) / v_S_ad_2
ggplot() +
 geom_line(aes(x = 0:n_t, y = v_prev_StageIV_1), size = 1.3, colour = "skyblue") +
 geom_line(aes(x = 0:n_t, y = v_prev_StageIV_2), size = 1.3, colour = "darkred", 
           alpha = 0.65) + 
 scale_y_continuous(labels = scales::percent) +
 xlab("Cycle") +
 ylab("Prevalence (%): Stage IV Cancer") +
 xlim(15, 85) +
 theme_bw(base_size = 14) +
 theme()

# End file ----------------------------------------------------------------