# ==========================================================================================
# Markov Model ------------------------------------------------------------
# ==========================================================================================
# Note that this code sits on top of the source code for the evidence synthesis model and
# probability matrix (03_evidence_synthesis_and_probability_array.R). 
source("R/03_evidence_synthesis_and_probability_array.R")

# The following code runs a Markov Model script for the HPV model developed by Sinanovic, E., 
# et al. 2009): "The potential cost-effectiveness of adding a human papillomavirus vaccine 
# to the cervical cancer screening programme in South Africa." 

# m_M_SoC is the Status Quo treatment model; m_M_NT is the vaccine treatment model.

# Initial state vector:
v_s_init <- c("Well" = 1,     "Infection" = 0, "LSIL" = 0, "HSIL" = 0, 
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

# Initialize cohort trace for age-dependent cSTMs:
m_M_SoC <- array(matrix(0, nrow = n_cycles + 1, ncol = n_states),
                dim = c(n_cycles + 1, n_states, n.sims), 
                dimnames = list(0:n_cycles, v_names_states, 1:n.sims))
m_M_NT <- m_M_SoC
# Store the initial state vector in the first row of the cohort trace
m_M_SoC[1, , ] <- v_s_init
m_M_NT[1, , ] <- v_s_init

# Create transition-dynamics array under SoC and New Treatment
a_A_SoC <- array(0,
             dim = c(n_states, n_states, (n_cycles + 1), n.sims),
             dimnames = list(v_names_states, v_names_states, 0:n_cycles))
# Store the initial state vector in the diagonal of the first slice of a_A
for (i in 1:n.sims) {
 diag(a_A_SoC[, , 1, i]) <- v_s_init
 }
# New Treatment
a_A_NT <- a_A_SoC

# Loop through the number of simulations
for (i in 1:n.sims) {
 # Loop through the number of cycles
 for(t in 1:n_cycles){ 
  # Iterative solution of age-dependent cSTM for SoC
   # estimate the state vector for the next cycle (t + 1)
  m_M_SoC[t + 1, , i] <- m_M_SoC[t, , i] %*% a_P_SoC[, , t, i]
   # estimate the transition dynamics at t + 1
  a_A_SoC[, , t + 1, i] <- diag(m_M_SoC[t, , i]) %*% a_P_SoC[, , t, i]
  # Iterative solution of age-dependent cSTM for New Treatment
   # estimate the state vector for the next cycle (t + 1)
  m_M_NT[t + 1, , i] <- m_M_NT[t, , i] %*% a_P_NT[, , t, i]
   # estimate the transition dynamics at t + 1
  a_A_NT[, , t + 1, i] <- diag(m_M_NT[t, , i]) %*% a_P_NT[, , t, i]
  }
}

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
          "Detected.Stage-IV Year 1" = "#FC56BD",
          "Detected.Stage-IV Year 2" = "#F051B3",
          "Detected.Stage-IV Year 3" = "#D649A1",
          "Detected.Stage-IV Year 4" = "#B03C83",
          "Detected.Stage-IV Year 5" = "#702654",
          "Cancer Survivor" = "#3EFA80",
          "Death" = "#B3BABA")
# Associates similar states with similar line types:
lty <-  c("Well" = 1,     "Infection" = 2, "LSIL" = 3,
          "HSIL" = 3, "Stage-I Cancer" = 4,
          "Stage-II Cancer" = 4, "Stage-III Cancer" = 4,
          "Stage-IV Cancer" = 4, "Detected.Stage-I Year 1" = 2,
          "Detected.Stage-I Year 2" = 2, "Detected.Stage-I Year 3" = 2,
          "Detected.Stage-I Year 4" = 2, "Detected.Stage-I Year 5" = 2,
          "Detected.Stage-II Year 1" = 3,
          "Detected.Stage-II Year 2" = 3,
          "Detected.Stage-II Year 3" = 3,
          "Detected.Stage-II Year 4" = 3,
          "Detected.Stage-II Year 5" = 3,
          "Detected.Stage-III Year 1" = 4,
          "Detected.Stage-III Year 2" = 4,
          "Detected.Stage-III Year 3" = 4,
          "Detected.Stage-III Year 4" = 4,
          "Detected.Stage-III Year 5" = 4,
          "Detected.Stage-IV Year 1" = 5,
          "Detected.Stage-IV Year 2" = 5,
          "Detected.Stage-IV Year 3" = 5,
          "Detected.Stage-IV Year 4" = 5,
          "Detected.Stage-IV Year 5" = 5,
          "Cancer Survivor" = 6,
          "Death" = 1)

# Visualisation of cohort proportions for Markov Model 1:
ggplot(melt(apply(m_M_SoC, c(1, 2), mean)), aes(x = Var1, y = value, 
                      color = Var2, linetype = Var2)) +
 geom_line(size = 0.5) +
 scale_colour_manual(name = "Health state", 
                     values = cols) +
  scale_linetype_manual(name = "Health state",
                       values = lty) +
  scale_y_continuous(labels = scales::percent) +
  xlab("Cycle") +
  ylab("Proportion of the cohort") +
  theme_light(base_size = 14) +
  theme(legend.position = "bottom", 
        legend.background = element_rect(fill = NA),
        text = element_text(size = 15))
# Survival curve for Markov Model 1:
v_S_ad_1 <- rowSums(apply(m_M_SoC[, -30, ], c(1, 2), mean))  # vector with survival curve
ggplot(data.frame(Cycle = 0:n_cycles, Survival = v_S_ad_1), 
       aes(x = Cycle, y = Survival)) +
  geom_line(size = 1.3) +
  scale_y_continuous(labels = scales::percent) +
  xlab("Cycle") +
  ylab("Proportion alive") +
  theme_bw(base_size = 14) +
  theme()

# Visualisation of cohort proportions for Markov Model 2:
ggplot(melt(apply(m_M_NT, c(1, 2), mean)), aes(x = Var1, y = value, 
                      color = Var2, linetype = Var2)) +
 geom_line(size = 0.5) +
 scale_colour_manual(name = "Health state", 
                     values = cols) +
  scale_linetype_manual(name = "Health state",
                       values = lty) +
  scale_y_continuous(labels = scales::percent) +
  xlab("Cycle") +
  ylab("Proportion of the cohort") +
  theme_light(base_size = 14) +
  theme(legend.position = "bottom", 
        legend.background = element_rect(fill = NA),
        text = element_text(size = 15))
# Survival curve for Markov Model 2:
v_S_ad_2 <- rowSums(apply(m_M_NT[, -30, ], c(1, 2), mean))  # vector with survival curve
ggplot(data.frame(Cycle = 0:n_cycles, Survival = v_S_ad_2), 
       aes(x = Cycle, y = Survival)) +
  geom_line(size = 1.3) +
  scale_y_continuous(labels = scales::percent) +
  xlab("Cycle") +
  ylab("Proportion alive") +
  theme_bw(base_size = 14) +
  theme()

# Life expectancy for average individual in Markov model 1 cohort:
le_ad_1 <- sum(v_S_ad_1)
le_ad_1
# Life expectancy for average individual in Markov model  cohort:
le_ad_2 <- sum(v_S_ad_2)
le_ad_2

# Prevalence of Infection in Model 1:
v_prev_Infection_1 <- apply(m_M_SoC[, "Infection", ], c(2, 1), mean)
v_prev_Infection_1 <- apply(v_prev_Infection_1, 2, mean) / v_S_ad_1
# Prevalence of Infection in Model 2:
v_prev_Infection_2 <- apply(m_M_NT[, "Infection", ], c(2, 1), mean)
v_prev_Infection_2 <- apply(v_prev_Infection_2, 2, mean) / v_S_ad_2
ggplot() +
 geom_line(aes(x = 0:n_cycles, y = v_prev_Infection_1), size = 1.3, colour = "skyblue", 
           na.rm = TRUE) +
 geom_line(aes(x = 0:n_cycles, y = v_prev_Infection_2), size = 1.3, colour = "darkred", 
           alpha = 0.65, na.rm = TRUE) + 
 scale_y_continuous(labels = scales::percent) +
 xlab("Cycle") +
 ylab("Prevalence (%): HPV Infection") +
 xlim(15, 85) +
 theme_bw(base_size = 14) +
 theme()

# Mortality in Model 1:
v_D_ad_1 <- rowSums(apply(m_M_SoC[ , "Death", ], c(1, 2), mean)) / n.sims
# Mortality in Model 1:
v_D_ad_2 <- rowSums(apply(m_M_NT[ , "Death", ], c(1, 2), mean)) / n.sims
ggplot() +
 geom_line(aes(x = 0:n_cycles, y = v_D_ad_1), size = 2, colour = "navyblue", 
           na.rm = TRUE, alpha = 1) +
 geom_point(aes(x = 0:n_cycles, y = v_D_ad_2), size = 1.5, colour = "red", 
           alpha = 0.95, na.rm = TRUE) + 
 scale_y_continuous(labels = scales::percent) +
 xlab("Cycle") +
 ylab("Probability of Death") +
 xlim(15, 85) +
 theme_bw(base_size = 14) +
 theme()

# Visualisation of cohort proportions for Markov Model 1:
barplot(apply(m_M_SoC, c(2, 1), mean), space = 1,
        ylab = "Proportion of patients in each state",
        main = "Status Quo Treatment", col = cols, cex.names = 1, cex.main = 1.5)
legend("bottomright", legend = c("Well", "Infection", "Death"),
       fill = c("Well" = "#BD1B00", "Infection" = "#FA7700", "Death" = "#B3BABA"),
       cex = .65, box.lwd = 1.85, x = 146.5, y = 1.02)
# Visualisation of cohort proportions for Markov Model 2:
barplot(apply(m_M_NT, c(2, 1), mean), space = 1,
        ylab = "Proportion of patients in each state",
        main = "Vaccine Treatment", col = cols, cex.names = 1, cex.main = 1.5)
legend("bottomright", legend = c("Well", "Infection", "Death"),
       fill = c("Well" = "#BD1B00", "Infection" = "#FA7700", "Death" = "#B3BABA"),
       cex = .65, box.lwd = 1.85, x = 146.5, y = 1.02)

# End file ----------------------------------------------------------------