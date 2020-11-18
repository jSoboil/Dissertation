source("Part_3_Evidence_Synthesis_and_Probability_Array.R")

# ==========================================================================================
# Markov Model ------------------------------------------------------------
# ==========================================================================================
# The following code runs a Markov Model script for the HPV model developed by Sinanovic, E., 
# et al. 2009): "The potential cost-effectiveness of adding a human papillomavirus vaccine to 
# the cervical cancer screening programme in South Africa." 

# Note that this code sits on top of the source code for the evidence synthesis model as well as 
# the probability matrix (Part 3). 

# m_M_ad_1 is the Status Quo treatment model; m_M_ad_2 is the vaccine treatment model.

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

# Iterative solution of age-dependent cSTM model 1:
for (t in 1:n_t) {
 for(i in 1:n.sims) {
  # Fill in cohort trace
  m_M_ad_1[t + 1, , i] <- m_M_ad_1[t, , i] %*% a_P_1[ , , t, i]
  }
}
# Iterative solution of age-dependent cSTM model 2:
for (t in 1:n_t) {
 for(i in 1:n.sims) {
  # Fill in cohort trace
  m_M_ad_2[t + 1, , i] <- m_M_ad_2[t, , i] %*% a_P_2[ , , t, i]
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
ggplot(melt(apply(m_M_ad_1, c(1, 2), mean)), aes(x = Var1, y = value, 
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
v_S_ad_1 <- rowSums(apply(m_M_ad_1[, -30, ], c(1, 2), mean))  # vector with survival curve
ggplot(data.frame(Cycle = 0:n_t, Survival = v_S_ad_1), 
       aes(x = Cycle, y = Survival)) +
  geom_line(size = 1.3) +
  scale_y_continuous(labels = scales::percent) +
  xlab("Cycle") +
  ylab("Proportion alive") +
  theme_bw(base_size = 14) +
  theme()

# Visualisation of cohort proportions for Markov Model 2:
ggplot(melt(apply(m_M_ad_2, c(1, 2), mean)), aes(x = Var1, y = value, 
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
v_S_ad_2 <- rowSums(apply(m_M_ad_2[, -30, ], c(1, 2), mean))  # vector with survival curve
ggplot(data.frame(Cycle = 0:n_t, Survival = v_S_ad_2), 
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
v_prev_Infection_1 <- apply(m_M_ad_1[, "Infection", ], c(2, 1), mean)
v_prev_Infection_1 <- apply(v_prev_Infection_1, 2, mean) / v_S_ad_1
# Prevalence of Infection in Model 2:
v_prev_Infection_2 <- apply(m_M_ad_2[, "Infection", ], c(2, 1), mean)
v_prev_Infection_2 <- apply(v_prev_Infection_2, 2, mean) / v_S_ad_2
ggplot() +
 geom_line(aes(x = 0:n_t, y = v_prev_Infection_1), size = 1.3, colour = "skyblue", na.rm = TRUE) +
 geom_line(aes(x = 0:n_t, y = v_prev_Infection_2), size = 1.3, colour = "darkred", 
           alpha = 0.65, na.rm = TRUE) + 
 scale_y_continuous(labels = scales::percent) +
 xlab("Cycle") +
 ylab("Prevalence (%): HPV Infection") +
 xlim(15, 85) +
 theme_bw(base_size = 14) +
 theme()

# Mortality in Model 1:
v_D_ad_1 <- rowSums(apply(m_M_ad_1[, 30, ], c(1, 2), mean))
# Mortality in Model 1:
v_D_ad_2 <- rowSums(apply(m_M_ad_2[, 30, ], c(1, 2), mean))
ggplot() +
 geom_line(aes(x = 0:n_t, y = v_D_ad_1), size = 2, colour = "navyblue", na.rm = TRUE, 
            alpha = 1) +
 geom_point(aes(x = 0:n_t, y = v_D_ad_2), size = 1.5, colour = "red", 
           alpha = 0.95, na.rm = TRUE) + 
 # scale_y_continuous(labels = scales::percent) +
 xlab("Cycle") +
 ylab("Comparative no. Deaths") +
 xlim(15, 85) +
 theme_bw(base_size = 14) +
 theme()

# Visualisation of cohort proportions for Markov Model 1:
barplot(apply(m_M_ad_1, c(2, 1), mean), space = 1,
        ylab = "Proportion of patients in each state",
        main = "Status Quo Treatment", col = cols, cex.names = 1, cex.main = 1.5)
legend("topright", legend = c("Well", "Infection"), 
       fill = c("Well" = "#BD1B00",     "Infection" = "#FA7700"), 
       cex = .5, box.lwd = 2,
       horiz = FALSE)
# Visualisation of cohort proportions for Markov Model 2:
barplot(apply(m_M_ad_2, c(2, 1), mean), space = 1,
        ylab = "Proportion of patients in each state",
        main = "Vaccine Treatment", col = cols, cex.names = 1, cex.main = 1.5)
legend("topright", legend = c("Well", "Infection"), 
       fill = c("Well" = "#BD1B00",     "Infection" = "#FA7700"), 
       cex = .5, box.lwd = 2,
       horiz = FALSE)

# End file ----------------------------------------------------------------