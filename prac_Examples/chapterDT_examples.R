library(matlib)

# ================================================================================
# An incoherent, Dutch Book bet -------------------------------------------
# ================================================================================
# The following illustrates an incoherent series of bets:
p1 <- .3
p2 <- .1
p3 <- .4
# The set of probabilities do not meet the criterion P(A_1 \cap, ...,\cap A_n) = 1
sum(p1 + p2 + p3)

# In matrix form:
R <- matrix(c(1 - p1,   - p2,   - p3,
                - p1, 1 - p2,   - p3,
                - p1,   - p2, 1 - p3), 
            nrow = 3)
# Notice matrix determinant > 0:
det(R)
# and we can invert the matrix
inv(R)

# Here we give some numerical states to each event with probability pn:
s <- c(5, 2, 1)
# Accordingly, if we solve for the system of linear equations, the gains for each bet
# equal:
G <- solve(a = inv(R), b = s)
# Notice the sum have positive gains:
round(sum(G))

# Because this series of bets is inconsistent, we can solve for s again:
s_solved <- solve(a = R, b = G)
s_solved

# ================================================================================
# A coherent bet ----------------------------------------------------------
# ================================================================================
# The following illustrates a coherent series of bets:
p1 <- .3
p2 <- .3
p3 <- .4
# The set of probabilities do  meet the criterion P(A_1 \cap, ...,\cap A_n) = 1
sum(p1 + p2 + p3)

# In matrix form:
R <- matrix(c(1 - p1,   - p2,   - p3,
                - p1, 1 - p2,   - p3,
                - p1,   - p2, 1 - p3), 
            nrow = 3)
# Notice matrix determinant is 0:

det(R)
# and we cannot invert the matrix
inv(R)

# Here we give some numerical states to each event with probability pn:
s <- c(5, 2, 1)
# Accordingly, if we solve for the system of linear equations, the gains for each bet
# equal:
G <- R%*%s
# Notice the net 0 gains:
round(sum(G))

# End file ----------------------------------------------------------------