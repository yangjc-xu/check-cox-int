rm(list = ls())
setwd("./check-cox-int")

# Figure 1
source("./simulation-code/prototype/cubic/main.R")
source("./simulation-code/prototype/exp/main.R")
source("./simulation-code/prototype/hr/main.R")
source("./simulation-code/prototype/hr_div/main.R")
source("./simulation-code/prototype/ind/main.R")
source("./simulation-code/prototype/log/main.R")
source("./simulation-code/prototype/quad/main.R")
source("./simulation-code/prototype/sqrt/main.R")

# Figure 2
source("./simulation-code/ff_n200_quad/main.R")
source("./simulation-code/ff_n400_quad/main.R")
source("./simulation-code/ff_n800_quad/main.R")
source("./simulation-code/ff_n200_twisted/main.R")
source("./simulation-code/ff_n400_twisted/main.R")
source("./simulation-code/ff_n800_twisted/main.R")
source("./simulation-code/omnibus_n200_quad/main.R")
source("./simulation-code/omnibus_n400_quad/main.R")
source("./simulation-code/omnibus_n800_quad/main.R")
source("./simulation-code/omnibus_n200_twisted/main.R")
source("./simulation-code/omnibus_n400_twisted/main.R")
source("./simulation-code/omnibus_n800_twisted/main.R")

# Figure 3
source("./simulation-code/prop_n200_mono/main.R")
source("./simulation-code/prop_n400_mono/main.R")
source("./simulation-code/prop_n800_mono/main.R")
source("./simulation-code/prop_n200_quad/main.R")
source("./simulation-code/prop_n400_quad/main.R")
source("./simulation-code/prop_n800_quad/main.R")



