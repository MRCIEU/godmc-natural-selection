library(tidyverse)
library(meta)
library(parallel)

load("tissue.RData")

blood <- clumped %>%
	select(cpg=cpg, snp=snp, id=id, beta_blood=Effect, se_blood=StdErr, cis=cis)

adipose <- a2 %>%
	select(id=MARKERNAME, beta_adipose=BETA2, se_adipose=SE)

brain <- b2 %>%
	select(id, beta_brain=BETA2, se_brain=SE)

dat <- left_join(blood, brain, by="id")
dat <- left_join(dat, adipose, by="id")

fn <- function(r)
{
	a <- metagen(c(r$beta_blood, r$beta_brain, r$beta_adipose), c(r$se_blood, r$se_brain, r$se_adipose))
	tibble(
		id = r$id,
		nstudies = a$k,
		Q = a$Q,
		Qpval = a$pval.Q,
		fe = a$TE.fixed,
		fe_se = a$seTE.fixed
	) %>% return()
}

out <- mclapply(1:nrow(dat), function(i) {
	message(i)
	fn(dat[i, ])
	}, mc.cores=16) %>%
	bind_rows() %>%
	inner_join(dat, ., by="id")


save(out, file="heterogeneity.rdata")

