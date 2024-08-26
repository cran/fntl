timing1_ex = function(R, n_levels)
{
	beta0_true = 0
	beta1_true = 1

	out = numeric(length(n_levels))
	for (k in seq_along(n_levels)) {
		n = n_levels[k]
		x = rnorm(n)
		pr_true = plogis(beta0_true + beta1_true * x)
		beta_hats = matrix(NA, R, 2)
		beta_trues = matrix(NA, R, 2)

		for (r in seq_len(R)) {
			# Generate a sample
			y = rbinom(n, size = 1, prob = pr_true)

			# A slow loglikelihood function
			loglik = function(par) {
				out = 0
				for (i in 1:n) {
					pr_i = plogis(par[1] + par[2] * x[i])
					out = out + y[i] * log(pr_i) + (1 - y[i]) * log1p(-pr_i)
				}
				return(out)
			}

			# Compute the MLE using L-BFGS-B
			init = c(0, 0)
			ctrl = list(fnscale = -1)
			optim_out = optim(init, loglik, method = "L-BFGS-B", control = ctrl)

			# Save estimates
			beta_hats[r,] = optim_out$par
			beta_trues[r,] = c(beta0_true, beta1_true)
		}

		out[k] = mean((beta_hats - beta_trues)^2)
	}

	return(out)
}
