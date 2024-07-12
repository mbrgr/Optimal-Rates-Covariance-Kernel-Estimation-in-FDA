##### Plot for covariate functions #####



bw_df = bw_df |> mutate(n = factor(n, levels = c(100, 200, 400)),
                p = factor(p, levels = c(10, 50, 100, 200)))


ggplot(bw_df) + geom_line(aes(x =h, y = sup.err, color = n, linetype = p)) + lims(y = c(0.02,0.61))
