
save_figures <- function(name, path = NULL, scale){
    
    # save .jpg (for presentation)
    ggsave(
        paste0(name, ".jpg"), path = path,
        scale = scale, width = 1.61803, height = 1, units = "cm"
    )
    
    # save .eps (for publication)
    ggsave(
        paste0(name, ".eps"), path = path,
        device = "eps", dpi = "retina",
        scale = scale, width = 1.61803, height = 1, units = "cm"
    )
}


proposal_vs_target = function(
        censor_type, age_screen, theta, t0,
        include_legend = FALSE,
        path_figure = "slurm/output/figures/misc/", scale = 5
        ) {
    
    # prepare data
    if(censor_type == "censored") {
        screen_detected = rep(0, length(age_screen))
        censor_time = max(age_screen)+1.5
        n_screen_positive = sum(screen_detected)
        
        endpoints = compute_endpoints_i(age_screen, censor_type, censor_time, t0)
        prob_tau = compute_prob_tau_i(censor_type, censor_time, endpoints, theta, t0)
    }else if(censor_type == "screen") {
        screen_detected = c(rep(0, length(age_screen) - 1), 1)
        censor_time = max(age_screen)
        n_screen_positive = sum(screen_detected)
        
        endpoints = compute_endpoints_i(age_screen, censor_type, censor_time, t0)
        prob_tau = compute_prob_tau_i(censor_type, censor_time, endpoints, theta, t0)
    }

    labels = c('Proposal', expression(paste('Target (', I[i], ' = 0)')), expression(paste('Target (', I[i], ' = 1)')))
    
    # line plot
    stepsize = 0.01
    g = tibble(
        taus=seq(censor_time-17.5, censor_time-stepsize, by = stepsize),
        ll_indolent0=taus%>%map_dbl(dlog_likelihood_i, indolent_i=0, censor_type, censor_time, age_screen, n_screen_positive, theta, t0),
        ll_indolent1=taus%>%map_dbl(dlog_likelihood_i, indolent_i=1, censor_type, censor_time, age_screen, n_screen_positive, theta, t0),
        log_prop=taus%>%map_dbl(dlog_prop_tau_HP_i, censor_type, censor_time, endpoints, prob_tau, theta, t0),
        Target0 = exp(ll_indolent0) / sum(exp(ll_indolent0)) / stepsize,
        Target1 = exp(ll_indolent1) / sum(exp(ll_indolent1)) / stepsize,
        Proposal =exp(log_prop) / sum(exp(log_prop)) / stepsize
    ) %>% 
        pivot_longer(cols = Target0:Proposal, names_to = "Distribution", values_to = "Density") %>%
        ggplot(aes(taus, Density, linetype=Distribution)) + 
        geom_line(alpha=0.5)+
        labs(x = expression(t[i]^HP), y = "Density") +
        scale_linetype_discrete(labels=labels) +
        theme(legend.text.align = 0)
    
    if(!include_legend)  g = g + theme(legend.position="none")
    
    print(g)
    
    # name of figure file
    name = paste0("proposal_target - censor_type=", censor_type)
    # save figure
    save_figures(name, path_figure, scale)
}


