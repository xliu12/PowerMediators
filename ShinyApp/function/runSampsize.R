runSampsize <- function(
        n = 100,
        steps = 20,
        TarPow = 0.8,
        max_n = 100,
        sig.adjust = c("no_adjust", "bonferroni", "modified_bon1", "modified_bon2"),
        mediation = c("IIE_M1", "IIE_M2"),
        effect = c("all", "IIE_M1(1,,1)", "IIE_M1(1,,0)", "IIE_M1(0,,1)",
                   "IIE_M1(0,,0)", "IIE_M2(1,1,)", "IIE_M2(1,0,)", "IIE_M2(0,1,)", "IIE_M2(0,0,)"),
        power = c("familywise", "per-test"),
        plot = TRUE,
        verbose = TRUE,
        ...
) {
    power <- match.arg(power)  # validate power argument
    nstart <- n                # store initial sample size
    power_df <- data.frame()   # to collect all iterations
    plot_obj <- NULL           # initialize plot
    
    repeat {
        if (verbose) message("Running n = ", n)
        
        df <- runPower(n = n, ...)
        
        # Filter for sig.adjust, mediation, and effect into cond_df
        cond_df <- df[df$sig.adjust %in% sig.adjust, ]
        if (effect[1] == "all") {
            cond_df <- cond_df[cond_df$mediation %in% mediation, ]
        } else {
            cond_df <- cond_df[cond_df$mediation %in% mediation & cond_df$effect %in% effect, ]
        }
        
        # Append to cumulative power_df
        power_df <- dplyr::bind_rows(power_df, cond_df)
        
        power_col <- if (power == "familywise") "power_FW" else "power_PT"
        
        if (nrow(cond_df) > 0 && all(cond_df[[power_col]] >= TarPow)) break
        
        if (n + steps > max_n) {
            n <- max_n
            if (verbose) message("Target power not achieved, the maximum sample size of ", max_n, " will be used.")
            
            df <- runPower(n = n, ...)
            cond_df <- df[df$sig.adjust %in% sig.adjust, ]
            if (effect[1] == "all") {
                cond_df <- cond_df[cond_df$mediation %in% mediation, ]
            } else {
                cond_df <- cond_df[cond_df$mediation %in% mediation & cond_df$effect %in% effect, ]
            }
            power_df <- dplyr::bind_rows(power_df, cond_df)
            break
        }
        
        n <- n + steps
    }
    
    # Construct result summary table
    if (power == "familywise") {
        power_target_table <- cond_df[, c("mediation", "sig.adjust", "n", "power_FW")]
    } else {
        power_target_table <- cond_df[, c("mediation", "effect", "sig.adjust", "n", "power_PT")]
    }
    power_target_table <- unique(power_target_table)
    
    # Create visualization if requested
    if (plot) {
        power_df$yval <- if (power == "familywise") power_df$power_FW else power_df$power_PT
        power_df$colval <- if (power == "familywise") power_df$sig.adjust else power_df$effect
        base_plot <- ggplot2::ggplot(
            power_df,
            aes(x = n, y = yval, color = colval, shape = sig.adjust)
        )
        
        plot_obj <- base_plot +
            ggplot2::geom_point(size = 2) +
            ggplot2::geom_line(linewidth = 0.5) +
            ggplot2::geom_hline(yintercept = TarPow, linetype = "dotted", color = "red") +
            ggplot2::annotate("text", x = nstart, y = TarPow,
                              label = paste0("Target Power = ", TarPow),
                              vjust = -0.5, hjust = 0,
                              color = "red", size = 3.5, fontface = "italic") +
            ggplot2::scale_x_continuous(breaks = seq(nstart, max_n, by = steps)) +
            ggplot2::scale_shape_manual(
                values = c(
                    "no_adjust" = 4,
                    "bonferroni" = 16,
                    "modified_bon1" = 17,
                    "modified_bon2" = 15
                ),
                labels = c(
                    "no_adjust" = "No Adjustment",
                    "bonferroni" = "Bonferroni",
                    "modified_bon1" = "Modified Bonferroni I",
                    "modified_bon2" = "Modified Bonferroni II"
                ),
                name = "Adjustment"
            ) +
            ggplot2::labs(
                title = if (power == "familywise") "Familywise Power Analysis" else "Per-test Power Analysis",
                x = "Sample Size",
                y = "Power",
                color = if (power == "familywise") "Adjustment" else "Effect",
                shape = "Adjustment"
            ) +
            facet_wrap(~ mediation) +
            ggplot2::theme_bw(base_size = 12) +
            ggplot2::theme(
                panel.grid.major = ggplot2::element_line(color = "grey85", size = 0.3),
                panel.grid.minor = ggplot2::element_blank(),
                axis.text = ggplot2::element_text(color = "black"),
                axis.title = ggplot2::element_text(face = "bold"),
                plot.title = ggplot2::element_text(face = "bold", size = 15, hjust = 0.5),
                legend.position = "right",
                legend.key = ggplot2::element_blank(),
                strip.background = ggplot2::element_blank(),
                strip.text = ggplot2::element_text(face = "bold")
            )
    }
    
    return(list(
        power_df = power_df,
        target_power_table = power_target_table,
        plot = plot_obj
    ))
}
