get_ESS_vs_power <- function(dataframe, gaussian = FALSE) {

  library(tidyverse)

  if (!gaussian) {

    extract_power_ess <- function(data) {

      data |>
        filter(N == 200) |>
        transmute(
          p_t_true = factor(
            p_t_true,
            levels = c(0.45, 0.60),
            labels = c("0.45", "0.60")
          ),

          HCP_power = HCP_power,
          HCP_ESS   = HCP_ESS,

          HW_power = HW_power,
          HW_ESS   = HW_ESS,

          GS_power = GS_power,
          GS_ESS   = GS_ESS,

          `SPRT(0.45)_power` = SPRT_power_0.45,
          `SPRT(0.45)_ESS`   = SPRT_ESS_0.45,

          `SPRT(0.6)_power` = SPRT_power_0.6,
          `SPRT(0.6)_ESS`   = SPRT_ESS_0.6,

          SPRT_adap_power = SPRT_power_adap,
          SPRT_adap_ESS   = SPRT_ESS_adap
        ) |>
        pivot_longer(
          cols = -p_t_true,
          names_to = c("method", ".value"),
          names_pattern = "(.+)_(power|ESS)"
        )
    }

    power_ess_N200 <- bind_rows(
      extract_power_ess(dataframe$res45),
      extract_power_ess(dataframe$res60)
    )

    shape_values <- c(
      "0.45" = 16,
      "0.60" = 17
    )

    method_labels <- power_ess_N200 |>
      group_by(method) |>
      summarise(
        ESS = mean(ESS),
        power = mean(power),
        .groups = "drop"
      )

    ggplot(
      power_ess_N200,
      aes(
        x = ESS,
        y = power,
        colour = method,
        shape = p_t_true
      )
    ) +
      geom_line(
        aes(group = method),
        linewidth = 0.6,
        alpha = 0.4
      ) +
      geom_point(size = 3.5) +
      geom_text(
        data = method_labels,
        aes(
          x = ESS,
          y = power,
          label = method,
          colour = method
        ),
        inherit.aes = FALSE,
        nudge_y = 0.012,
        size = 3.5,
        show.legend = FALSE
      ) +
      scale_colour_manual(
        values = c(
          "GS"         = "darkgreen",
          "HCP"        = "firebrick",
          "HW"         = "steelblue",
          "SPRT(0.45)" = "orange3",
          "SPRT(0.6)"  = "goldenrod",
          "SPRT(0.8)"  = "orange",
          "SPRT_adap"  = "#009E73"
        ),
        guide = "none"
      ) +
      scale_shape_manual(
        values = shape_values,
        name = expression("True " * p[t])
      ) +
      coord_cartesian(xlim = c(0, 200)) +
      labs(
        x = "Expected sample size",
        y = "Power"
      ) +
      theme_bw(base_size = 14) +
      theme(
        panel.grid.minor = element_blank()
      )


  } else {

    extract_power_ess <- function(data, p_true) {
      data |>
        filter(N == 200) |>
        transmute(
          p_t_true = factor(
            p_true,
            levels = c(0.6, 0.8),
            labels = c("0.6", "0.8")
          ),

          HCP_power = HCP_power,
          HCP_ESS   = HCP_ESS,

          HW_power = HW_power,
          HW_ESS   = HW_ESS,

          GS_power = GS_power,
          GS_ESS   = GS_ESS,

          UIE_power = UIE_power,
          UIE_ESS = UIE_ESS

        ) |>
        pivot_longer(
          cols = -p_t_true,
          names_to = c("method", ".value"),
          names_pattern = "(.+)_(power|ESS)"
        )
    }

    power_ess_N200 <- bind_rows(
      extract_power_ess(dataframe$res6, 0.6),
      extract_power_ess(dataframe$res8, 0.8)
    )

    shape_values <- c(
      "0.6" = 16,
      "0.8" = 17
    )

    method_labels <- power_ess_N200 |>
      group_by(method) |>
      summarise(
        ESS = mean(ESS),
        power = mean(power),
        .groups = "drop"
      )

    ggplot(
      power_ess_N200,
      aes(
        x = ESS,
        y = power,
        colour = method,
        shape = p_t_true
      )
    ) +
      geom_line(
        aes(group = method),
        linewidth = 0.6,
        alpha = 0.4
      ) +
      geom_point(size = 3.5) +
      geom_text(
        data = method_labels,
        aes(
          x = ESS,
          y = power,
          label = method,
          colour = method
        ),
        inherit.aes = FALSE,
        nudge_y = 0.012,
        size = 3.5,
        show.legend = FALSE
      ) +
      scale_colour_manual(
        values = c(
          "GS"         = "darkgreen",
          "HCP"        = "firebrick",
          "HW"         = "steelblue",
          "UIE"        = "goldenrod"
        ),
        guide = "none"
      ) +
      scale_shape_manual(
        values = shape_values,
        name = expression("True " * p[t])
      ) +
      coord_cartesian(xlim = c(0, 200)) +
      labs(
        x = "Expected sample size",
        y = "Power"
      ) +
      theme_bw(base_size = 14) +
      theme(
        panel.grid.minor = element_blank()
      )

  }


  }
