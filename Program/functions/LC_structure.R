# Helper: add “missing” transitions with value = 0
# Some land-use transitions are not observed in the data but must be explicitly 
# included (with value = 0) to ensure a complete transition matrix for the model.
add_zero_transitions <- function(df, ns_vec, from_set, to_set, Ts, value = 0) {
  df %>%
    dplyr::bind_rows(
      tidyr::expand_grid(
        ns      = ns_vec,
        lu.from = from_set,
        lu.to   = to_set,
        Ts      = Ts
      ) %>%
        dplyr::mutate(value = value)
    )
}

# Main function
make_country_luc <- function(LandCoverChange,
                             map_HILDA_LUC,
                             Ts = 2015,
                             drop_from = c("ocean", "water", "not relevant"),
                             drop_to   = c("ocean", "water", "not relevant"),
                             other_lu  = c("cropland", "newforest", "otherland", "pasture", "urban")) {
  
  ns_vec <- unique(LandCoverChange$id_c)
  
  out <- LandCoverChange %>%
    tidyr::pivot_longer(-id_c, values_to = "AreaPerCode") %>%
    dplyr::mutate(name = stringr::str_remove(name, "X")) %>%
    dplyr::left_join(
      map_HILDA_LUC %>% dplyr::mutate(code = as.character(code)),
      by = c("name" = "code")
    ) %>%
    dplyr::filter(!(from %in% drop_from)) %>%
    dplyr::filter(!(to %in% drop_to)) %>%
    dplyr::mutate(AreaPerCode = as.numeric(AreaPerCode)) %>%
    dplyr::group_by(id_c) %>%
    dplyr::mutate(TotalArea = sum(AreaPerCode, na.rm = TRUE)) %>%
    dplyr::mutate(value = AreaPerCode / TotalArea) %>%
    dplyr::group_by(id_c, from, to) %>%
    dplyr::summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::rename(ns = id_c, lu.from = from, lu.to = to) %>%
    dplyr::mutate(Ts = Ts) %>%
    dplyr::select(ns, lu.from, lu.to, Ts, value) %>%
    dplyr::mutate(value = dplyr::if_else(is.na(value), 0, value)) %>%
    # FABLE: new forest definition
    dplyr::mutate(lu.to = dplyr::if_else(lu.from != "forest" & lu.to == "forest", "newforest", lu.to))
  
  # Add impossible/unused transitions with zeros
  out <- add_zero_transitions(out, ns_vec, other_lu, "forest", Ts, value = 0)
  out <- add_zero_transitions(out, ns_vec, "newforest", other_lu, Ts, value = 0)
  out <- add_zero_transitions(out, ns_vec, "forest", "newforest", Ts, value = 0)
  
  out
}
