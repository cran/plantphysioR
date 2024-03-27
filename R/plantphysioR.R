#' Chlorophyll'a' Concentration by Arnon method
#' @description Calculates Chlorophyll a Concentration according to Arnon(1949) method
#' @param A663 Absorbance at 663nm
#' @param A645 Absorbance at 645nm
#' @param v Final volume of solvent used in ml
#' @param w Fresh weight of the sample used in grams
#' @return Chlorophyll a in mg/g of fresh weight
#' @export
#' @references
#' Arnon, D. I. (1949). Copper enzymes in isolated chloroplasts. Polyphenoloxidase in Beta vulgaris. Plant physiology, 24(1), 1. \doi{10.1104/pp.24.1.1}
#' @examples chl_a(0.025, 0.041, 15, 1)
chl_a <- function(A663, A645, v, w) {
  ((12.7 * A663) - (2.69 * A645)) * v / (1000 * w)
} # mg/g fw

#' Chlorophyll b concentration
#' @description Calculates Chlorophyll b Concentration according to Arnon(1949) method
#' @param A645 Absorbance at 645nm
#' @param A663 Absorbance at 663nm
#' @param v Final volume of solvent used in ml
#' @param w Fresh weight of the sample used in grams
#' @return Chlorophyll b in mg/g of fresh weight
#' @export
#' @references Arnon, D. I. (1949). Copper enzymes in isolated chloroplasts. Polyphenoloxidase in Beta vulgaris. Plant physiology, 24(1), 1. \doi{10.1104/pp.24.1.1}
#' @examples chl_b(0.041, 0.025, 15, 1)
chl_b <- function(A645, A663, v, w) {
  ((22.9 * A645) - (4.68 * A663)) * v / (1000 * w)
} # mg/g fw

#' Total chlorophyll (a+b) concentration
#' @description  Calculate Total chlorophyll (a+b) concentration using method by Arnon (1949)
#' @param A645 Absorbance at 645nm
#' @param A663 Absorbance at 663nm
#' @param v Final volume of solvent used in ml
#' @param w Fresh weight of the sample used in grams
#' @return Total chlorophyll (a+b) in mg/g of fresh weight
#' @export
#' @references  Arnon, D. I. (1949). Copper enzymes in isolated chloroplasts. Polyphenoloxidase in Beta vulgaris. Plant physiology, 24(1), 1. \doi{10.1104/pp.24.1.1}
#' @examples chl_total(0.041, 0.025, 15, 1)
chl_total <- function(A645, A663, v, w) {
  ((20.2 * A645) + (8.02 * A663)) * v / (1000 * w)
} # mg/g fw


#' Total carotenoids content
#' @description  Calculate total carotenoids using Method by Lichtenthaler (1987)
#' @param A470 Absorbance at 470nm
#' @param A663 Absorbance at 663nm
#' @param A646 Absorbance at 646nm
#' @param fresh_weight Fresh weight of the sample used in grams
#' @return Carotenoids concentration in µg/ml
#' @export
#' @references Lichtenthaler, H. K. (1987). Chlorophylls and carotenoids: pigments of photosynthetic biomembranes. In Methods in enzymology (Vol. 148, pp. 350-382). Academic Press.
#' @examples caro_total(0.7, 0.041, 0.025, 1)
caro_total <- function(A470, A663, A646, fresh_weight) {
  carotenoid_concentration_ug_ml <- (1000 * A470 - 1.82 * (A663 + A646) + 85.02) / 198
}

#' Water potential of Polyethylene glycol (PEG) 6000
#' @description Calculate the corresponding water potential of PEG6000 when dissolved in 1l of water
#' @param peg Amount PEG600O in grams
#' @param C Temperature of the solution in degree centigrade
#' @return Water potential in bars
#' @export
#' @references Michel, B. E., & Kaufmann, M. R. (1973). The osmotic potential of polyethylene glycol 6000. Plant physiology, 51(5), 914-916.
#' @examples peg_6000(20, 25)
peg_6000 <- function(peg, C) {
  -(1.18 * 10^-2) * peg - (1.18 * 10^-4) * (peg^2) + (2.67 * 10^-4) * (peg * C) + (8.39 * 10^-7) * (peg^2) * (C)
}

#' Calculate Polyethylene glycol (PEG) 6000 requirement
#' @description Calculate Amount of PEG6000 required to reach desired water potential at given temperature
#' @param C Temperature of solution in degree centigrade
#' @param bar Water potential in bars
#'
#' @return PEG6000 required
#' @export
#' @references Michel, B. E., & Kaufmann, M. R. (1973). The osmotic potential of polyethylene glycol 6000. Plant physiology, 51(5), 914-916.
#' @examples calculate_PEG_6000(25, -4)
calculate_PEG_6000 <- function(C, bar) {
  PEG <- ((0.0118 - 0.000267 * C) - sqrt((0.000267 * C - 0.0118)^2 + 4 * bar * (0.000000839 * C - 0.000118))) / (2 * (0.000000839 * C - 0.000118))
  return(list("gram/liter of water"= PEG, "gram/cc of water" = PEG / 1000))
}

#' Calculate Stress tolerance index (STI) suggested by Fernandez (1992)
#' @description The genotype with high STI values will be tolerant to drought
#' @param Yp Yield under control condition
#' @param Ys Yield under stress condition
#'
#' @return STI
#' @export
#' @references Fernandez, G. C. (1993). Effective selection criteria for assessing plant stress tolerance.
#' @examples st_index(500, 350)
st_index <- function(Yp, Ys) {
  (Ys * Yp) / (Yp)^2
} # Fernandez,1992

#' Mean productivity Index (MPI)- by Hossain et al., (1990)
#' @description The genotype with high values of this index will be more desirable
#' @param Yp Yield under control condition
#' @param Ys Yield under stress condition
#'
#' @return Mean productivity Index
#' @references Hossain, A. B. S., Sears, R. G., Cox, T. S., & Paulsen, G. M. (1990). Desiccation tolerance and its relationship to assimilate partitioning in winter wheat. Crop Science, 30(3), 622-627.
mp_index <- function(Yp, Ys) {
  ((Ys + Yp) / 2)
}


#' Geometric mean productivity (GMP) by Fernandez (1992)
#' @description The genotype with high values of this index will be more desirable
#' @param Yp Yield under control condition
#' @param Ys Yield under stress condition
#'
#' @return GMP
#' @export
#' @references Fernandez, G. C. (1993). Effective selection criteria for assessing plant stress tolerance.
#' @examples
#' gmp(5, 3)
gmp <- function(Yp, Ys) {
  sqrt(Ys * Yp)
} # Fernandez,1992


#' Tolerance index -TOL by Hossain et al., (1990)
#' @description Higher the TOL value indicates the genotype is tolerant to stress
#' @param Yp Yield under control condition
#' @param Ys Yield under stress condition
#'
#' @return TOL
#' @export
#' @references Hossain, A. B. S., Sears, R. G., Cox, T. S., & Paulsen, G. M. (1990). Desiccation tolerance and its relationship to assimilate partitioning in winter wheat. Crop Science, 30(3), 622-627.
#' @examples tol_index(500, 350)
tol_index <- function(Yp, Ys) {
  (Yp - Ys)
} # Hossain etal., 1990


#' Stress susceptibility index (SSI) by Fischer and Maurer (1978)
#' @description The genotype with high SSI < 1 are more resistant to drought stress conditions
#' @param Yp Yield under control condition
#' @param Ys Yield under Stress condition
#' @param Ms Mean yield of all the genotypes under control condition
#' @param Mp Mean yield of all the genotypes under stress condition
#'
#' @return SSI
#' @export
#' @references Fischer, R. A., & Maurer, R. (1978). Drought resistance in spring wheat cultivars. I. Grain yield responses. Australian Journal of Agricultural Research, 29(5), 897-912.
#' @examples ss_index(500, 350, 450, 370)
ss_index <- function(Yp, Ys, Ms, Mp) {
  ((1 - Ys) / Yp) / ((1 - Ms) / Mp)} # Fischer and Maurer,1978

#' Yield reduction index or Yield Stability Index (YSI)
#' @description Higer YSI value depicts that particular genotype is stable under both normal and stressed conditions
#' @param Yp Yield under control condition
#' @param Ys Yield under stress condition
#'
#' @return YSI
#' @export
#' @references Bouslama, M., & Schapaugh Jr, W. T. (1984). Stress tolerance in soybeans. I. Evaluation of three screening techniques for heat and drought tolerance 1. Crop science, 24(5), 933-937.
#' @examples YSI(500, 350)
YSI <- function(Yp, Ys) {
  Ys / Yp
} # Bouslama and Schapaugh,1984

#' Yield reduction ratio (YR)
#' @description lesser the YR value more stable under stress conditions
#' @param Yp Yield under control condition
#' @param Ys Yield under stress condition
#'
#' @return YR
#' @export
#'
#' @examples YR_ratio(500, 350)
YR_ratio <- function(Yp, Ys) {
  1 - (Ys / Yp)
}

#' Drought resistant index (DRI)
#' @description The genotype with high values of this index will be more suitable for drought stress condition
#' @param Yp Yield under control condition
#' @param Ys Yield under stress condition
#'
#' @return DRI
#' @export
#'
#' @examples DRI(500, 350)
DRI <- function(Yp, Ys) {
  (Ys * (Ys / Yp) / Yp)
}
#' Harmonic Mean
#'
#' @param Yp Yield under control condition
#' @param Ys Yield under stress condition
#'
#' @return Harmonic mean
#' @export
#' @references Schneider, K. A., Rosales‐Serna, R., Ibarra‐Perez, F., Cazares‐Enriquez, B., Acosta‐Gallegos, J. A., Ramirez‐Vallejo, P., ... & Kelly, J. D. (1997). Improving common bean performance under drought stress. Crop science, 37(1), 43-50.
#' @examples HAM(500, 350)
HAM <- function(Yp, Ys) {
  2 * (Yp * Ys) / (Yp + Ys)}
#' Yield index (YI)
#'
#' @param Ys Yield under stress condition
#' @param Ms Mean Yield of all the genotypes under stress Condition
#'
#' @return Yield Index
#' @export
#' @references Gavuzzi, P., Rizza, F., Palumbo, M., Campanile, R. G., Ricciardi, G. L., & Borghi, B. (1997). Evaluation of field and laboratory predictors of drought and heat tolerance in winter cereals. Canadian journal of plant science, 77(4), 523-531.
#' @examples Y_index(500, 300)
Y_index <- function(Ys, Ms) {Ys / Ms}
#' Yield Reduction
#' @description  Claculate percent yield reduction over control
#' @param Yp Yield under control condition
#' @param Ys Yield under stress condition
#' @return YR
#' @export
#' @examples yield_reduction(500, 350)
yield_reduction <- function(Yp, Ys) {
  ((Yp - Ys) / Ys) * 100}
#' Relative Drought Index (RDI)
#' @description Calculates relative drought index according to Fisher and Wood (1979)
#' @param Yp Yield under control condition
#' @param Ys Yield under stress condition
#' @param Mp Mean Yield of all the genotypes under control Condition
#' @param Ms Mean Yield of all the genotypes under stress Condition
#'
#' @return RDI
#' @export
#' @references Fischer RA, Wood JT (1979) Drought resistance in spring wheat cultivars III. Yield association with morphological traits. Aust J Agr Res. 30: 1001-1020
#' @examples R_drought_index(500, 350, 400, 300)
R_drought_index <- function(Yp, Ys, Mp, Ms) {
  (Ys / Ys) / (Mp / Ms)
}

#' Golden Mean (GM)
#' @description  Calculates Golden mean value using Moradi et al.,(2012)
#' @param Yp Yield under control condition
#' @param Ys Yield under stress condition
#' @return GM
#' @export
#' @references Moradi H, Akbari GA, Khorasani SK, Ramshini HA (2012) Evaluation of drought tolerance in corn (Zea Mays L.) new hybrids with using stress tolerance indices. Eur J Sustain Dev 1. (3): 543-560
#' @examples
#' Golden_mean(500, 350)
Golden_mean <- function(Yp, Ys) {
  (Yp + Ys) / (Yp - Ys)
}
#' Abiotic Tolerance Index (ATI)
#' @description  Calculate abiotic tolerance index according to Moosavi et al. (2008)
#' @param Yp Yield under control condition
#' @param Ys Yield under Stress condition
#' @param Ms Mean yield of all the genotypes under control condition
#' @param Mp Mean yield of all the genotypes under stress condition
#' @return ATI
#' @export
#' @references Moosavi SS, Samadi YB, Naghavi MR, Zali AA, Dashti H, Pourshahbazi A (2008) Introduction of new indices to identify relative drought tolerance and resistance in wheat genotypes. Desert. 12: 165-178.
#' @examples ATI(500, 350, 400, 300)
ATI <- function(Yp, Ys, Mp, Ms) {
  ((Yp - Ys) / (Mp / Ms)) * sqrt(Yp * Ys)}
#' All indices combined
#' @description Function to all the indices related to biomass/ yield under different growth conditions
#' @param Yp Yield under control condition
#' @param Ys Yield under stress condition
#' @param Mp Mean yield of all the genotypes under control condition
#' @param Ms Mean yiels of all the genotyps under Stress condition
#' @return Indices Combined
#' @export
#' @examples Mp <- mean(yield_data$Yp)
#' Ms <- mean(yield_data$Ys)
#' Yp <- yield_data$Yp
#' Ys <- yield_data$Ys
#' all_indices(Yp, Ys, Mp, Ms)
all_indices <- function(Yp, Ys, Mp, Ms) {
  # Sub-functions
  st_index <- function(Yp, Ys) {
    (Ys * Yp) / (Yp)^2
  }
  mp_index <- function(Yp, Ys) {
    (Ys + Yp) / 2
  }

  gmp <- function(Yp, Ys) {
    sqrt(Ys * Yp)
  }
  tol_index <- function(Yp, Ys) {
    (Yp - Ys)
  }

  ss_index <- function(Yp, Ys, Mp, Ms) {
    ((1 - Ys) / Yp) / ((1 - Ms) / Mp)
  }

  YSI <- function(Yp, Ys) {
    Ys / Yp
  }

  YR_ratio <- function(Yp, Ys) {
    1 - (Ys / Yp)
  }

  DRI_index <- function(Yp, Ys) {
    Ys(Ys / Yp) / Yp
  }

  HAM <- function(Yp, Ys) {
    2 * (Yp * Ys) / (Yp + Ys)
  }

  Y_index <- function(Ys, Ms) {
    Ys / Ms
  }
  yield_reduction <- function(Yp, Ys) {
    ((Yp - Ys) / Ys) * 100
  }
  R_drought_index <- function(Yp, Ys, Mp, Ms) {
    (Ys / Ys) / (Mp / Ms)
  }
  Golden_mean <- function(Yp, Ys, Mp, Ms) {
    (Yp + Ys) / (Yp - Ys)
  }
  ATI <- function(Yp, Ys, Mp, Ms) {
    ((Yp - Ys) / (Mp / Ms)) * sqrt(Yp * Ys)
  }
  # Call sub-functions and return their results
  results <- list(
    StressToleranceIndex = st_index(Yp, Ys),
    MeanProductivityIndex = mp_index(Yp, Ys),
    GeometricMeanProductivity = gmp(Yp, Ys),
    Tolerance_Index = tol_index(Yp, Ys),
    StressSusceptibilityIndex = ss_index(Yp, Ys, Mp, Ms),
    YSI = YSI(Yp, Ys),
    YR_Ratio = YR_ratio(Yp, Ys),
    DRI = DRI(Yp, Ys),
    HAM = HAM(Yp, Ys),
    Y_Index = Y_index(Ys, Ms),
    yield_reduction = yield_reduction(Yp, Ys),
    R_drought_index = R_drought_index(Yp, Ys, Mp, Ms),
    Golden_mean = Golden_mean(Yp, Ys),
    ATI = ATI(Yp, Ys, Mp, Ms)
  )

  return(results)
}
## code to prepare `yield_data` data set goes here
#' @title Example data
#' @description Yield data of rice in kg/ha under two different growth conditions
#' @format A data frame with 50 rows and 3 variables:
#' \describe{
#'   \item{\code{Genotype}}{character Genotype}
#'   \item{\code{Yp}}{integer Yield under control condition}
#'   \item{\code{Ys}}{integer Yield under drought condition}
#' }
#' @source Simulated data, no external source were used
#' @references No external reference
"yield_data"
