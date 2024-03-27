test_that("chl_a function calculates chlorophyll-a content correctly", {
  # Define test inputs
  A663 <- 0.5
  A645 <- 0.3
  v <- 1.5
  w <- 2.0

  expected_result <- ((12.7 * A663) - (2.69 * A645)) * v / (1000 * w)

  result <- chl_a(A663, A645, v, w)
  expect_equal(result, expected_result, tolerance = 1e-10)
})

test_that("chl_b function calculates chlorophyll-b content correctly", {
  A645 <- 0.3
  A663 <- 0.5
  v <- 1.5
  w <- 2.0
  expected_result <- ((22.9 * A645) - (4.68 * A663)) * v / (1000 * w)
  result <- chl_b(A645, A663, v, w)
  expect_equal(result, expected_result, tolerance = 1e-10)
})
test_that("chl_total function calculates total chlorophyll content correctly", {
  A645 <- 0.3
  A663 <- 0.5
  v <- 1.5
  w <- 2.0
  expected_result <- ((20.2 * A645) + (8.02 * A663)) * v / (1000 * w)
  result <- chl_total(A645, A663, v, w)
  expect_equal(result, expected_result, tolerance = 1e-10)
})

test_that("caro_total function calculates total carotenoid concentration correctly", {
  A470 <- 0.4
  A663 <- 0.6
  A646 <- 0.5
  fresh_weight <- 10
  expected_result <- (1000 * A470 - 1.82 * (A663 + A646) + 85.02) / 198
  result <- caro_total(A470, A663, A646, fresh_weight)
  expect_equal(result, expected_result, tolerance = 1e-10)
})



test_that("peg_6000 function calculates PEG6000 requirement correctly", {
  peg <- 100
  C <- 25
  expected_result <- -(1.18 * 10^-2) * peg - (1.18 * 10^-4) * (peg^2) + (2.67 * 10^-4) * (peg * C) + (8.39 * 10^-7) * (peg^2) * (C)
  result <- peg_6000(peg, C)
  expect_equal(result, expected_result, tolerance = 1e-10)
})

test_that("calculate_PEG_6000 function calculates PEG6000 requirement correctly", {
  C <- 25
  bar <- -10
  result <- calculate_PEG_6000(C, bar)
  expected_PEG <- PEG <- ((0.0118 - 0.000267 * C) - ((0.000267 * C - 0.0118) ^ 2 + 4 * bar * (0.000000839 * C - 0.000118)) ^ 0.5) / (2 * (0.000000839 * C - 0.000118))
  expected_result <- list("gram/liter of water" = expected_PEG, "gram/cc of water" = expected_PEG / 1000)
  expect_equal(result, expected_result, tolerance = 1e-10)
})


test_that("Test st_index function", {
  expect_equal(st_index(0.5, 0.3), 0.6)
})

test_that("Test mp_index function", {
  expect_equal(mp_index(0.5, 0.3), 0.4)
})

test_that("Test gmp function", {
  expect_equal(gmp(0.5, 0.3), sqrt(0.5 * 0.3))
})

test_that("Test tol_index function", {
  expect_equal(tol_index(0.5, 0.3), 0.2)
})

test_that("Test ss_index function", {
  expect_equal(ss_index(0.5, 0.3, 0.2, 0.4), ((1 - 0.3) / 0.5) / ((1 - 0.2) / 0.4))
})

test_that("Test YSI function", {
  expect_equal(YSI(0.5, 0.3), 0.6)
})

test_that("Test YR_ratio function", {
  expect_equal(YR_ratio(0.5, 0.3), 1 - (0.3 / 0.5))
})

test_that("Test DRI function", {
  expect_equal(DRI(0.5, 0.3), (0.3 * (0.3 / 0.5) / 0.5))
})

test_that("Test HAM function", {
  expect_equal(HAM(0.5, 0.3), 2 * (0.5 * 0.3) / (0.5 + 0.3))
})

test_that("Test Y_index function", {
  expect_equal(Y_index(0.3, 0.2), 0.3 / 0.2)
})

test_that("Test yield_reduction function", {
  expect_equal(yield_reduction(0.5, 0.3), ((0.5 - 0.3) / 0.3) * 100)
})

test_that("Test R_drought_index function", {
  expect_equal(R_drought_index(0.5, 0.3, 0.4, 0.2), (0.3 / 0.3) / (0.4 / 0.2))
})

test_that("Test Golden_mean function", {
  expect_equal(Golden_mean(0.5, 0.3), (0.5 + 0.3) / (0.5 - 0.3))
})

test_that("Test ATI function", {
  expect_equal(ATI(0.5, 0.3, 0.4, 0.2), ((0.5 - 0.3) / (0.4 / 0.2)) * sqrt(0.5 * 0.3))
})

test_that("Test all_indices function", {
  Yp <- 0.5
  Ys <- 0.3
  Mp <- 0.4
  Ms <- 0.2

  # Call the all_indices function
  result <- all_indices(Yp, Ys, Mp, Ms)

  # Compare the result with the expected values
  expected_results <- list(
    StressToleranceIndex = (Ys * Yp) / (Yp)^2,
    MeanProductivityIndex = (Ys + Yp) / 2,
    GeometricMeanProductivity = sqrt(Ys * Yp),
    Tolerance_Index = (Yp - Ys),
    StressSusceptibilityIndex = ((1 - Ys) / Yp) / ((1 - Ms) / Mp),
    YSI = Ys / Yp,
    YR_Ratio = 1 - (Ys / Yp),
    DRI = Ys * (Ys / Yp) / Yp,
    HAM = 2 * (Yp * Ys) / (Yp + Ys),
    Y_Index = Ys / Ms,
    yield_reduction = ((Yp - Ys) / Ys) * 100,
    R_drought_index = (Ys / Ys) / (Mp / Ms),
    Golden_mean = (Yp + Ys) / (Yp - Ys),
    ATI = ((Yp - Ys) / (Mp / Ms)) * sqrt(Yp * Ys)
  )

  # Compare each result with the expected value
  expect_equal(result$StressToleranceIndex, expected_results$StressToleranceIndex)
  expect_equal(result$MeanProductivityIndex, expected_results$MeanProductivityIndex)
  expect_equal(result$GeometricMeanProductivity, expected_results$GeometricMeanProductivity)
  expect_equal(result$Tolerance_Index, expected_results$Tolerance_Index)
  expect_equal(result$StressSusceptibilityIndex, expected_results$StressSusceptibilityIndex)
  expect_equal(result$YSI, expected_results$YSI)
  expect_equal(result$YR_Ratio, expected_results$YR_Ratio)
  expect_equal(result$DRI, expected_results$DRI)
  expect_equal(result$HAM, expected_results$HAM)
  expect_equal(result$Y_Index, expected_results$Y_Index)
  expect_equal(result$yield_reduction, expected_results$yield_reduction)
  expect_equal(result$R_drought_index, expected_results$R_drought_index)
  expect_equal(result$Golden_mean, expected_results$Golden_mean)
  expect_equal(result$ATI, expected_results$ATI)
})
