
test_that("Formula matches q != 1", {


  for(i in 1:15)
  {

    q  <-  runif(1, 0.001,10)
    if(q == 1)
      q  <- q + 1e-10

    # Three cluster
    ab_1 <- runif(4,1,100)
    ab_2 <- runif(6,1,100)
    ab_3 <- runif(15,1,100)

    # By formula
    total <- sum(ab_1) + sum(ab_2) + sum(ab_3)
    P_1 <- ab_1 / sum(ab_1)
    P_2 <- ab_2 / sum(ab_2)
    P_3 <- ab_3 / sum(ab_3)

    numenator <-  (sum(ab_1)**q)*sum(P_1**q) + (sum(ab_2)**q)*sum(P_2**q)  + (sum(ab_3)**q)*sum(P_3**q)
    denominator <- sum(ab_1)**q + sum(ab_2)**q  + sum(ab_3)**q

    by_formula <- (numenator / denominator)**(1/(1-q))

    ab <- c(ab_1, ab_2, ab_3)
    clust <- c(rep(1, length(ab_1)), rep(2, length(ab_2)), rep(3, length(ab_3)))
    by_implementation  <- multiplicity.inventory(ab, clust, q)

    expect_equal(by_formula, by_implementation)

  }



})


test_that("Formula matches q = 1", {


  q <- 1
  for(i in 1:15)
  {

    # Three cluster
    ab_1 <- runif(4,1,100)
    ab_2 <- runif(6,1,100)
    ab_3 <- runif(15,1,100)

    # By formula
    N <- sum(ab_1) + sum(ab_2) + sum(ab_3)
    N_1 <- sum(ab_1)
    N_2 <- sum(ab_2)
    N_3 <- sum(ab_3)

    ab <- c(ab_1, ab_2, ab_3)

    first <- (1/N)*(N_1*log(N_1) + N_2*log(N_2) + N_3*log(N_3))
    second <- (1/N)*sum(ab*log(ab))


    by_formula <- exp((first-second))

    clust <- c(rep(1, length(ab_1)), rep(2, length(ab_2)), rep(3, length(ab_3)))
    by_implementation  <- multiplicity.inventory(ab, clust, q)

    expect_equal(by_formula, by_implementation)

  }


})



test_that("High Multlicity Works", {

  ab <- rep(10,30)
  clust <- c(rep(1,10), rep(2,10), rep(3,10))
  # Multiplicity should be 10
  expect_equal(multiplicity.inventory(ab, clust), 10)
})


test_that("Low Multlicity Works", {

  ab <- rep(10,3)
  clust <- c(1,2,3)
  # Multiplicity should be 1
  expect_equal(multiplicity.inventory(ab, clust), 1)
})


test_that("Abundance with zero works", {

  ab_1 <- c(10,3,3,5,7,1,2,3)
  clust_1 <- c(1,2,3,2,3,2,1,2)

  ab_2 <- c(10,3,3,5,7,1,2,3,0,0,0,0)
  clust_2 <- c(1,2,3,2,3,2,1,2,1,2,3,1)

  # Multiplicity should be equal
  expect_equal(multiplicity.inventory(ab_1, clust_1), multiplicity.inventory(ab_2, clust_2))
})


