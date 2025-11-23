test_that("Test cluster distance matrix", {
    for (n_clust in 2:10) {
        # Constructs Matrix
        diss_clust <- matrix(
            runif(min = 0.3, max = 1, n_clust * n_clust),
            nrow = n_clust,
            ncol = n_clust
        )

        # Converts to distance
        diss_clust <- (diss_clust + t(diss_clust)) / 2
        diag(diss_clust) <- 0
        diss_clust <- round(diss_clust, 2)

        as.dist(diss_clust)

        clust_ids_order <- seq_len(n_clust)

        # Generates the cluster ids
        clust <- c()
        for (i in clust_ids_order) {
            clust <- c(clust, rep(i, runif(1, min = 1, max = 10)))
        }

        n_elem <- length(clust)
        n_clust <- length(clust_ids_order)

        # Generates the distance
        diss_avg <- matrix(0, nrow = n_elem, ncol = n_elem)
        diss_max <- matrix(0, nrow = n_elem, ncol = n_elem)
        diss_min <- matrix(2, nrow = n_elem, ncol = n_elem)
        diag(diss_min) <- 0

        for (i in seq_len(n_clust - 1)) {
            for (j in (i + 1):n_clust) {
                # Average
                diss_avg[which(i == clust), which(j == clust)] <- diss_clust[
                    i,
                    j
                ]
                diss_avg[which(j == clust), which(i == clust)] <- diss_clust[
                    i,
                    j
                ]

                # Distorts
                if (sum(i == clust) > 1 && sum(j == clust) > 1) {
                    noise <- runif(1, min = 0.1, max = 0.3)
                    diss_avg[
                        which(i == clust)[1],
                        which(j == clust)[1]
                    ] <- diss_clust[i, j] - noise
                    diss_avg[
                        which(j == clust)[1],
                        which(i == clust)[1]
                    ] <- diss_clust[i, j] - noise
                    diss_avg[
                        which(i == clust)[2],
                        which(j == clust)[2]
                    ] <- diss_clust[i, j] + noise
                    diss_avg[
                        which(j == clust)[2],
                        which(i == clust)[2]
                    ] <- diss_clust[i, j] + noise
                }

                # Max
                diss_max[
                    which(i == clust)[1],
                    which(j == clust)[1]
                ] <- diss_clust[
                    i,
                    j
                ]
                diss_max[
                    which(j == clust)[1],
                    which(i == clust)[1]
                ] <- diss_clust[
                    i,
                    j
                ]

                # Min
                diss_min[
                    which(i == clust)[1],
                    which(j == clust)[1]
                ] <- diss_clust[
                    i,
                    j
                ]
                diss_min[
                    which(j == clust)[1],
                    which(i == clust)[1]
                ] <- diss_clust[
                    i,
                    j
                ]
            }
        }

        diss_avg <- as.dist(diss_avg)
        diss_min <- as.dist(diss_min)
        diss_max <- as.dist(diss_max)

        expect_true(all(
            as.dist(diss_clust) ==
                cluster_distance_matrix(
                    diss_avg,
                    clust,
                    clust_ids_order,
                    method = "average"
                )
        ))
        expect_true(all(
            as.dist(diss_clust) ==
                cluster_distance_matrix(
                    diss_min,
                    clust,
                    clust_ids_order,
                    method = "min"
                )
        ))
        expect_true(all(
            as.dist(diss_clust) ==
                cluster_distance_matrix(
                    diss_max,
                    clust,
                    clust_ids_order,
                    method = "max"
                )
        ))
    }
})
