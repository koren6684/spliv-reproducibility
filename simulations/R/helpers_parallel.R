spliv_sim_apply <- function(tasks, worker, n_cores = 1L) {
  n_cores <- max(1L, as.integer(n_cores))
  if (length(tasks) == 0L) {
    return(list())
  }

  if (n_cores <= 1L || .Platform$OS.type == "windows") {
    return(lapply(tasks, worker))
  }

  parallel::mclapply(
    tasks,
    worker,
    mc.cores = n_cores
  )
}
