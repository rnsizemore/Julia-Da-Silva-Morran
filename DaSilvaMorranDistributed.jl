using Distributed
ntasks = parse(Int,ARGS[1]) # This is the number of cores over which the number of replicates will be distributed

# Add worker processes. If on the HPC, use slurm, else use just regular CPU thread workers
if "SLURM_JOB_ID" in keys(ENV) && ENV["SLURM_SUBMIT_HOST"] != "kimura"
    @info "Adding slurm procs"
    using ClusterManagers
    addprocs_slurm(ntasks, unbuffered="")
else
    @info "Adding local procs"
    addprocs(ntasks)
end

# Create queues for jobs and results
const jobqueue = RemoteChannel(() -> Channel{Tuple}(ntasks))
const resultsqueue = RemoteChannel(() -> Channel{NamedTuple}(ntasks))

# This is the main loop for the workers. They take a job from the queue and runs
# it. The result is put in the results queue. take! operations will block if
# there are no jobs available, put! operations block if the results queue is
# full.
@everywhere function dowork(jobs, results)
    while true
        # fn, args... = take!(jobs)     # only available in 1.6
        job = take!(jobs)               # Take the next available job
        fn, args = job[1], job[2:end]
        result = eval(fn)(args...)      # Run the job
        put!(results, result)           # Send the result
    end
end

# Starts event loop on all workers
start_workers() = foreach(pid -> remote_do(dowork, pid, jobqueue, resultsqueue), workers())

# Gracefully stops workers
stop_workers() = rmprocs(workers())

# Submits a single job. Blocks if jobqueue is full
submit_job(job) = put!(jobqueue, job)

# Launches an async task to submit many jobs. Does not block the main thread
#submit_jobs(jobs) = @async foreach(submit_job, jobs)
submit_jobs(jobs) = Threads.@spawn foreach(submit_job, jobs)