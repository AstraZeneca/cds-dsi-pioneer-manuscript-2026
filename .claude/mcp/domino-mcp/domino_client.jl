# Domino API Client
# HTTP client for Domino Data Lab API v4

using HTTP
using JSON3

# Get configuration from environment
function get_config()
    api_key = get(ENV, "DOMINO_USER_API_KEY", "")
    host = get(ENV, "DOMINO_HOST", "")

    if isempty(api_key)
        error("DOMINO_USER_API_KEY environment variable is not set")
    end
    if isempty(host)
        error("DOMINO_HOST environment variable is not set")
    end

    # Ensure host doesn't have trailing slash
    host = rstrip(host, '/')

    return (api_key=api_key, host=host)
end

# Base HTTP request function
function domino_request(method::Symbol, endpoint::String; body=nothing, query=nothing)
    config = get_config()

    url = "$(config.host)$(endpoint)"

    headers = [
        "X-Domino-Api-Key" => config.api_key,
        "Content-Type" => "application/json",
        "Accept" => "application/json"
    ]

    try
        if method == :GET
            if !isnothing(query)
                response = HTTP.get(url, headers; query=query)
            else
                response = HTTP.get(url, headers)
            end
        elseif method == :POST
            body_str = isnothing(body) ? "{}" : JSON3.write(body)
            response = HTTP.post(url, headers, body_str)
        else
            error("Unsupported HTTP method: $method")
        end

        return (
            success = true,
            status = response.status,
            data = JSON3.read(String(response.body))
        )
    catch e
        if e isa HTTP.StatusError
            error_body = try
                JSON3.read(String(e.response.body))
            catch
                String(e.response.body)
            end
            return (
                success = false,
                status = e.status,
                error = "HTTP $(e.status): $(error_body)"
            )
        else
            return (
                success = false,
                status = 0,
                error = string(e)
            )
        end
    end
end

# Start a new job
function start_job(project_id::String, command::String;
                   title::Union{String,Nothing}=nothing,
                   hardware_tier_id::Union{String,Nothing}=nothing,
                   environment_id::Union{String,Nothing}=nothing,
                   git_ref::Union{String,Nothing}=nothing)

    body = Dict{String,Any}(
        "projectId" => project_id,
        "commandToRun" => command
    )

    if !isnothing(title)
        body["title"] = title
    end
    if !isnothing(hardware_tier_id)
        # Use overrideHardwareTierId (not hardwareTierId) - this is what Domino UI uses
        body["overrideHardwareTierId"] = hardware_tier_id
    end
    if !isnothing(environment_id)
        body["environmentId"] = environment_id
        body["environmentRevisionSpec"] = "ActiveRevision"
    end
    if !isnothing(git_ref)
        # Format git ref as nested object with type and value
        body["mainRepoGitRef"] = Dict(
            "type" => "branches",
            "value" => git_ref
        )
    end

    result = domino_request(:POST, "/v4/jobs/start"; body=body)

    if result.success
        return (
            success = true,
            job_id = get(result.data, :id, get(result.data, :jobId, nothing)),
            data = result.data
        )
    else
        return result
    end
end

# Stop a running job
function stop_job(job_id::String; commit_results::Bool=true)
    body = Dict{String,Any}(
        "jobId" => job_id,
        "commitResults" => commit_results
    )

    result = domino_request(:POST, "/v4/jobs/stop"; body=body)

    if result.success
        return (
            success = true,
            message = "Job $job_id stop requested",
            data = result.data
        )
    else
        return result
    end
end

# Get job status/details
function get_job_status(job_id::String)
    result = domino_request(:GET, "/v4/jobs/$job_id")

    if result.success
        data = result.data
        statuses = get(data, :statuses, Dict())
        return (
            success = true,
            job_id = job_id,
            status = get(statuses, :executionStatus, "Unknown"),
            data = data
        )
    else
        return result
    end
end

# List jobs in a project
function list_jobs(project_id::String;
                   status::Union{String,Nothing}=nothing,
                   page_size::Int=20,
                   page_no::Int=1)

    query = Dict{String,String}(
        "projectId" => project_id,
        "pageSize" => string(page_size),
        "pageNo" => string(page_no)
    )

    if !isnothing(status)
        query["status"] = status
    end

    result = domino_request(:GET, "/v4/jobs"; query=query)

    if result.success
        raw_jobs = get(result.data, :jobs, result.data)

        # Extract only essential fields to reduce token usage
        jobs = map(raw_jobs) do job
            statuses = get(job, :statuses, Dict())
            stage_time = get(job, :stageTime, Dict())
            started_by = get(job, :startedBy, Dict())
            git_ref = get(job, :mainRepoGitRef, Dict())
            hardware_tier = get(job, :hardwareTier, Dict())

            Dict(
                :number => get(job, :number, 0),
                :id => get(job, :id, ""),
                :command => get(job, :jobRunCommand, ""),
                :title => get(job, :title, nothing),
                :status => get(statuses, :executionStatus, ""),
                :is_completed => get(statuses, :isCompleted, false),
                :branch => get(git_ref, :value, ""),
                :started_by => get(started_by, :username, ""),
                :submission_time => get(stage_time, :submissionTime, nothing),
                :run_start_time => get(stage_time, :runStartTime, nothing),
                :completed_time => get(stage_time, :completedTime, nothing),
                :hardware_tier_name => get(hardware_tier, :name, ""),
                :hardware_tier_id => get(hardware_tier, :id, "")
            )
        end

        return (
            success = true,
            count = length(jobs),
            jobs = jobs
        )
    else
        return result
    end
end

# List hardware tiers available for a project
function list_hardware_tiers(project_id::String)
    result = domino_request(:GET, "/v4/projects/$project_id/hardwareTiers")

    if result.success
        # The response is typically an array directly, or nested under a key
        raw_tiers = if result.data isa AbstractArray
            result.data
        else
            get(result.data, :hardwareTiers, get(result.data, :data, []))
        end

        # Extract only essential fields to reduce token usage
        hardware_tiers = map(raw_tiers) do tier_obj
            hwt = get(tier_obj, :hardwareTier, tier_obj)
            resources = get(hwt, :hwtResources, Dict())
            flags = get(hwt, :hwtFlags, Dict())
            gpu_config = get(hwt, :gpuConfiguration, Dict())
            memory = get(resources, :memory, Dict())

            Dict(
                :id => get(hwt, :id, ""),
                :name => get(hwt, :name, ""),
                :cores => get(resources, :cores, 0),
                :memory_gib => get(memory, :value, 0),
                :gpus => get(gpu_config, :numberOfGpus, 0),
                :cost_per_minute => get(hwt, :centsPerMinute, 0) / 100.0,
                :is_default => get(flags, :isDefault, false)
            )
        end

        return (
            success = true,
            count = length(hardware_tiers),
            hardware_tiers = hardware_tiers
        )
    else
        return result
    end
end
