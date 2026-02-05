# HTTP client for communicating with R plumber backend

using HTTP
using JSON3

const PLUMBER_HOST = "127.0.0.1"
const PLUMBER_PORT = 8484
const PLUMBER_URL = "http://$PLUMBER_HOST:$PLUMBER_PORT"

# Start the plumber server if not running
function ensure_plumber_running()
    # Check if server is healthy
    try
        response = HTTP.get("$PLUMBER_URL/health"; connect_timeout=2, readtimeout=2)
        if response.status == 200
            return true
        end
    catch e
        # Server not responding, need to start it
    end

    # Start plumber in background
    script_dir = @__DIR__
    start_script = joinpath(script_dir, "start_server.R")

    # Run in background
    run(pipeline(`Rscript $start_script $PLUMBER_PORT`; stdout=devnull, stderr=devnull); wait=false)

    # Wait for server to start
    max_attempts = 20
    for i in 1:max_attempts
        sleep(0.5)
        try
            response = HTTP.get("$PLUMBER_URL/health"; connect_timeout=1, readtimeout=1)
            if response.status == 200
                return true
            end
        catch
            # Keep waiting
        end
    end

    error("Failed to start plumber server after $max_attempts attempts")
end

# Endpoints that use GET (read-only, no body parameters)
const GET_ENDPOINTS_LOCAL = ["list_pipelines", "list_sessions", "health"]

# Call a plumber endpoint
function call_plumber(endpoint::String, params::Dict=Dict())
    url = "$PLUMBER_URL/$endpoint"

    try
        is_get = endpoint in GET_ENDPOINTS_LOCAL
        response = if is_get
            HTTP.get(url; connect_timeout=30, readtimeout=120)
        else
            HTTP.post(url,
                ["Content-Type" => "application/json"],
                JSON3.write(params);
                connect_timeout=30,
                readtimeout=120
            )
        end

        body = String(response.body)
        JSON3.read(body, Dict)
    catch e
        Dict("success" => false, "error" => sprint(showerror, e))
    end
end
