#!/usr/bin/env julia
# Domino MCP Server
# MCP (Model Context Protocol) server for Domino Data Lab API

# CRITICAL: Handle initialize request IMMEDIATELY before any heavy loading
# MCP has a very short timeout for server startup

# Check if first message is initialize and respond immediately
first_line = readline(stdin)
if occursin(r"\"method\"\s*:\s*\"initialize\"", first_line)
    # Fast path: respond to initialize immediately (hardcoded response)
    response = """{"jsonrpc":"2.0","id":0,"result":{"protocolVersion":"2024-11-05","capabilities":{"tools":{}},"serverInfo":{"name":"domino-mcp","version":"0.1.0"}}}"""
    println(stdout, response)
    flush(stdout)
    global INIT_HANDLED = true
else
    # Not initialize - save for later processing
    global FIRST_LINE = first_line
    global INIT_HANDLED = false
end

# Now do heavy initialization
using JSON3

# Get script directory
const SCRIPT_DIR = @__DIR__

# Include modules
include(joinpath(SCRIPT_DIR, "tools.jl"))
include(joinpath(SCRIPT_DIR, "domino_client.jl"))

# Handle MCP initialize request
function handle_initialize(params)
    Dict(
        "protocolVersion" => "2024-11-05",
        "capabilities" => Dict(
            "tools" => Dict()
        ),
        "serverInfo" => SERVER_INFO
    )
end

# Handle tools/list request
function handle_tools_list()
    Dict("tools" => TOOLS)
end

# Helper to get value from NamedTuple or Dict
function safe_get(x, key::Symbol, default)
    if x isa NamedTuple
        haskey(x, key) ? getfield(x, key) : default
    elseif x isa Dict
        get(x, key, get(x, string(key), default))
    else
        default
    end
end

# Handle tools/call request
function handle_tools_call(params)
    tool_name = get(params, "name", "")
    tool_args = get(params, "arguments", Dict())

    result = try
        call_tool(tool_name, tool_args)
    catch e
        Dict(
            "success" => false,
            "error" => sprint(showerror, e)
        )
    end

    # Format result for MCP
    if safe_get(result, :success, false)
        Dict(
            "content" => [
                Dict(
                    "type" => "text",
                    "text" => JSON3.write(result)
                )
            ]
        )
    else
        Dict(
            "content" => [
                Dict(
                    "type" => "text",
                    "text" => "Error: $(safe_get(result, :error, "Unknown error"))"
                )
            ],
            "isError" => true
        )
    end
end

# Dispatch tool calls to appropriate functions
function call_tool(tool_name::String, args::Dict)
    if tool_name == "start_job"
        project_id = args["project_id"]
        command = args["command"]
        title = get(args, "title", nothing)
        hardware_tier_id = get(args, "hardware_tier_id", nothing)
        environment_id = get(args, "environment_id", nothing)
        git_ref = get(args, "git_ref", nothing)

        start_job(project_id, command;
                  title=title,
                  hardware_tier_id=hardware_tier_id,
                  environment_id=environment_id,
                  git_ref=git_ref)

    elseif tool_name == "stop_job"
        job_id = args["job_id"]
        commit_results = get(args, "commit_results", true)

        stop_job(job_id; commit_results=commit_results)

    elseif tool_name == "get_job_status"
        job_id = args["job_id"]

        get_job_status(job_id)

    elseif tool_name == "list_jobs"
        project_id = args["project_id"]
        status = get(args, "status", nothing)
        page_size = get(args, "page_size", 20)
        page_no = get(args, "page_no", 1)

        list_jobs(project_id;
                  status=status,
                  page_size=page_size,
                  page_no=page_no)

    elseif tool_name == "list_hardware_tiers"
        project_id = args["project_id"]

        list_hardware_tiers(project_id)

    else
        error("Unknown tool: $tool_name")
    end
end

# Process a single JSON-RPC request
function process_request(request::Dict)
    method = get(request, "method", "")
    params = get(request, "params", Dict())
    id = get(request, "id", nothing)

    result = try
        if method == "initialize"
            handle_initialize(params)
        elseif method == "tools/list"
            handle_tools_list()
        elseif method == "tools/call"
            handle_tools_call(params)
        elseif method == "notifications/initialized"
            nothing  # No response needed for notifications
        elseif method == "ping"
            Dict()  # Empty response for ping
        else
            Dict(
                "code" => -32601,
                "message" => "Unknown method: $method"
            )
        end
    catch e
        Dict(
            "code" => -32603,
            "message" => sprint(showerror, e)
        )
    end

    # Don't send response for notifications
    if isnothing(id)
        return nothing
    end

    # Build response
    if haskey(result, "code")
        # Error response
        Dict(
            "jsonrpc" => "2.0",
            "id" => id,
            "error" => result
        )
    else
        # Success response
        Dict(
            "jsonrpc" => "2.0",
            "id" => id,
            "result" => result
        )
    end
end

# Main loop - read from stdin, write to stdout
function main()
    # Check if we have a buffered first line that wasn't initialize
    pending_line = @isdefined(FIRST_LINE) ? FIRST_LINE : nothing

    while true
        # Get next line (from buffer or stdin)
        line = if !isnothing(pending_line)
            tmp = pending_line
            pending_line = nothing
            tmp
        elseif !eof(stdin)
            readline(stdin)
        else
            break
        end

        if isempty(strip(line))
            continue
        end

        # Parse JSON-RPC request
        request = try
            JSON3.read(line, Dict)
        catch e
            @error "Failed to parse JSON" exception=e line=line
            continue
        end

        # Skip initialize if we already handled it at startup
        method = get(request, "method", "")
        if method == "initialize" && @isdefined(INIT_HANDLED) && INIT_HANDLED
            continue
        end

        # Process request
        response = process_request(request)

        if !isnothing(response)
            # Write response to stdout
            println(stdout, JSON3.write(response))
            flush(stdout)
        end
    end
end

# Run main
main()
