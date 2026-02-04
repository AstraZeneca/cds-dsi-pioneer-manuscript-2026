# MCP Tool Definitions for Domino API

# Tool definitions following MCP schema
const TOOLS = [
    Dict(
        "name" => "start_job",
        "description" => "Start a new job in Domino. Jobs execute commands in the project's environment.",
        "inputSchema" => Dict(
            "type" => "object",
            "properties" => Dict(
                "project_id" => Dict(
                    "type" => "string",
                    "description" => "The Domino project ID (UUID format)"
                ),
                "command" => Dict(
                    "type" => "string",
                    "description" => "Command to execute (e.g., 'python train.py', 'Rscript analysis.R')"
                ),
                "title" => Dict(
                    "type" => "string",
                    "description" => "Optional title for the job"
                ),
                "hardware_tier_id" => Dict(
                    "type" => "string",
                    "description" => "Optional hardware tier ID for compute resources (uses overrideHardwareTierId in API)"
                ),
                "environment_id" => Dict(
                    "type" => "string",
                    "description" => "Optional environment ID to use"
                ),
                "git_ref" => Dict(
                    "type" => "string",
                    "description" => "Optional git branch name to use (e.g., 'main', 'karim/process-noise')"
                )
            ),
            "required" => ["project_id", "command"]
        )
    ),
    Dict(
        "name" => "stop_job",
        "description" => "Stop a running job in Domino. Can optionally commit results before stopping.",
        "inputSchema" => Dict(
            "type" => "object",
            "properties" => Dict(
                "job_id" => Dict(
                    "type" => "string",
                    "description" => "The job ID to stop"
                ),
                "commit_results" => Dict(
                    "type" => "boolean",
                    "description" => "Whether to commit results before stopping (default: true)"
                )
            ),
            "required" => ["job_id"]
        )
    ),
    Dict(
        "name" => "get_job_status",
        "description" => "Get the status and details of a specific job. Returns execution status, timestamps, and other metadata.",
        "inputSchema" => Dict(
            "type" => "object",
            "properties" => Dict(
                "job_id" => Dict(
                    "type" => "string",
                    "description" => "The job ID to check"
                )
            ),
            "required" => ["job_id"]
        )
    ),
    Dict(
        "name" => "list_jobs",
        "description" => "List jobs in a Domino project. Can filter by status and paginate results.",
        "inputSchema" => Dict(
            "type" => "object",
            "properties" => Dict(
                "project_id" => Dict(
                    "type" => "string",
                    "description" => "The Domino project ID"
                ),
                "status" => Dict(
                    "type" => "string",
                    "description" => "Filter by job status. IMPORTANT: Only 'running' (lowercase) is supported. Other status values will cause API errors."
                ),
                "page_size" => Dict(
                    "type" => "integer",
                    "description" => "Number of jobs per page (default: 20)"
                ),
                "page_no" => Dict(
                    "type" => "integer",
                    "description" => "Page number (default: 1)"
                )
            ),
            "required" => ["project_id"]
        )
    ),
    Dict(
        "name" => "list_hardware_tiers",
        "description" => "List available hardware tiers for a Domino project. Hardware tiers define compute resources (CPU, memory, GPU) available for jobs.",
        "inputSchema" => Dict(
            "type" => "object",
            "properties" => Dict(
                "project_id" => Dict(
                    "type" => "string",
                    "description" => "The Domino project ID (UUID format)"
                )
            ),
            "required" => ["project_id"]
        )
    )
]

# Server info for MCP
const SERVER_INFO = Dict(
    "name" => "domino-mcp",
    "version" => "0.1.0"
)
