# Domino MCP Server

MCP (Model Context Protocol) server for Domino Data Lab API, providing job management tools for AI agents.

## Features

- **start_job** - Start a new job in Domino
- **stop_job** - Stop a running job
- **get_job_status** - Get job status and details
- **list_jobs** - List jobs in a project with filtering
- **list_hardware_tiers** - List available hardware tiers for a project

## Prerequisites

- Julia 1.12+ installed
- Domino Data Lab API key
- Access to a Domino Data Lab instance

## Setup

### 1. Install Dependencies

```bash
cd .claude/mcp/domino-mcp
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### 2. Configure Environment Variables

Export your Domino credentials:

```bash
export DOMINO_USER_API_KEY="your-api-key-here"
export DOMINO_HOST="https://your-domino-instance.com"
```

To get your API key:
1. Log into Domino
2. Go to Account Settings
3. Navigate to API Keys section
4. Generate a new key or copy existing one

### 3. Add to Claude Code

The server is already configured in `.mcp.json`. When you restart Claude Code, it will automatically connect to the Domino MCP server.

## Usage

Once configured, you can use natural language to interact with Domino:

- "Start a job in project abc123 that runs 'python train.py'"
- "What's the status of job xyz789?"
- "List all running jobs in project abc123"
- "Stop job xyz789"
- "What hardware tiers are available for project abc123?"

## Testing

Test the server manually:

```bash
# Initialize
echo '{"jsonrpc":"2.0","id":1,"method":"initialize","params":{}}' | julia --project=. server.jl

# List tools
echo '{"jsonrpc":"2.0","id":2,"method":"tools/list","params":{}}' | julia --project=. server.jl
```

## Architecture

Single-process Julia application:
- `server.jl` - MCP protocol handler (stdin/stdout JSON-RPC)
- `domino_client.jl` - HTTP client for Domino API v4
- `tools.jl` - MCP tool definitions

## Troubleshooting

**Server not starting:**
- Check Julia is installed: `julia --version`
- Verify dependencies installed: `julia --project=. -e 'using HTTP, JSON3'`

**Authentication errors:**
- Verify `DOMINO_USER_API_KEY` and `DOMINO_HOST` are set
- Check API key is valid in Domino UI
- Ensure host URL includes protocol (https://) but no trailing slash

**Job operations failing:**
- Verify project IDs are correct UUIDs
- Check you have permissions in Domino
- Review error messages in MCP tool responses
