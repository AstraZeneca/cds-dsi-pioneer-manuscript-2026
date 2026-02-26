# GitHub Project Management

Issues across all PIONEER repos are tracked in the **PIONEER 2026** GitHub Project.

## Project IDs (reference)
- **Project number**: 56
- **Project node ID**: `PVT_kwDOCLUTRM4A7VlR`
- **Owner**: `azu-oncology-rd`

## Field IDs
| Field  | Field ID | Options |
|--------|----------|---------|
| Status | `PVTSSF_lADOCLUTRM4A7VlRzgvqnuE` | Todo: `f75ad846`, In Progress: `47fc9ee4`, Done: `98236657`, Blocked: `262596ac` |
| Trial  | `PVTSSF_lADOCLUTRM4A7VlRzg76fW0` | SCLC-01: `a7d09f93`, Trial-Lung-01: `48351c43`, TRIAL-C: `b01a613f`, TRIAL-D: `706eb7ec`, Pioneer: `22aaa53c` |

## Creating an Issue and Adding to Project

```bash
# 1. Create the issue in the appropriate repo (assigned to current user)
gh issue create \
  --repo azu-oncology-rd/REPO_NAME \
  --assignee kmjq089_azu \
  --title "Issue title" \
  --body "Issue body"

# 2. Add issue to the project
gh project item-add 56 --owner azu-oncology-rd --url ISSUE_URL

# 3. Find the project item ID
gh project item-list 56 --owner azu-oncology-rd --limit 200 --format json \
  | python3 -c "
import json, sys
data = json.load(sys.stdin)
for item in data.get('items', []):
    if 'SEARCH_TERM' in item.get('title', ''):
        print(item['id']); break
"

# 4. Set Trial field (use option IDs from table above)
gh project item-edit \
  --project-id PVT_kwDOCLUTRM4A7VlR \
  --id ITEM_ID \
  --field-id PVTSSF_lADOCLUTRM4A7VlRzg76fW0 \
  --single-select-option-id TRIAL_OPTION_ID

# 5. Set Status field (use option IDs from table above)
gh project item-edit \
  --project-id PVT_kwDOCLUTRM4A7VlR \
  --id ITEM_ID \
  --field-id PVTSSF_lADOCLUTRM4A7VlRzgvqnuE \
  --single-select-option-id STATUS_OPTION_ID
```

## Auth Requirements

The `gh` token needs the `project` scope. If missing, run:
```bash
gh auth refresh -h github.com -s project
```

## Querying Field Options

If new trials or statuses are added, query the current options:
```bash
gh api graphql -f query='
{
  node(id: "PVT_kwDOCLUTRM4A7VlR") {
    ... on ProjectV2 {
      field(name: "FIELD_NAME") {
        ... on ProjectV2SingleSelectField {
          options { id name }
        }
      }
    }
  }
}'
```
