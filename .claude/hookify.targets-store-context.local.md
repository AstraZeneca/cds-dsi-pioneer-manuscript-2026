---
name: targets-store-context
enabled: true
event: prompt
conditions:
  - field: user_prompt
    operator: regex_match
    pattern: store|tar_make|tar_progress|tar_meta|tar_outdated|pipeline.*progress|progress.*pipeline|targets.*store|store.*targets|main store|dco\d|TAR_BRANCH|TAR_PROJECT
---

**Targets store context check required.**

Before using any targets MCP tools (`tar_progress`, `tar_meta`, `tar_outdated`, `tar_read`, etc.), always verify:

1. **Which project?** — Do NOT default to sclc. Ask or infer from context:
   - `TAR_PROJECT=sclc` → store at `/mnt/data/analysis-results/<user>/sclc/<TAR_BRANCH>/`
   - `TAR_PROJECT=pioneer` → store at `/mnt/data/analysis-results/<user>/pioneer/<TAR_BRANCH>/`

2. **Which store (TAR_BRANCH)?** — If the user says "main store", "dco3", etc., map it to the correct path. If unclear, run:
   ```bash
   ls /mnt/data/analysis-results/karim_naguib/<project>/
   ```
   to list available stores before guessing.

3. **Never hardcode sclc** — the codebase supports multiple projects. Always confirm project before constructing a store path.
