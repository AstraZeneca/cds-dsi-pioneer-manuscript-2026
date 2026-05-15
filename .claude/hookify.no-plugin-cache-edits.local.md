---
name: no-plugin-cache-edits
enabled: true
event: file
pattern: /home/ubuntu/\.claude/plugins/
action: block
---

⛔ **Plugin cache is read-only**

You're trying to edit a file inside the plugin cache:
`/home/ubuntu/.claude/plugins/`

This directory is managed by the plugin manager and gets overwritten on reinstall.

**Make the edit in the marketplace repo instead:**
`/home/ubuntu/claude-marketplace/plugins/<plugin-name>/...`

Then commit and push from `/home/ubuntu/claude-marketplace/`.
The plugin manager will pick up changes on next install/update.
