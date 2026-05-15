# Publishing to RStudio Connect

**Server URL**: `https://rstudio-connect.seml.scp.astrazeneca.net/connect/`
**Account**: `kmjq089`

**IMPORTANT**: Always use the R `rsconnect` package for publishing. The `quarto publish connect` CLI command does not work reliably in non-interactive/remote environments because it requires browser-based SSO authentication.

## Step-by-Step Publishing Instructions

**IMPORTANT**: Always render the site locally first to ensure everything works. Run from the project root:
```bash
# From /mnt/code (project root)
quarto render quarto/website
```

### Method 1: Using R rsconnect (Recommended)

This is the most reliable method. Follow these steps:

**Step 1: Get your API key**
1. Go to `https://rstudio-connect.seml.scp.astrazeneca.net/connect/` in your browser
2. Sign in with SSO
3. Click your name (top right) → "API Keys"
4. Click "New API Key" and copy it

**Step 2: Configure credentials (one-time setup)**
```r
library(rsconnect)

# Add the server
rsconnect::addConnectServer(
  url = "https://rstudio-connect.seml.scp.astrazeneca.net",
  name = "az-connect"
)

# Add your API key (paste your actual key here)
rsconnect::connectApiUser(
  account = "kmjq089",
  server = "az-connect",
  apiKey = "YOUR_API_KEY_HERE"
)
```

This stores credentials in `~/.rsconnect/` (outside the git repo).

**Step 3: Deploy the website**

**IMPORTANT**: Always run from the project root directory (where renv is configured), not from inside the website directory.

```r
# From /mnt/code directory
library(rsconnect)
rsconnect::deploySite(
  siteDir = "quarto/website",
  server = "az-connect",
  account = "kmjq089"
)
```

**Step 4: Find your published site**
After deployment completes, look for the URL in the output or go to:
`https://rstudio-connect.seml.scp.astrazeneca.net/connect/#/content/listing?q=owner:kmjq089`

### Method 2: Using Quarto CLI

This requires browser authentication (SSO) and may not work in remote CLI environments.

```bash
# From /mnt/code (project root)
quarto publish connect quarto/website --server https://rstudio-connect.seml.scp.astrazeneca.net/connect/
```

This will open your browser for SSO authentication.

## Updating an Existing Deployment

Once you've published the first time, subsequent deployments are simple:

```r
# From /mnt/code directory
library(rsconnect)
rsconnect::deploySite(
  siteDir = "quarto/website",
  server = "az-connect",
  account = "kmjq089"
)
```

Or from bash:
```bash
# From /mnt/code (project root)
quarto publish connect quarto/website
```

## Troubleshooting

- **API key expired**: Generate a new one and re-run `rsconnect::connectApiUser()`
- **Publishing fails**: Make sure `quarto render quarto/website` (from project root) completes successfully first
- **Missing plots**: Delete `_freeze/` cache and re-render from project root
- **IMPORTANT**: Always run quarto commands from `/mnt/code` (project root), not from inside `quarto/website/`
- **Git tracking warnings**: Confirm `rsconnect/` and `_publish.yml` are in `.gitignore`
- **Can't find account/server**: Make sure you're using the correct server name (`az-connect`). Run `rsconnect::accounts()` to check configured accounts.
- **Package not found errors**: Always run deployment from `/mnt/code` (project root) where renv is configured

**SECURITY**: API keys should NEVER be committed to git. The `rsconnect/` directory and `_publish.yml` are already in `.gitignore`.

## Publishing Presentations

For Quarto revealjs presentations (e.g., `quarto/presentations/pioneer-gng/`), use `rsconnect::deployDoc()`:

```r
library(rsconnect)

# First render the presentation
# quarto render quarto/presentations/pioneer-gng/pioneer-gng.qmd

# Then deploy the rendered HTML
rsconnect::deployDoc(
  doc = "quarto/presentations/pioneer-gng/pioneer-gng.html",
  server = "az-connect",
  account = "kmjq089",
  appName = "pioneer-gng-presentation"  # Choose a unique name
)
```

**Published sites:**
- SCLC-01: `rsconnect::deploySite(siteDir = "quarto/sclc/website", server = "az-connect", account = "kmjq089")`
- Pioneer: `rsconnect::deploySite(siteDir = "quarto/pioneer/website", server = "az-connect", account = "kmjq089")`
  - URL: https://rstudio-connect.seml.scp.astrazeneca.net/content/4ec50122-5d33-43fa-831c-3df0243084f7/
  - siteName: `pioneer-pioneer`

**Published presentations:**
- PIONEER Go/No-Go: https://rstudio-connect.seml.scp.astrazeneca.net/content/71d3bcc6-c677-485b-b8fc-562b0f580c9a/
