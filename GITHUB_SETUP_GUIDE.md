# GitHub Setup Guide — Step by Step for R Projects

This guide walks you through linking your R project to GitHub from scratch.

---

## STEP 1 — Install Git on your computer

Download from: https://git-scm.com/downloads

Check it worked — open your terminal (Mac/Linux) or Git Bash (Windows):
```bash
git --version
# Should return something like: git version 2.43.0
```

---

## STEP 2 — Configure Git with your identity

Only do this once on your computer:
```bash
git config --global user.name "Your Name"
git config --global user.email "your@email.com"
```

---

## STEP 3 — Create a GitHub account & new repository

1. Go to https://github.com and sign in
2. Click the **+** icon (top right) → **New repository**
3. Name it: `fracture-risk-nhanes`
4. Set to **Public** (employers need to see it)
5. ✅ Tick "Add a README file" — NO, skip this (you already have one)
6. Click **Create repository**
7. Copy the URL shown — looks like: `https://github.com/yourusername/fracture-risk-nhanes.git`

---

## STEP 4 — Initialise Git in your local project folder

Open your terminal and navigate to your project:
```bash
cd path/to/fracture-risk-nhanes    # change this to where your folder is

git init                            # initialises git tracking
git add .                           # stages all files
git commit -m "Initial commit: data prep and T-score scripts"
```

---

## STEP 5 — Connect local repo to GitHub and push

```bash
git remote add origin https://github.com/yourusername/fracture-risk-nhanes.git
git branch -M main
git push -u origin main
```

Now your code is live on GitHub.

---

## STEP 6 — Your daily Git workflow going forward

Every time you make changes and want to save them to GitHub:

```bash
# Check what files have changed
git status

# Stage the files you want to commit (or use . for all)
git add scripts/01_data_prep.R
git add scripts/02_t_score_calc.R

# Commit with a meaningful message
git commit -m "Add rename() step for readable column names"

# Push to GitHub
git push
```

**Good commit message examples:**
```bash
git commit -m "Fix: drop NA target variable before factor conversion"
git commit -m "Add: purrr::reduce() for scalable dataset merging"
git commit -m "Refactor: consolidate library imports, remove duplicates"
git commit -m "Docs: add clinical justification comments to T-score script"
```

---

## STEP 7 — Using RStudio's built-in Git panel (easier!)

RStudio has a Git tab that does all of the above with clicks:

1. Open your project in RStudio
2. Go to **Tools** → **Version Control** → **Project Setup**
3. Select **Git** → click Yes to restart RStudio
4. You'll now see a **Git tab** in the top-right panel
5. Use the checkboxes to stage files, write commit messages, and click **Push**

To connect to your GitHub remote the first time via RStudio:
- Go to **Tools** → **Shell** and run the `git remote add origin` command from Step 5

---

## STEP 8 — Keep raw data out of GitHub

Your `.gitignore` file already handles this. Verify it's working:
```bash
git status
# data/raw/ and *.XPT files should NOT appear in the list
```

Instead, document where to get the data in your README (already done ✅).

---

## STEP 9 — Make your repo look professional

Things to do on the GitHub website after pushing:

1. **Add a description** — click the gear icon on your repo homepage
   - Description: *"Fracture risk prediction using NHANES 2017-2020 | R | DEXA | WHO T-scores"*
   - Website: link to RPubs report if you publish one
   - Topics (tags): `r`, `nhanes`, `public-health`, `data-science`, `osteoporosis`, `machine-learning`

2. **Pin the repo** to your GitHub profile
   - Go to your profile → click **Customize your pins**

3. **Write a good profile README** (optional but impressive)
   - Create a repo named exactly the same as your username
   - Add a README.md — it shows on your GitHub profile page

---

## STEP 10 — Publish your R Markdown report (optional but very impressive)

Once you write `03_eda.Rmd` (an R Markdown file), you can publish it:

**Option A — RPubs (easiest):**
- Click **Knit** in RStudio → HTML
- Click **Publish** → RPubs
- Free, instant, shareable link

**Option B — GitHub Pages:**
```bash
# Knit your .Rmd to HTML first, then:
git add reports/eda_report.html
git commit -m "Add: rendered EDA report"
git push
```
Then in GitHub → Settings → Pages → set source to main branch → `/reports` folder.

---

## Quick Reference Card

| Action | Command |
|---|---|
| Check status | `git status` |
| Stage all files | `git add .` |
| Commit | `git commit -m "message"` |
| Push to GitHub | `git push` |
| Pull latest from GitHub | `git pull` |
| See commit history | `git log --oneline` |
| Undo last commit (keep files) | `git reset --soft HEAD~1` |
