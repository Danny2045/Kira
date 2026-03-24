#!/bin/bash
# ============================================================
# PUSH KIRA TO GITHUB
# ============================================================
# Run this from ~/kira after:
#   1. Creating a new repository on GitHub called "kira"
#      (go to github.com -> New repository -> name it "kira"
#       -> do NOT initialize with README -> Create repository)
#   2. Setting up GitHub authentication (see below)
#
# AUTHENTICATION SETUP (do this once):
#   Option A: GitHub CLI (recommended)
#     brew install gh
#     gh auth login
#     (follow prompts, choose HTTPS, authenticate via browser)
#
#   Option B: Personal Access Token
#     1. Go to github.com -> Settings -> Developer settings
#        -> Personal access tokens -> Tokens (classic) -> Generate new token
#     2. Check "repo" scope, generate, copy the token
#     3. When git asks for password, paste the token instead
#
# Then run: bash push_to_github.sh
# ============================================================

set -e

echo "========================================="
echo "  Pushing Kira to GitHub"
echo "========================================="

cd ~/kira

# Check if remote already exists
if git remote get-url origin 2>/dev/null; then
    echo "Remote 'origin' already configured:"
    git remote get-url origin
else
    echo "Adding GitHub remote..."
    # CHANGE THIS to your actual GitHub username if different
    git remote add origin https://github.com/DanielNgabonziza/kira.git
    echo "Remote added: https://github.com/DanielNgabonziza/kira.git"
fi

echo ""
echo "Current branch: $(git branch --show-current)"
echo "Commits: $(git log --oneline | wc -l | tr -d ' ')"
echo ""

# Rename branch to main if it's master (GitHub default is main)
BRANCH=$(git branch --show-current)
if [ "$BRANCH" = "master" ]; then
    echo "Renaming branch master -> main (GitHub convention)..."
    git branch -M main
fi

echo "Pushing to GitHub..."
git push -u origin main

echo ""
echo "========================================="
echo "  Done! Your repo is live at:"
echo "  https://github.com/DanielNgabonziza/kira"
echo "========================================="
