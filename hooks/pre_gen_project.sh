#!/bin/bash

set -e

echo "🔧 Checking and setting up required tools..."

# Check if GitHub CLI is installed
if ! command -v gh &> /dev/null; then
    echo "📦 GitHub CLI not found. Installing via pixi..."
    if command -v pixi &> /dev/null; then
        pixi global install gh
    else
        echo "❌ pixi not found. Please install pixi first or install GitHub CLI manually."
        echo "   You can install pixi from: https://pixi.sh/latest/installation/"
        exit 1
    fi
fi

# Check if user is logged into GitHub
if ! gh auth status &> /dev/null; then
    echo "🔐 Not logged into GitHub. Please authenticate..."
    gh auth login
fi

# Check if pre-commit is installed
if ! command -v pre-commit &> /dev/null; then
    echo "📦 pre-commit not found. Installing via uv..."
    if command -v uv &> /dev/null; then
        uv tool install pre-commit
    else
        echo "❌ uv not found. Please install uv first or install pre-commit manually."
        echo "   You can install uv from: https://docs.astral.sh/uv/getting-started/installation/"
        exit 1
    fi
else
    echo "🔄 Upgrading pre-commit..."
    uv tool install --upgrade pre-commit
fi

echo "✅ Pre-generation setup complete!"