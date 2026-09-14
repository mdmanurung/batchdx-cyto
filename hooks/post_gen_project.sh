#!/bin/bash

set -e

echo "🚀 Initializing your new Python project..."

# Cookiecutter runs this hook with cwd already at the generated project root
# (the folder named "{{ cookiecutter.__repo_name }}"). Only cd if we are still
# in the parent directory (older/alternate Cookiecutter behavior).
REPO="{{ cookiecutter.__repo_name }}"
if [ "$(basename "$PWD")" != "$REPO" ] && [ -d "$REPO" ]; then
  cd "$REPO"
fi

# Create .env file
echo "📝 Creating .env file..."
cat > .env << 'EOF'
# Environment variables for {{ cookiecutter.project_name }}
# NOTE: This file is _never_ committed into the git repository!
#       It might contain secrets (e.g. API keys) that should never be exposed publicly.
# export ENV_VAR="some_value"
EOF

# Install pixi environment
echo "📦 Installing pixi environment (this might take a few moments!)..."
if command -v pixi &> /dev/null; then
    pixi install --manifest-path pyproject.toml
else
    echo "❌ pixi not found. Please install pixi first."
    echo "   You can install pixi from: https://pixi.sh"
    exit 1
fi

# Configure Git
echo "🔧 Configuring git..."

# Infer GitHub username from gh CLI
if ! command -v gh &> /dev/null; then
    echo "❌ GitHub CLI not found. Please install it first."
    exit 1
fi

GITHUB_USERNAME=$(gh api user --jq .login 2>/dev/null || echo "")
if [ -z "$GITHUB_USERNAME" ]; then
    echo "❌ Could not infer GitHub username. Please run 'gh auth login' first."
    exit 1
fi

echo "🔍 GitHub account: ${GITHUB_USERNAME}"
gh auth status

# Infer full name and email from git config
FULL_NAME=$(git config user.name 2>/dev/null || echo "")
if [ -z "$FULL_NAME" ]; then
    echo "❌ Could not infer full name. Please set it with: git config --global user.name \"Your Name\""
    exit 1
fi

EMAIL=$(git config user.email 2>/dev/null || echo "")
if [ -z "$EMAIL" ]; then
    echo "❌ Could not infer email. Please set it with: git config --global user.email \"you@example.com\""
    exit 1
fi

echo "📋 Full name: ${FULL_NAME}"
echo "📋 Email: ${EMAIL}"

# Replace placeholders in generated files
echo "🔧 Filling in inferred values..."
sed -i.bak "s|__GITHUB_USERNAME__|${GITHUB_USERNAME}|g" README.md mkdocs.yaml docs/index.md
sed -i.bak "s|__FULL_NAME__|${FULL_NAME}|g" README.md
find . -name "*.bak" -delete

# Initialize git repository
git init -b main

# Prompt user about creating GitHub repository
REPO_NAME="{{ cookiecutter.__repo_name }}"
echo ""

if [ -t 0 ]; then
    echo "🤔 Would you like to create a GitHub repository for this project?"
    echo "   Repository name: ${REPO_NAME}"
    echo "   GitHub username: ${GITHUB_USERNAME}"
    echo ""
    read -p "Create GitHub repository? (y/N): " -n 1 -r
    echo ""
else
    REPLY=""
    echo "ℹ️  Non-interactive mode: skipping GitHub repository creation."
fi

# Install pre-commit hooks
echo "🔧 Installing pre-commit hooks..."
if command -v pre-commit &> /dev/null; then
    pre-commit install
    echo "🔄 Updating pre-commit hooks to latest versions..."
    pre-commit autoupdate
else
    echo "❌ pre-commit not found. Please install it first."
    echo "   You can install it with: uv tool install pre-commit"
    exit 1
fi

# Make initial commit (with safety mechanism for pre-commit hooks)
echo "📝 Making initial commit..."

# Function to attempt git commit with retry logic
commit_with_retry() {
    local max_attempts=3
    local attempt=1

    while [ $attempt -le $max_attempts ]; do
        echo "Attempt $attempt of $max_attempts..."

        if git add . && git commit -m "Initial commit"; then
            echo "✅ Initial commit successful!"
            return 0
        else
            echo ""
            echo "❌ Commit failed on attempt $attempt"
            echo ""
            echo "This is likely due to pre-commit hooks failing (e.g., pydoclint, ruff, etc.)"
            echo "The error output above shows which checks failed."
            echo ""

            if [ $attempt -lt $max_attempts ]; then
                echo "Would you like to:"
                echo "1) Fix the issues and retry"
                echo "2) Skip this commit and handle manually later"
                echo "3) Force commit (not recommended)"
                echo ""
                read -p "Choose option (1/2/3): " -n 1 -r
                echo ""

                case $REPLY in
                    1)
                        echo "Please fix the issues above and press Enter to retry..."
                        read -r
                        ;;
                    2)
                        echo "⏭️  Skipping commit. You can run 'git commit -m \"Initial commit\"' manually later."
                        return 0
                        ;;
                    3)
                        echo "⚠️  Force committing (bypassing pre-commit hooks)..."
                        git commit -m "Initial commit" --no-verify
                        echo "✅ Force commit successful!"
                        return 0
                        ;;
                    *)
                        echo "Invalid option. Retrying..."
                        ;;
                esac
            else
                echo "❌ All attempts failed. Please fix the issues manually and run:"
                echo "   git add ."
                echo "   git commit -m \"Initial commit\""
                return 1
            fi

            attempt=$((attempt + 1))
        fi
    done
}

commit_with_retry

if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "📦 Creating GitHub repository..."
    gh repo create "${GITHUB_USERNAME}/${REPO_NAME}" \
        --private \
        --description "{{ cookiecutter.short_description }}" \
        --source . \
        --remote origin \
        --push
    echo "✅ GitHub repository created and pushed!"
else
    echo "📝 Adding remote origin manually..."
    GIT_SSH_URL="git@github.com:${GITHUB_USERNAME}/${REPO_NAME}"
    git remote add origin "$GIT_SSH_URL"
    echo "✅ Remote origin added. You can create the repository manually on GitHub later."
fi

echo "🎉 Your project has been successfully initialized!"
echo ""
echo "Next steps:"
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "✅ Repository already pushed to GitHub!"
else
    echo "1. Create repository on GitHub and push: git push -u origin main"
fi
echo "2. Activate your environment: pixi shell"
echo "3. Start coding!"