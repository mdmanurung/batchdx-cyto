"""Tests for post-generation hook functionality."""

from pathlib import Path


def test_post_gen_hook_exists():
    """Test that the post-generation hook file exists.

    This test verifies that the post_gen_project.sh hook file exists
    in the hooks directory.
    """
    hook_file = Path(__file__).parent.parent / "hooks" / "post_gen_project.sh"
    assert hook_file.exists(), "Post-generation hook should exist"


def test_post_gen_hook_contains_env_creation():
    """Test that the post-generation hook contains .env file creation logic.

    This test verifies that the post_gen_project.sh hook contains the
    necessary code to create a .env file.
    """
    hook_file = Path(__file__).parent.parent / "hooks" / "post_gen_project.sh"
    hook_content = hook_file.read_text()

    # Check for .env file creation
    assert "# Create .env file" in hook_content, (
        "Hook should contain .env file creation comment"
    )
    assert ".env" in hook_content, "Hook should reference .env file"
    assert "Environment variables for" in hook_content, (
        "Hook should create environment variables file"
    )


def test_post_gen_hook_contains_git_setup():
    """Test that the post-generation hook contains git setup logic.

    This test verifies that the post_gen_project.sh hook contains the
    necessary code to initialize a git repository.
    """
    hook_file = Path(__file__).parent.parent / "hooks" / "post_gen_project.sh"
    hook_content = hook_file.read_text()

    # Check for git initialization
    assert "git init" in hook_content, "Hook should initialize git repository"
    assert "git add" in hook_content, "Hook should add files to git"
    assert "git commit" in hook_content, "Hook should make initial commit"


def test_post_gen_hook_contains_pixi_setup():
    """Test that the post-generation hook contains pixi setup logic.

    This test verifies that the post_gen_project.sh hook contains the
    necessary code to install pixi environment.
    """
    hook_file = Path(__file__).parent.parent / "hooks" / "post_gen_project.sh"
    hook_content = hook_file.read_text()

    # Check for pixi installation
    assert "pixi install" in hook_content, "Hook should install pixi environment"
    assert "pixi shell" in hook_content, "Hook should mention pixi shell activation"


def test_post_gen_hook_contains_github_setup():
    """Test that the post-generation hook contains GitHub setup logic.

    This test verifies that the post_gen_project.sh hook contains the
    necessary code to set up GitHub repository.
    """
    hook_file = Path(__file__).parent.parent / "hooks" / "post_gen_project.sh"
    hook_content = hook_file.read_text()

    # Check for GitHub CLI usage
    assert "gh repo create" in hook_content, "Hook should create GitHub repository"
    assert "gh auth" in hook_content, "Hook should check GitHub authentication"


def test_post_gen_hook_commits_before_github_push() -> None:
    """Ensure first commit happens before `gh repo create --push`.

    This guards against a regression where the hook attempted to push to GitHub
    before any commit existed, which makes `gh repo create --push` fail.
    """
    hook_file = Path(__file__).parent.parent / "hooks" / "post_gen_project.sh"
    hook_content = hook_file.read_text()

    github_create_idx = hook_content.index("gh repo create")
    github_push_idx = hook_content.index("--push")
    initial_commit_idx = hook_content.index('git commit -m "Initial commit"')

    assert initial_commit_idx < github_create_idx, (
        "Initial commit must happen before creating/pushing GitHub repository"
    )
    assert initial_commit_idx < github_push_idx, (
        "Initial commit must happen before invoking --push"
    )
