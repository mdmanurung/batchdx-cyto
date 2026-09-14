# Cookiecuter Python Project Template

A modern, opinionated Python project template for data science and scientific computing projects.

## Overview

This template creates a fully-featured Python project with modern tooling and best practices built-in. It's designed for data science, scientific computing, and research projects that need robust development workflows.

## Features

### 🏗️ Project Structure
- **Modular design** with separate modules for CLI, models, preprocessing, schemas, and utilities
- **Source layout** with proper package structure
- **Comprehensive testing** setup with pytest and coverage
- **Documentation** with MkDocs and Material theme
- **Marimo notebooks** for reactive exploratory work

### 🛠️ Modern Tooling
- **Pixi** for fast, reliable dependency management
- **Typer** for beautiful CLI interfaces
- **Pre-commit** for code quality hooks
- **Ruff** for formatting, linting, and import sorting
- **Interrogate** for documentation coverage

### 📦 Pre-configured Dependencies
- **Data Science Stack**: pandas, numpy, scikit-learn, matplotlib, seaborn
- **Scientific Computing**: scipy, pymc, jax
- **Development Tools**: typer, python-dotenv, pyprojroot
- **Optional CUDA Support** for GPU acceleration

### 🔧 Development Workflow
- **Automatic Git setup** with GitHub integration
- **Pre-commit hooks** for code quality
- **Environment variable management**
- **Multiple Pixi environments** for different use cases

## Usage

### Prerequisites

Before using this template, ensure you have:

1. **Pixi** installed: https://pixi.sh
2. **GitHub CLI** installed: `pixi global install gh`
3. **GitHub authentication**: `gh auth login`

### Quick Start

```bash
# Install cookiecutter if you haven't already
pip install cookiecutter

# Generate a new project
cookiecutter https://github.com/your-username/cookiecutter-python-project

# Follow the prompts to customize your project
```

### Template Variables

When you run the template, you'll be prompted for:

- **project_name**: The full name of your project (e.g., "My Awesome Analysis")
- **short_description**: A brief description of what your project does

The following are inferred automatically from your environment (no prompts):

- **github_username**: From `gh api user`
- **full_name**: From `git config user.name`
- **email**: From `git config user.email`

The template automatically generates:
- **__project_kebabcase**: URL-friendly version (e.g., "my-awesome-analysis")
- **__project_snakecase**: Python-friendly version (e.g., "my_awesome_analysis")
- **__module_name**: Python module name
- **__package_name**: PyPI package name
- **__repo_name**: GitHub repository name

## Generated Project Structure

```
your-project/
├── src/
│   └── your_project/
│       ├── __init__.py
│       ├── cli.py          # Command-line interface
│       ├── models.py       # Data models and ML models
│       ├── preprocessing.py # Data preprocessing utilities
│       ├── schemas.py      # Data validation schemas
│       └── utils.py        # General utilities
├── tests/
│   ├── __init__.py
│   ├── test_cli.py
│   ├── test_models.py
│   └── ...
├── docs/
│   ├── index.md
│   ├── api.md
│   └── ...
├── pyproject.toml          # Project configuration (Hatchling build + Ruff)
├── README.md              # Project documentation
├── .env                   # Environment variables (gitignored)
└── .gitignore
```

## Development Workflow

### Getting Started

After project generation, the template automatically:

1. **Creates a `.env` file** for environment variables
2. **Installs Pixi environment** with all dependencies
3. **Sets up Git repository** with proper remote
4. **Installs pre-commit hooks** for code quality
5. **Makes initial commit**

### Available Pixi Environments

```bash
# Default environment (tests, devtools, notebook, setup)
pixi shell

# Documentation environment
pixi shell --feature docs

# Testing environment
pixi shell --feature tests

# CUDA environment (for GPU acceleration)
pixi shell --feature cuda
```

### Common Commands

```bash
# Run tests
pixi run test

# Format code
pixi run lint

# Build documentation
pixi run build-docs

# Serve documentation locally
pixi run serve-docs

# Update pre-commit hooks
pixi run update
```

### CLI Development

The template includes a basic CLI structure using Typer:

```python
# src/your_project/cli.py
import typer

app = typer.Typer()

@app.command()
def hello():
    """Echo the project's name."""
    typer.echo("This project's name is Your Project")

@app.command()
def describe():
    """Describe the project."""
    typer.echo("Your project description")
```

Your CLI will be available as `your-project-name` after installation.

## Configuration

### Customizing Dependencies

Edit `pyproject.toml` to add or remove dependencies:

```toml
[project]
dependencies = [
    "pandas",
    "numpy",
    # Add your dependencies here
]

[tool.pixi.dependencies]
# Add conda-forge packages here
```

### Environment Variables

Use the `.env` file for configuration:

```bash
# .env
export API_KEY="your-secret-key"
export DATABASE_URL="your-database-url"
```

### Pre-commit Hooks

The template includes pre-commit hooks for:
- Formatting and linting (Ruff, including import sorting via the `I` rules)
- Documentation coverage (Interrogate)

## Publishing

### To PyPI

```bash
# Build and publish
python -m build
python -m twine upload dist/*
```

### Documentation

Documentation is automatically built and can be deployed to GitHub Pages:

```bash
# Build docs
pixi run build-docs

# Deploy to GitHub Pages (if configured)
mkdocs gh-deploy
```

## Customization

### Adding New Modules

1. Create new files in `src/your_project/`
2. Add imports to `__init__.py`
3. Update tests in `tests/`

### Modifying Templates

The template uses Jinja2 templating. Key variables:
- `{{ cookiecutter.project_name }}`: Full project name
- `{{ cookiecutter.__module_name }}`: Python module name
- `{{ cookiecutter.__package_name }}`: Package name

### Hooks

The template includes pre and post-generation hooks:

- **pre_gen_project.sh**: Checks for required tools (GitHub CLI, pre-commit)
- **post_gen_project.sh**: Sets up Git, Pixi environment, and initial commit

## Best Practices

### Code Organization

- Keep models in `models.py`
- Put data preprocessing in `preprocessing.py`
- Define schemas in `schemas.py`
- General utilities go in `utils.py`
- CLI commands in `cli.py`

### Testing

- Write tests for all public functions
- Use pytest fixtures for common setup
- Aim for high test coverage
- Use hypothesis for property-based testing

### Documentation

- Document all public functions and classes
- Keep README.md up to date
- Use type hints for better IDE support
- Write docstrings in Google or NumPy style

## Troubleshooting

### Common Issues

**Pixi not found**: Install from https://pixi.sh

**GitHub CLI not authenticated**: Run `gh auth login`

**Pre-commit hooks failing**: Run `pixi run update` to update hooks

**CUDA not working**: Ensure you have CUDA 12 installed and use `pixi shell --feature cuda`

### Getting Help

- Check the [Pixi documentation](https://pixi.sh)
- Review [Typer documentation](https://typer.tiangolo.com)
- See [MkDocs Material documentation](https://squidfunk.github.io/mkdocs-material/)

## Contributing

To contribute to this template:

1. Fork the repository
2. Make your changes
3. Test with a sample project
4. Submit a pull request

## License

This template is licensed under the MIT License.