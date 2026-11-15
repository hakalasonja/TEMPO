# Contributing to TEMPO

We welcome contributions of all kinds, including bug reports, questions, feature suggestions, documentation improvements, new examples, and code changes.

## Getting help, reporting bugs, and suggesting features

All bug reports, questions, and feature requests are handled through [GitHub Issues](https://github.com/hakalasonja/TEMPO/issues). Before opening a new issue, please search the existing issues to see whether your topic has already been raised.

To help us respond quickly, start your issue title with the type of request:

- **[Bug]** for something that is broken or gives unexpected results
- **[Question]** for help using TEMPO
- **[Feature]** for a suggested new feature or improvement

### Reporting a bug

Please include:

- Your TEMPO, QuTiP, and Python versions, and your operating system. You can print the package versions with:
  ```
  python -c "import tempo, qutip; print(tempo.__version__, qutip.__version__)"
  ```
- A short, self-contained code example that reproduces the problem
- What you were trying to achieve and what happened instead
- The full error message or traceback, if there is one

### Suggesting a feature

Describe the use case you have in mind and what you would like TEMPO to achieve. If you have an idea of how it could be implemented, feel free to include it.

### Private inquiries

For matters that are not suitable for a public issue, contact jjoon@umd.edu.

## Contributing changes

Small changes, such as typo fixes, documentation improvements, or small bug fixes, can be submitted directly as a pull request. For larger changes, such as new features or changes to existing behavior, please open an issue first so we can discuss the approach before you invest time in it.

### Workflow

1. Fork the repository on GitHub and clone your fork.
2. Create a new branch for your changes.
3. Set up a development environment (see below).
4. Make your changes and add tests where relevant.
5. Run the test suite and make sure it passes.
6. Push your branch to your fork and open a pull request against the `main` branch.

In the pull request description, briefly explain what you changed and why. If the pull request addresses an existing issue, reference it (for example, `Closes #12`).

### Development setup

```
git clone https://github.com/<your-username>/TEMPO.git
cd TEMPO
python -m venv venv
source venv/bin/activate
pip install -e ".[dev]"
```

This installs TEMPO in editable mode along with the development tools.

To work on the documentation, also install the documentation dependencies and build the docs locally:

```
pip install -e ".[docs]"
cd docs
make html
```

The built documentation can then be found in `docs/build/html`.

### Running tests

From the root of the repository, run:

```
pytest tempo/tests.py
```

The test suite also runs automatically on every pull request through GitHub Actions.

### Pull request checklist

Before requesting a review, please check that:

- [ ] All tests pass
- [ ] New features and bug fixes include tests in `tempo/tests.py`
- [ ] Docstrings and documentation are updated where relevant
- [ ] Code is formatted with [black](https://black.readthedocs.io/) (`black tempo/`)

A maintainer will review your pull request and may suggest changes before it is merged.

## License

By contributing to TEMPO, you agree that your contributions will be licensed under the [GNU General Public License v3.0 or later](COPYING.txt), the same license that covers the project.