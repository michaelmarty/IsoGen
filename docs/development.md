# Development and testing

## Editable installation

Clone the repository and install the test dependencies:

```shell
python -m pip install -e ".[test]"
```

Run the test suite:

```shell
python -m pytest
```

Pyteomics is intentionally an optional test dependency. The tests use it as an
independent reference for peptide fragment masses and elemental-composition
masses, but the installed IsoGen runtime does not import or require it.

## Documentation

Install the documentation tools:

```shell
python -m pip install -e ".[docs]"
```

Preview changes:

```shell
python -m mkdocs serve
```

Validate the production site:

```shell
python -m mkdocs build --strict
```

### Enable GitHub Pages

Before the documentation workflow can deploy for the first time, a repository
owner must open **Settings → Pages** and set **Build and deployment → Source**
to **GitHub Actions**. GitHub does not allow the workflow's built-in
`GITHUB_TOKEN` to enable Pages automatically.

After that one-time setting is saved, pushes to `main` that change the
documentation or its workflow build and deploy the site automatically. Pull
requests perform the strict documentation build without attempting a Pages
deployment.

## Native libraries

The Python wrapper loads the platform library from the repository-level `bin`
directory during source development and from `isogen/bin` in an installed
wheel. Linux development environments must provide the FFTW 3 runtime.

Pretrained neural-network artifacts are stored in `isogen/models`. The shared
model loader uses that directory by default, and model training saves new
`.pth` and `.bin` artifacts there unless an explicit `working_dir` is supplied.
