# l-hyphen Documentation (Sphinx)

This directory contains the Sphinx documentation source for l-hyphen.

## Building the Documentation

### Prerequisites

- Python 3.6+
- Sphinx: `pip install sphinx sphinx-rtd-theme`

### Build Commands

```bash
# Build HTML documentation
make html

# Clean and rebuild
make clean html

# View the documentation
open build/html/index.html
```

### Source Files

- **source/index.rst** — Main documentation index
- **source/quickStart.rst** — Quick-start guide
- **source/model.rst** — Mathematical model and theory
- **source/syntax.rst** — Input file syntax reference
- **source/prepro.rst** — Pre-processing with cellPrepro
- **source/images/** — Documentation images

### Output

The compiled HTML documentation is in `build/html/`. Open `index.html` in a web browser to view it.

### Documentation Structure

```
Contents:
├── Quick Start Guide
│   ├── What is l-hyphen?
│   ├── Compilation instructions
│   ├── Running simulations
│   └── Visualizing with see2
├── Model (Theory)
│   ├── DEM overview
│   ├── Cell structure
│   ├── Forces and interactions
│   └── Contact model details
├── Pre-processing
│   ├── Image conversion
│   └── cellPrepro workflow
└── Input File Syntax
    ├── Timing parameters
    ├── Dissipation
    ├── Contact parameters
    ├── Cohesion (glue)
    ├── Cell properties
    ├── Neighbor detection
    └── Output settings
```

## Additional Resources

- **Cheatsheet** (PDF): Quick reference for all configuration commands — `cheatsheets/lhyphen_cheatsheet.pdf`
- **DOCUMENTATION.md** : Comprehensive user guide with examples and troubleshooting
- **GitHub**: https://github.com/richefeu/l-hyphen

## Notes

- Documentation is written in reStructuredText (RST) format
- The theme uses Read the Docs (RTD) Sphinx theme
- Mathematical equations are rendered using LaTeX via MathJax
- French comments in some sections (model.rst, prepro.rst) may be updated for consistency
