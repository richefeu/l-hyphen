# Makefile for L-HYPHEN Cheatsheet
# Compiles Typst source to PDF

TYPST := typst
SOURCE := lhyphen_cheatsheet.typ
PDF := lhyphen_cheatsheet.pdf

.PHONY: all build clean help watch

all: build

build: $(PDF)

$(PDF): $(SOURCE)
	@echo "Compiling cheatsheet..."
	$(TYPST) compile $(SOURCE) $(PDF)
	@echo "✓ Cheatsheet compiled: $(PDF)"

watch:
	@echo "Watching $(SOURCE) for changes..."
	$(TYPST) watch $(SOURCE) $(PDF)

clean:
	@echo "Cleaning..."
	rm -f $(PDF)
	@echo "✓ Cleaned"

help:
	@echo "L-HYPHEN Cheatsheet Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  all      - Build the PDF (default)"
	@echo "  build    - Compile cheatsheet to PDF"
	@echo "  watch    - Watch source file and recompile on changes"
	@echo "  clean    - Remove generated PDF"
	@echo "  help     - Show this help message"
	@echo ""
	@echo "Requirements:"
	@echo "  - typst (https://typst.app)"
	@echo ""
	@echo "Installation:"
	@echo "  brew install typst  # macOS"
	@echo "  cargo install typst  # or via Rust"
